"""Actual Qt widgets, worker updates, groups, point masks and PNG/SVG export."""
from dataclasses import replace
import os
import json
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
from pathlib import Path
from threading import Event
from unittest.mock import patch
import tempfile,time,unittest
import numpy as np
from PySide6 import QtCore,QtGui,QtWidgets,QtTest
import pyqtgraph as pg
from ui import Window, GROUP_COLUMN, INTEGRAL_COLUMN, ILT_COLUMN
from spectrum_settings import SpectrumAppearance, SpectrumSettingsDialog
from core import default_options
from test_dosy_lab import STE,synthetic
from generate_demo import demo_bytes
from storage import DATA, PROJECTS, PROCESSING, INTEGRALS, PROJECT_FILTER, PROFILE_FILTER, INTEGRALS_FILTER
from jeol import RawData
from core import process_fids

APP=QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
# The Windows offscreen Qt plugin has no system font database by default.
if not QtGui.QFontDatabase.families():
    for name in ('segoeui.ttf','segoeuib.ttf','arial.ttf'):
        font=Path('C:/Windows/Fonts')/name
        if font.exists():QtGui.QFontDatabase.addApplicationFont(str(font))
APP.setFont(QtGui.QFont('Segoe UI',10))
pg.setConfigOptions(background='#101827',foreground='#dbe4ef')


def settle(window,seconds=30):
    end=time.monotonic()+seconds
    while time.monotonic()<end:
        APP.processEvents();QtTest.QTest.qWait(15)
        if window.worker is None and not window.timer.isActive() and not window.pending:
            if not window.valid:raise AssertionError(window.statusBar().currentMessage())
            return
    raise AssertionError('Worker did not finish: '+window.statusBar().currentMessage())


def prepared_window(settings=None,source=STE):
    w=Window(settings);w.show();w.open_path(source);settle(w)
    options=default_options(w.raw);options.update(lb_Hz=.3,zero_fill_factor=2)
    if source==STE:
        options.update(phase1_deg=7200.,phase0_deg=42.29972276350504)
    w.processing.set_options(options,len(w.raw.fid));w.request_calculation();settle(w)
    return w


def drag(viewport,start,end,button=QtCore.Qt.MouseButton.LeftButton):
    QtTest.QTest.mousePress(viewport,button,pos=start)
    for t in (.25,.5,.75,1.):
        pos=QtCore.QPoint(round(start.x()+t*(end.x()-start.x())),round(start.y()+t*(end.y()-start.y())))
        QtTest.QTest.mouseMove(viewport,pos,15)
    QtTest.QTest.mouseRelease(viewport,button,pos=end)
    APP.processEvents()


class InterfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = STE
        if not STE.is_file():
            directory = cls.enterClassContext(tempfile.TemporaryDirectory())
            cls.source = Path(directory)/'DOSY_demo.jdf'
            cls.source.write_bytes(demo_bytes())

    def setUp(self):
        directory = self.enterContext(tempfile.TemporaryDirectory())
        self.settings_path = str(Path(directory) / 'preferences.ini')
        settings = QtCore.QSettings(self.settings_path, QtCore.QSettings.Format.IniFormat)
        self.w = prepared_window(settings,self.source)
    def tearDown(self):
        if self.w.worker is not None:settle(self.w)
        self.w.close();APP.processEvents()

    def test_phase_navigation_shortcuts_and_double_click_reset(self):
        self.w.show_phase()
        dialog = self.w.phase_window
        dialog.activateWindow()
        dialog.plot.setFocus()
        APP.processEvents()
        self.assertTrue(dialog.zoom_button.text().endswith('[Z]'))
        self.assertTrue(dialog.pan_button.text().endswith('[P]'))
        self.assertTrue(any(b.text().endswith('[F]') for b in dialog.findChildren(QtWidgets.QPushButton)))
        QtTest.QTest.keyClicks(dialog.plot, 'z')
        self.assertTrue(dialog.zoom_button.isChecked())
        self.assertFalse(dialog.view.phasing)
        QtTest.QTest.keyClicks(dialog.plot, 'p')
        self.assertTrue(dialog.pan_button.isChecked())
        self.assertFalse(dialog.zoom_button.isChecked())
        self.assertEqual(dialog.view.mode, 'pan')
        QtTest.QTest.keyClicks(dialog.plot, 'p')
        self.assertTrue(dialog.view.phasing)
        dialog.reset_view()
        expected = dialog.view.viewRange()[0]
        dialog.plot.setXRange(3.6, 3.7, padding=0.)
        QtTest.QTest.keyClicks(dialog.plot, 'f')
        np.testing.assert_allclose(dialog.view.viewRange()[0], expected)
        dialog.plot.setXRange(3.6, 3.7, padding=0.)
        view = dialog.plot.viewport()
        QtTest.QTest.mouseDClick(view, QtCore.Qt.MouseButton.LeftButton, pos=view.rect().center())
        np.testing.assert_allclose(dialog.view.viewRange()[0], expected)
        dialog.reject()

    def test_shift_scales_only_new_mouse_movement_for_both_phases(self):
        self.w.show_phase()
        APP.processEvents()
        dialog = self.w.phase_window
        viewport = dialog.plot.viewport()
        start = QtCore.QPoint(round(viewport.width() * .3), round(viewport.height() * .55))
        shift = QtCore.Qt.KeyboardModifier.ShiftModifier
        normal = QtCore.Qt.KeyboardModifier.NoModifier
        for button, key, speed, fine_speed in (
            (QtCore.Qt.MouseButton.LeftButton, 'phase0_deg', .5, .025),
            (QtCore.Qt.MouseButton.RightButton, 'phase1_deg', 10., .25),
        ):
            with self.subTest(key=key):
                control = dialog.controls[key]
                before = control.value()
                QtTest.QTest.mousePress(viewport, button, pos=start)
                for offset, modifiers, expected in (
                    (20, normal, 20 * speed),
                    (40, shift, 20 * speed + 20 * fine_speed),
                    (60, normal, 40 * speed + 20 * fine_speed),
                ):
                    pos = start + QtCore.QPoint(0, -offset)
                    QtTest.QTest.qWait(20)
                    event = QtGui.QMouseEvent(
                        QtCore.QEvent.Type.MouseMove, QtCore.QPointF(pos),
                        QtCore.QPointF(viewport.mapToGlobal(pos)),
                        QtCore.Qt.MouseButton.NoButton, button, modifiers)
                    APP.sendEvent(viewport, event)
                    self.assertAlmostEqual(control.value() - before, expected, places=6)
                QtTest.QTest.mouseRelease(viewport, button, pos=pos)
        dialog.reject()

    def test_spectrum_settings_update_views_and_persist(self):
        w = self.w
        w.add_roi([3.60, 3.68])
        settle(w)
        dataset = w.dataset
        revision = w.revision
        old_range = w.view.viewRange()

        def edit(dialog):
            with patch.object(QtWidgets.QColorDialog, 'getColor', return_value=QtGui.QColor('#ab45ef')):
                dialog.color_button.click()
            dialog.width.setValue(2.4)
            dialog.grid.setChecked(False)
            return QtWidgets.QDialog.DialogCode.Accepted

        with patch.object(SpectrumSettingsDialog, 'exec', edit):
            w.edit_spectrum_settings()
        self.assertIs(w.dataset, dataset)
        self.assertEqual(w.revision, revision)
        np.testing.assert_allclose(w.view.viewRange(), old_range)
        w.show_phase()
        APP.processEvents()
        for plot in (w.spectrum, w.roi_overlay, w.phase_window.plot):
            pen = plot.listDataItems()[0].opts['pen']
            self.assertEqual(pen.color().name(), '#ab45ef')
            self.assertAlmostEqual(pen.widthF(), 2.4)
            self.assertFalse(plot.getAxis('bottom').grid)
            self.assertFalse(plot.getAxis('left').grid)
        w.phase_window.reject()
        settings = QtCore.QSettings(self.settings_path, QtCore.QSettings.Format.IniFormat)
        self.assertEqual(SpectrumAppearance.load(settings), SpectrumAppearance('#ab45ef', 2.4, False))

    @unittest.skipUnless(STE.is_file(),'Lokalny JDF sacharozy jest opcjonalny.')
    def test_fresh_JDF_load_applies_filter_ramp_before_manual_phasing(self):
        w=self.w;w.open_path(STE);settle(w)
        self.assertEqual(w.processing.controls['filter_delay_samples'].value(),19.6875)
        self.assertEqual(w.dataset.options['phase1_deg'],0)
        w.info_button.click()
        self.assertIn('model JDF',w.info_window.processing.toPlainText())
        self.assertIn('Skany na warstwę',w.info_window.experiment.toPlainText())
        self.assertIn('relaxation_delay',w.info_window.parameters.toPlainText())
        w.info_window.close()
        x=w.dataset.ppm;mask=(x>=5.65)&(x<=6.35)
        self.assertLess(np.std(w.dataset.spectrum[0,mask].real),1.)

    def test_mouse_ROI_then_type_and_reuse_group(self):
        w=self.w;view=w.spectrum.getViewBox();viewport=w.spectrum.viewport()
        w.roi_button.click()
        a=w.spectrum.mapFromScene(view.mapViewToScene(QtCore.QPointF(3.60,0)))
        b=w.spectrum.mapFromScene(view.mapViewToScene(QtCore.QPointF(3.68,0)))
        QtTest.QTest.mousePress(viewport,QtCore.Qt.MouseButton.LeftButton,pos=a)
        QtTest.QTest.mouseMove(viewport,b,60);QtTest.QTest.mouseRelease(viewport,QtCore.Qt.MouseButton.LeftButton,pos=b)
        settle(w)
        self.assertEqual(len(w.entries),1)
        self.assertTrue(w.results[0]['fit']['ok'])
        table=w.roi_table
        table.editItem(table.item(0,GROUP_COLUMN))
        editor=table.findChild(QtWidgets.QComboBox)
        self.assertTrue(editor.isEditable())
        QtTest.QTest.keyClicks(editor.lineEdit(),'Sacharoza')
        QtTest.QTest.keyClick(editor.lineEdit(),QtCore.Qt.Key.Key_Return)
        settle(w)
        self.assertEqual(w.entries[0]['group'],'Sacharoza')
        w.region_items[w.active_id].setRegion((3.605,3.675));settle(w)
        np.testing.assert_allclose(w.entries[0]['bounds'],[3.605,3.675])
        w.add_roi([5.26,5.33]);settle(w)
        table.editItem(table.item(1,GROUP_COLUMN))
        editor=table.findChild(QtWidgets.QComboBox)
        self.assertEqual([editor.itemText(i) for i in range(editor.count())],['Sacharoza'])
        editor.setCurrentIndex(0)
        editor.activated.emit(0)
        settle(w)
        self.assertEqual(w.entries[1]['group'],'Sacharoza')
        self.assertEqual([table.item(i,0).text() for i in range(2)],['S1','S2'])

    def test_modify_and_split_integral_with_live_preview(self):
        w=self.w
        w.add_roi([3.60,3.68],name='H-test',group='Sacharoza')
        settle(w)
        w.spectrum.setXRange(3.5,3.9,padding=0)
        y_low,y_high=w.view.viewRange()[1]
        w.spectrum.setYRange(min(y_low,-.2*(y_high-y_low)),y_high,padding=0)
        view=w.view
        viewport=w.spectrum.viewport()
        y=sum(view.viewRange()[1])/2
        def point(x):
            return w.spectrum.mapFromScene(view.mapViewToScene(QtCore.QPointF(x,y)))
        w.mode_buttons['modify'].click()
        drag(viewport,point(3.64),point(3.68))
        settle(w)
        parent=w.entries[0]
        lo,hi=parent['bounds']
        np.testing.assert_allclose([lo,hi],[3.64,3.72],atol=.002)
        self.assertLess(w.integral_labels[parent['id']].pos().y(),0)
        measured=w.results[0]['result']['corrected_integral'].copy()
        w.mode_buttons['split'].click()
        cut=(lo+hi)/2
        start,end=point(hi),point(cut)
        QtTest.QTest.mousePress(viewport,QtCore.Qt.MouseButton.LeftButton,pos=start)
        try:
            QtTest.QTest.mouseMove(viewport,end,70)
            APP.processEvents()
            self.assertEqual(w.entries,[parent])
            self.assertTrue(view.selection_preview.isVisible())
            self.assertEqual(w.region_items[parent['id']].getRegion(),(lo,hi))
        finally:
            QtTest.QTest.mouseRelease(viewport,QtCore.Qt.MouseButton.LeftButton,pos=end)
        settle(w)
        self.assertFalse(view.selection_preview.isVisible())
        first,second=w.entries
        self.assertEqual(first['bounds'][1],second['bounds'][0])
        self.assertAlmostEqual(first['bounds'][1],cut,delta=.002)
        self.assertEqual((first['bounds'][0],second['bounds'][1]),(lo,hi))
        self.assertEqual([first['group'],second['group']],['Sacharoza','Sacharoza'])
        self.assertEqual(len({first['id'],second['id'],parent['id']}),3)
        total=sum(result['result']['corrected_integral'] for result in w.results)
        np.testing.assert_allclose(total,measured,atol=1e-10)

    def test_clear_baseline_restores_spectrum_after_subtraction(self):
        w=self.w
        before=w.dataset.spectrum.copy()
        w.processing.add_region([2.65,2.95])
        w.processing.add_region([5.65,6.1])
        w.processing.controls['baseline_enabled'].setChecked(True)
        settle(w)
        self.assertTrue(w.clear_baseline_button.isEnabled())
        self.assertFalse(np.array_equal(before,w.dataset.spectrum))
        w.mode_buttons['baseline'].click()
        w.clear_baseline_button.click()
        settle(w)
        self.assertEqual(w.processing.regions,[])
        self.assertEqual(w.base_items,[])
        self.assertFalse(w.dataset.options['baseline_enabled'])
        self.assertIsNone(w.dataset.diagnostics['baseline'])
        self.assertEqual(w.view.mode,'zoom')
        np.testing.assert_array_equal(w.dataset.spectrum,before)

    def test_dialogs_start_in_workspace_folders_and_filter_extensions(self):
        w=self.w
        stem=self.source.stem
        with patch.object(QtWidgets.QFileDialog,'getOpenFileName',return_value=('','')) as dialog:
            w.choose_file()
            self.assertEqual(dialog.call_args.args[2:],(str(DATA),'JEOL (*.jdf)'))
            w.choose_project()
            self.assertEqual(dialog.call_args.args[2:],(str(PROJECTS),PROJECT_FILTER))
            w.processing.load()
            self.assertEqual(dialog.call_args.args[2:],(str(PROCESSING),PROFILE_FILTER))
        with patch('ui.choose_save_file',return_value=None) as dialog:
            w.save_project()
            self.assertEqual(dialog.call_args.args[2:],(PROJECTS,f'{stem}.proj',PROJECT_FILTER))
        with patch('processing_panel.choose_save_file',return_value=None) as dialog:
            w.processing.save()
            self.assertEqual(dialog.call_args.args[2:],(PROCESSING,f'{stem}.list',PROFILE_FILTER))
        w.show_map()
        w.show_statistics()
        s=w.statistics_window
        with patch.object(QtWidgets.QFileDialog,'getSaveFileName',return_value=('','')) as dialog:
            w.choose_spectrum_png()
            self.assertEqual(dialog.call_args.args[2],f'{stem}.png')
            w.map_window.choose_save()
            self.assertEqual(dialog.call_args.args[2],f'{stem}.png')
            s.choose_plot_export('png')
            self.assertEqual(dialog.call_args.args[2],f'{stem}.png')

            s.pin_current()
            saved_id=next(key for key in s.experiments if key != s.live_id)
            s.experiments[saved_id]['source_name']='Pomiar.02.JDF'
            s.experiment_combo.setCurrentIndex(s.experiment_combo.findData(saved_id))
            s.choose_plot_export('svg')
            self.assertEqual(dialog.call_args.args[2],'Pomiar.02.svg')

    def test_integral_template_preserves_definitions_and_recalculates_current_measurement(self):
        w = self.w
        w.add_roi([3.602031234567, 3.680912345678], group='Związek α')
        w.next_id = 7
        w.add_roi([5.26, 5.33], group='Druga grupa')
        settle(w)
        w.entries[0]['fit_mask'][0] = False
        w.integral_scale = 3.
        expected = [dict(id=e['id'], group=e['group'], bounds=e['bounds']) for e in w.entries]
        path = Path(self.settings_path).with_name('ranges.ints')
        with patch('ui.choose_save_file', return_value=str(path)) as dialog:
            w.export_integrals()
            self.assertEqual(dialog.call_args.args[2:], (INTEGRALS, f'{self.source.stem}.ints', INTEGRALS_FILTER))
        packet = json.loads(path.read_text(encoding='utf-8'))
        self.assertEqual(packet['integrals'], expected)

        w.open_path(self.source)
        settle(w)
        w.add_roi([3.50, 3.56], group='Inny')
        settle(w)
        with patch.object(QtWidgets.QFileDialog, 'getOpenFileName', return_value=(str(path), '')) as dialog:
            with patch.object(QtWidgets.QMessageBox, 'question', return_value=QtWidgets.QMessageBox.StandardButton.No):
                w.import_integrals()
            self.assertEqual(dialog.call_args.args[2:], (str(INTEGRALS), INTEGRALS_FILTER))
            self.assertEqual(w.entries[0]['group'], 'Inny')

            malformed = dict(packet, integrals=[expected[0], expected[0]])
            path.write_text(json.dumps(malformed), encoding='utf-8')
            with patch.object(QtWidgets.QMessageBox, 'warning') as warning:
                w.import_integrals()
            self.assertIn('Powtórzone ID', warning.call_args.args[2])
            self.assertEqual(w.entries[0]['group'], 'Inny')

            path.write_text(json.dumps(packet), encoding='utf-8')
            with patch.object(QtWidgets.QMessageBox, 'question', return_value=QtWidgets.QMessageBox.StandardButton.Yes):
                w.import_integrals()
        settle(w)
        self.assertEqual([dict(id=e['id'], group=e['group'], bounds=e['bounds']) for e in w.entries], expected)
        self.assertTrue(all(all(e['fit_mask']) and not e['weighted'] for e in w.entries))
        self.assertEqual(w.integral_scale, 1.)
        self.assertEqual(w.next_id, 8)
        self.assertTrue(all(len(e['result']['corrected_integral']) == len(w.raw.fid) for e in w.results))

    def test_switch_language_preserves_project_and_relabels_open_windows(self):
        w=self.w
        w.add_roi([3.60,3.68],group='Przesuwaj')
        settle(w)
        w.show_info();w.show_statistics();w.show_map()
        APP.processEvents()
        w.tabs.setCurrentIndex(1)
        w.spectrum.setXRange(3.6,3.7,padding=0)
        original=(w.dataset,w.results,w.view.viewRange())
        w.language_button.click()
        APP.processEvents()
        self.assertEqual(w.open_button.text(),'Open DOSY')
        self.assertEqual([w.tabs.tabText(i) for i in range(2)],['Decays','Residuals'])
        self.assertTrue(w.residual.isVisible())
        self.assertFalse(w.decay.isVisible())
        self.assertEqual(w.roi_table.horizontalHeaderItem(1).text(),'Group / compound')
        self.assertEqual(w.roi_table.item(0,1).text(),'Przesuwaj')
        self.assertIn('.',w.roi_table.item(0,INTEGRAL_COLUMN).text())
        self.assertIn('Source of b',w.info_window.experiment.toPlainText())
        self.assertEqual(w.statistics_window.group_combo.currentText(),'Przesuwaj')
        self.assertEqual(w.statistics_window.export_kind.itemText(0),'D summary')
        self.assertIn('2D DOSY map',w.map_window.windowTitle())
        self.assertIs(w.dataset,original[0])
        self.assertIs(w.results,original[1])
        np.testing.assert_array_equal(w.view.viewRange(),original[2])
        self.assertEqual(w.settings.value('language'),'en')
        w.language_button.click()
        self.assertEqual(w.open_button.text(),'Otwórz DOSY')
        self.assertEqual(w.tabs.tabText(1),'Reszty')
        self.assertIn(',',w.roi_table.item(0,INTEGRAL_COLUMN).text())
        self.assertEqual(w.settings.value('language'),'pl')

    def test_load_selects_lowest_b_and_preview_layer_keeps_common_phase(self):
        w=self.w;r=synthetic();order=np.r_[3,0,2,1,np.arange(4,16)]
        raw=RawData(r.fid[order].copy(),r.gradients[order].copy(),r.b[order].copy(),dict(r.metadata))
        o=default_options(raw);o.update(phase0_deg=24.,phase1_deg=135.,pivot_ppm=4.3)
        d=process_fids(raw,o);w.loaded((raw,d,[],[],1.))
        self.assertEqual(w.layer.currentIndex(),1);self.assertEqual(w.processing.options()['reference_layer'],1)
        before=w.dataset.spectrum.copy();w.layer.setCurrentIndex(7);APP.processEvents()
        self.assertEqual(w.processing.options()['phase0_deg'],24.)
        self.assertEqual(w.processing.options()['phase1_deg'],135.)
        self.assertEqual(w.processing.options()['pivot_ppm'],4.3)
        np.testing.assert_array_equal(w.dataset.spectrum,before)

    def test_phase_processing_and_point_exclusion_refresh_results(self):
        w=self.w;raw=w.raw.fid.copy();w.add_roi([3.60,3.68],group='Sacharoza');settle(w)
        before=w.results[0]['fit']['D_m2_s'];w.point_table.item(0,0).setCheckState(QtCore.Qt.CheckState.Unchecked);settle(w)
        self.assertFalse(w.results[0]['fit']['fit_mask'][0]);self.assertEqual(w.results[0]['fit']['n_active'],15)
        self.assertEqual(w.results[0]['previous_D'],before)
        w.processing.controls['phase0_deg'].setValue(43);settle(w)
        self.assertEqual(w.dataset.options['phase0_deg'],43)
        w.processing.controls['phase1_deg'].setValue(7210)
        w.processing.controls['pivot_ppm'].setValue(4.51);settle(w)
        expected=process_fids(w.raw,w.processing.options())
        np.testing.assert_array_equal(w.dataset.phased,expected.phased)
        np.testing.assert_array_equal(w.raw.fid,raw)

    def test_ilt_button_uses_cached_signed_data_and_exports_plots(self):
        w = self.w
        w.add_roi([3.60, 3.68], group='Sacharoza')
        settle(w)
        w.point_table.item(0, 0).setCheckState(QtCore.Qt.CheckState.Unchecked)
        settle(w)
        dataset, results = w.dataset, w.results
        fit = results[0]['fit']

        def wait_for(predicate):
            end = time.monotonic() + 10
            while not predicate() and time.monotonic() < end:
                QtTest.QTest.qWait(15)
            self.assertTrue(predicate(), dialog.status.text())

        self.assertIsNone(w.ilt_window)
        with patch('ui.analyze', side_effect=AssertionError('ILT must not redo the mono fit')):
            w.roi_table.cellWidget(0, ILT_COLUMN).click()
            dialog = w.ilt_window
            self.assertTrue(dialog.isVisible())
            self.check_ilt_help(dialog, 'ILT pojedynczej całki', 'λ', 'Druga pochodna')
            wait_for(lambda: dialog.result is not None)
            np.testing.assert_array_equal(dialog.b, dataset.b)
            np.testing.assert_array_equal(dialog.y, results[0]['result']['corrected_integral'])
            np.testing.assert_array_equal(dialog.mask, fit['fit_mask'])
            self.assertFalse(dialog.mask[0])
            self.assertIs(dialog.fit, fit)
            dialog.slider.setValue(-63)
            dialog.slider.setValue(-55)
            dialog.lambda_edit.setText('0.000024137')
            dialog.lambda_changed()
            dialog.points.setCurrentIndex(dialog.points.findData(50))
            wait_for(lambda: dialog.solved_parameters.regularization == .000024137)
            self.assertEqual(dialog.solved_parameters.points, 50)
            self.assertTrue(dialog.distribution.getAxis('bottom').logMode)
            result = dialog.result
            dialog.log_axis.setChecked(False)
            self.assertFalse(dialog.distribution.getAxis('bottom').logMode)
            self.assertIs(dialog.result, result)
            dialog.log_axis.setChecked(True)
            dialog.stability.setChecked(True)
            wait_for(dialog.stability_summary.isVisible)
            self.assertIn('90%', dialog.stability_summary.text())
            self.assertIn('λ=', dialog.stability_summary.text())
            np.testing.assert_array_equal(dialog.result.weights,result.weights)
            self.assertIs(dialog.fit,fit)
            dialog.distribution.fullscreen_controller.open()
            self.assertIsNotNone(dialog.distribution.fullscreen_controller.window)
            dialog.distribution.fullscreen_controller.close()
            APP.processEvents()

            with tempfile.TemporaryDirectory() as folder:
                single = Path(folder)/'distribution.png'
                combined = Path(folder)/'ilt.png'
                with patch.object(QtWidgets.QFileDialog, 'getSaveFileName', side_effect=[(str(single), ''), (str(combined), '')]):
                    dialog.export_png((dialog.distribution,), 'distribution')
                    dialog.combined_export.click()
                self.assertGreaterEqual(QtGui.QImage(str(single)).width(), 1600)
                self.assertGreater(QtGui.QImage(str(combined)).height(), QtGui.QImage(str(single)).height())
            dialog.close()
            self.assertIsNone(w.ilt_window)
            self.assertIs(w.dataset, dataset)
            self.assertIs(w.results, results)
            self.assertTrue(w.valid)

    def test_global_ilt_selected_rows_comparison_and_export(self):
        w = self.w
        for bounds in ([3.60, 3.68], [5.26, 5.33], [3.74, 3.82]):
            w.add_roi(bounds)
        settle(w)
        table = w.roi_table
        for row, modifier in [(0, QtCore.Qt.KeyboardModifier.NoModifier), (1, QtCore.Qt.KeyboardModifier.ControlModifier)]:
            pos = table.visualItemRect(table.item(row, 0)).center()
            QtTest.QTest.mouseClick(table.viewport(), QtCore.Qt.MouseButton.LeftButton, modifier, pos)
        self.assertEqual(w.selected_roi_ids(), ['S1', 'S2'])
        w.render_roi_table()
        self.assertEqual(w.selected_roi_ids(), ['S1', 'S2'])
        self.assertEqual(w.global_ilt_button.text(), 'Global ILT (2)')
        results = w.results

        def wait_for(predicate):
            end = time.monotonic()+10
            while not predicate() and time.monotonic() < end:
                QtTest.QTest.qWait(15)
            self.assertTrue(predicate(), dialog.status.text())

        with patch('ui.analyze', side_effect=AssertionError('Global ILT must not redo the mono fit')):
            w.global_ilt_button.click()
            dialog = w.global_ilt_window
            self.assertIsNotNone(dialog)
            self.assertEqual(dialog.ids, ['S1', 'S2'])
            self.check_ilt_help(dialog, 'Global ILT dla wielu całek', 'λ — regularyzacja Tichonowa', 'μ — wspólna rzadkość')
            wait_for(lambda: dialog.result is not None)
            np.testing.assert_array_equal(dialog.b, w.dataset.b)
            np.testing.assert_array_equal(dialog.y[:,1], results[1]['result']['corrected_integral'])
            self.assertIs(dialog.fits[0], results[0]['fit'])
            self.assertEqual(dialog.table.rowCount(), 2)
            self.assertTrue(dialog.result.components)
            relative = dialog.heatmap_values()
            np.testing.assert_allclose(relative.sum(axis=0), [100, 100])
            dialog.absolute.setChecked(True)
            np.testing.assert_array_equal(dialog.heatmap_values(), dialog.result.contributions)
            dialog.visible_rois[0].setChecked(False)
            self.assertEqual([item.name() for item in dialog.roi_plot.listDataItems()], ['S2'])
            dialog.log_axis.setChecked(False)
            self.assertFalse(dialog.roi_plot.getAxis('bottom').logMode)
            dialog.penalties[1][1].setText('0')
            dialog.compare.setChecked(True)
            dialog.points.setCurrentIndex(dialog.points.findData(50))
            dialog.recalculate()
            wait_for(lambda: len(dialog.independent) == 2)
            self.assertEqual(dialog.solved_parameters.sparsity, 0.)
            self.assertEqual(dialog.result.weights.shape, (50, 2))
            self.assertTrue(all(dialog.independent))
            dialog.roi.setCurrentIndex(1)
            self.assertTrue(dialog.result.converged)
            self.assertIs(w.map_ilt[1], dialog.result)
            w.show_map()
            panel = w.map_window.alignment
            self.assertTrue(panel.mode.model().item(3).isEnabled())
            panel.mode.setCurrentIndex(3)
            wait_for(lambda: panel.ready)
            self.assertEqual(set(panel.profiles) | set(panel.rejected), {'S1','S2'})
            with tempfile.TemporaryDirectory() as folder:
                single, combined = [Path(folder)/name for name in ('heatmap.png', 'global.png')]
                with (patch.object(QtWidgets.QFileDialog, 'getSaveFileName', side_effect=[(str(single), ''), (str(combined), '')]),
                      patch.object(QtWidgets.QMessageBox, 'warning', side_effect=AssertionError('PNG export failed'))):
                    dialog.export_png((dialog.heatmap,), 'heatmap')
                    dialog.export_png(dialog.plots, 'combined')
                image = QtGui.QImage(str(single))
                self.assertEqual(image.width(), 1600)
                self.assertEqual(image.dotsPerMeterX(), 11811)
                montage = QtGui.QImage(str(combined))
                self.assertGreater(montage.height(), image.height())
                for exported in (image, montage):
                    rgb = exported.convertToFormat(QtGui.QImage.Format.Format_RGB32)
                    pixels = np.frombuffer(rgb.constBits(), dtype=np.uint32)
                    self.assertGreater(np.count_nonzero(pixels == 0xff60a5fa), 20)
            dialog.close()
            self.assertIsNone(w.global_ilt_window)
            self.assertIs(w.results, results)
            self.assertTrue(w.valid)

    def test_global_ilt_cancels_an_obsolete_running_job(self):
        from global_ilt import SolveCancelled, solve_global_ilt

        w = self.w
        w.add_roi([3.60, 3.68]); w.add_roi([5.26, 5.33])
        settle(w)
        w.roi_table.selectAll()
        started, cancelled = Event(), Event()
        seen = []

        def solve(b, y, parameters, **kwargs):
            seen.append(parameters.sparsity)
            if len(seen) == 1:
                started.set()
                deadline = time.monotonic()+5
                while not kwargs['cancelled']() and time.monotonic() < deadline:
                    time.sleep(.01)
                if not kwargs['cancelled']():
                    raise RuntimeError('Expected cancellation')
                cancelled.set()
                raise SolveCancelled('Obsolete')
            return solve_global_ilt(b, y, parameters, **kwargs)

        with patch('global_ilt_window.solve_global_ilt', side_effect=solve):
            w.global_ilt_button.click()
            dialog = w.global_ilt_window
            deadline = time.monotonic()+8
            while not started.is_set() and time.monotonic() < deadline:
                QtTest.QTest.qWait(15)
            self.assertTrue(started.is_set())
            dialog.penalties[1][0].setValue(-55)
            dialog.penalties[1][0].setValue(-60)
            while dialog.result is None and time.monotonic() < deadline:
                QtTest.QTest.qWait(15)
            self.assertTrue(cancelled.is_set())
            self.assertIsNotNone(dialog.result, dialog.status.text())
            self.assertEqual(seen, [1e-4, 1e-6])
            self.assertEqual(dialog.solved_parameters.sparsity, 1e-6)
            dialog.close()

    def check_ilt_help(self, dialog, *phrases):
        self.assertTrue(dialog.windowFlags() & QtCore.Qt.WindowType.WindowCloseButtonHint)
        self.assertTrue(dialog.windowFlags() & QtCore.Qt.WindowType.WindowSystemMenuHint)
        info = dialog.findChild(QtWidgets.QPushButton, 'ilt_info_button')
        info.click()
        guide = dialog.findChild(QtWidgets.QDialog, 'ilt_info_window')
        self.assertTrue(guide.isVisible())
        text = guide.findChild(QtWidgets.QTextBrowser).toPlainText()
        for phrase in phrases:
            self.assertIn(phrase, text)
        info.click()
        self.assertEqual(len(dialog.findChildren(QtWidgets.QDialog, 'ilt_info_window')), 1)
        guide.close()
        self.assertTrue(dialog.isVisible())

    def test_ilt_discards_running_solution_when_parameters_change(self):
        from ilt import solve_ilt

        w = self.w
        w.add_roi([3.60, 3.68])
        settle(w)
        started, release = Event(), Event()
        solved_parameters = []

        def solve(b, intensity, parameters, **kwargs):
            solved_parameters.append(parameters)
            started.set()
            if not release.wait(10):
                raise RuntimeError('Test did not release worker')
            return solve_ilt(b, intensity, parameters, **kwargs)

        with patch('ilt_window.solve_ilt', side_effect=solve):
            w.roi_table.cellWidget(0, ILT_COLUMN).click()
            dialog = w.ilt_window
            try:
                end = time.monotonic()+8
                while not started.is_set() and time.monotonic() < end:
                    QtTest.QTest.qWait(15)
                self.assertTrue(started.is_set())
                dialog.slider.setValue(-60)
                dialog.slider.setValue(-70)
                self.assertFalse(dialog.combined_export.isEnabled())
            finally:
                release.set()
            end = time.monotonic()+8
            while dialog.result is None and time.monotonic() < end:
                QtTest.QTest.qWait(15)
            self.assertIsNotNone(dialog.result)
            self.assertEqual([p.regularization for p in solved_parameters], [1e-5, 1e-7])
            self.assertEqual(dialog.solved_parameters.regularization, 1e-7)
            dialog.close()

    def test_reference_click_assigns_ppm_and_preserves_phase_and_integrals(self):
        w = self.w
        panel = w.processing
        panel.show()
        panel.controls['ppm_shift'].setValue(.25)
        settle(w)
        w.add_roi([3.85, 3.93])
        panel.add_region([2.90, 3.20])
        panel.add_region([5.90, 6.35])
        panel.controls['baseline_enabled'].setChecked(True)
        settle(w)
        before = w.dataset
        integrals = w.results[0]['result']['corrected_integral'].copy()
        w.spectrum.setXRange(3.8, 4.0, padding=0.)
        APP.processEvents()
        x_range = np.array(w.view.viewRange()[0])
        y = np.mean(w.view.viewRange()[1])
        pos = w.spectrum.mapFromScene(w.view.mapViewToScene(QtCore.QPointF(3.89, y)))
        clicked = w.view.mapSceneToView(w.spectrum.mapToScene(pos)).x()
        viewport = w.spectrum.viewport()
        self.assertIn('Docelowa wartość referencji', [label.text() for label in panel.findChildren(QtWidgets.QLabel)])
        self.assertFalse(w.reference_guide.isVisible())

        panel.reference_target.setValue(1.25)
        self.assertTrue(w.valid)
        panel.reference_button.click()
        self.assertEqual(w.view.mode, 'reference')
        for point in (pos + QtCore.QPoint(30, 0), pos):
            QtTest.QTest.mouseMove(viewport, point, 35)
            APP.processEvents()
            self.assertTrue(w.reference_guide.isVisible())
            ppm = w.view.mapSceneToView(w.spectrum.mapToScene(point)).x()
            self.assertAlmostEqual(w.reference_guide.value(), ppm)
        np.testing.assert_array_equal(w.view.viewRange()[0], x_range)
        APP.sendEvent(viewport, QtCore.QEvent(QtCore.QEvent.Type.Leave))
        self.assertFalse(w.reference_guide.isVisible())
        QtTest.QTest.mouseClick(viewport, QtCore.Qt.MouseButton.LeftButton, pos=pos)
        settle(w)
        self.assertEqual(w.view.mode, 'zoom')
        self.assertFalse(panel.reference_button.isChecked())
        self.assertFalse(w.reference_guide.isVisible())
        expected_shift = round(.25 + 1.25 - clicked, 6)
        self.assertEqual(panel.controls['ppm_shift'].value(), expected_shift)
        delta = expected_shift - .25
        self.assertAlmostEqual(clicked + delta, 1.25, places=6)
        np.testing.assert_allclose(w.dataset.ppm, before.ppm + delta, atol=1e-12)
        np.testing.assert_allclose(w.view.viewRange()[0], x_range + delta, atol=1e-9)
        np.testing.assert_allclose(w.dataset.spectrum, before.spectrum, atol=1e-9)
        np.testing.assert_allclose(w.entries[0]['bounds'], np.array([3.85, 3.93]) + delta)
        np.testing.assert_allclose(w.dataset.options['baseline_regions'], np.array(before.options['baseline_regions']) + delta)
        np.testing.assert_allclose(w.results[0]['result']['corrected_integral'], integrals, atol=1e-9)
        w.language_button.click()
        self.assertEqual(panel.reference_button.text(), 'Reference from spectrum')
        self.assertIn('Target reference value', [label.text() for label in panel.findChildren(QtWidgets.QLabel)])

    def test_baseline_preview_apply_and_undo(self):
        w=self.w;w.add_roi([3.60,3.68],group='Sacharoza');settle(w);original=w.dataset.spectrum.copy()
        self.assertEqual(w.processing.controls['baseline_degree'].value(),1)
        per_layer=w.processing.controls['baseline_per_layer']
        self.assertTrue(per_layer.isChecked())
        w.processing.add_region([2.65,2.95]);w.processing.add_region([5.65,6.1]);settle(w)
        self.assertIsNotNone(w.dataset.diagnostics['baseline']);np.testing.assert_array_equal(w.dataset.spectrum,original)
        w.overlay_layers.setText('1,3');w.overlay.setChecked(True)
        previews=[curve.getData()[1] for curve in w.spectrum.listDataItems()
                  if curve.opts['pen'].color().name()=='#fbbf24']
        np.testing.assert_allclose(previews,w.dataset.baseline[[0,2]])
        per_layer.setChecked(False);settle(w)
        np.testing.assert_array_equal(w.dataset.baseline[0],w.dataset.baseline[2])
        w.processing.controls['baseline_enabled'].setChecked(True);settle(w)
        self.assertTrue(per_layer.isChecked())
        self.assertFalse(per_layer.isEnabled())
        self.assertFalse(np.array_equal(w.dataset.baseline[0],w.dataset.baseline[2]))
        self.assertFalse(np.array_equal(w.dataset.spectrum,original))
        w.processing.undo_baseline();settle(w)
        self.assertTrue(per_layer.isEnabled())
        np.testing.assert_array_equal(w.dataset.spectrum,original)

    def test_statistics_groups_selection_and_diagnostics_never_refit(self):
        w=self.w;w.add_roi([3.60,3.68],name='H1',group='A');w.add_roi([5.26,5.33],name='H2',group='B');settle(w)
        w.show_statistics();s=w.statistics_window;APP.processEvents()
        self.assertEqual(s.summary['N'],0)
        self.assertEqual(s.scatter.points().size,0)
        s.group_combo.setCurrentIndex(s.group_combo.findData('A'));self.assertEqual(s.summary['N'],1)
        self.assertEqual(s.d_ticks(),[(0,'S1')])
        self.assertEqual([p.data() for p in s.scatter.points()],['S1'])
        with patch('core.fit_decay',side_effect=AssertionError('Unexpected refit')):
            s.select_roi('S1');s.table.item(0,0).setCheckState(QtCore.Qt.CheckState.Unchecked)
            self.assertEqual(s.summary['N'],0)
            self.assertEqual(s.scatter.points().size,0)
            self.assertEqual(s.highlight.points().size,0)
            s.group_combo.setCurrentIndex(s.group_combo.findData('B'))
            self.assertEqual(s.d_ticks(),[(0,'S2')])
            self.assertEqual([p.data() for p in s.scatter.points()],['S2'])
            s.group_combo.setCurrentIndex(s.group_combo.findData('A'))
            s.table.item(0,0).setCheckState(QtCore.Qt.CheckState.Checked)
            self.assertEqual(s.summary['N'],1);s.show_point(2)
        w.processing.controls['lb_Hz'].setValue(.4);settle(w)
        self.assertEqual(s.current['statistics_group'],'A');self.assertEqual(s.summary['N'],1)
        self.assertEqual(s.current['processing']['lb_Hz'],.4)

    def test_light_spectrum_and_statistics_exports(self):
        w=self.w
        w.add_roi([3.60,3.68],group='A')
        w.add_roi([5.26,5.33],group='A')
        settle(w)
        w.spectrum.setXRange(3.50,3.85,padding=0)
        before=np.array(w.view.viewRange())
        live_pen=w.spectrum.listDataItems()[0].opts['pen'].color()
        w.show_statistics()
        s=w.statistics_window
        APP.processEvents()
        s.d_grid.click()
        self.assertFalse(s.d_plot.getAxis('left').grid)
        self.assertTrue(s.decay_plot.getAxis('left').grid)
        s.detail_grid.click()
        self.assertFalse(s.decay_plot.getAxis('left').grid)
        self.assertFalse(s.residual_plot.getAxis('left').grid)
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory)
            png=folder/'spectrum.png'
            with patch.object(QtWidgets.QFileDialog,'getSaveFileName',return_value=(str(png),'')):
                w.save_spectrum_button.click()
            for kind in range(3):
                s.export_plot(folder/f'stats_{kind}.png',kind)
                svg=folder/f'stats_{kind}.svg'
                s.export_plot(svg,kind)
                markup=svg.read_text(encoding='utf-8')
                self.assertIn('fill="#ffffff"',markup)
                self.assertIn('#142b5c',markup)
            for path in folder.glob('*.png'):
                image=QtGui.QImage(str(path)).convertToFormat(QtGui.QImage.Format.Format_RGB888)
                self.assertFalse(image.isNull(), path)
                rows=np.frombuffer(image.constBits(),dtype=np.uint8).reshape(image.height(),image.bytesPerLine())
                colors=rows[:,:3*image.width()].reshape(-1,3)
                self.assertGreater(np.mean(np.all(colors>=210,axis=1)),.90,path)
                self.assertTrue(np.any(np.all(colors==255,axis=1)),path)
                self.assertTrue(np.any(np.all(colors==0,axis=1)),path)
                self.assertTrue(np.any(np.all(colors==[20,43,92],axis=1)),path)
            self.assertTrue(png.is_file())
        self.assertEqual(w.spectrum.listDataItems()[0].opts['pen'].color(),live_pen)
        np.testing.assert_allclose(w.view.viewRange(),before)

    def test_statistics_plots_zoom_reset_and_open_fullscreen_by_button(self):
        w = self.w
        w.add_roi([3.60,3.68],group='Sacharoza')
        w.add_roi([5.26,5.33],group='Sacharoza')
        settle(w)
        w.show_statistics()
        s = w.statistics_window
        APP.processEvents()
        self.assertGreater(s.group_combo.mapTo(s,QtCore.QPoint()).y(),
                           s.experiment_combo.mapTo(s,QtCore.QPoint()).y())
        for plot,button,reset in (
            (s.d_plot,s.d_fullscreen,s.reset_overview),
            (s.decay_plot,s.decay_fullscreen,s.reset_diagnostics),
            (s.residual_plot,s.residual_fullscreen,s.reset_diagnostics),
        ):
            with self.subTest(plot=plot):
                reset()
                APP.processEvents()
                view = plot.getViewBox()
                before = np.array(view.viewRange())
                span = np.ptp(before,axis=1)
                start = plot.mapFromScene(view.mapViewToScene(QtCore.QPointF(*(before[:,0]+.2*span))))
                end = plot.mapFromScene(view.mapViewToScene(QtCore.QPointF(*(before[:,0]+.7*span))))
                viewport = plot.viewport()
                QtTest.QTest.mousePress(viewport,QtCore.Qt.MouseButton.LeftButton,pos=start)
                try:
                    QtTest.QTest.mouseMove(viewport,end,70)
                    APP.processEvents()
                    self.assertTrue(view.rbScaleBox.isVisible())
                finally:
                    QtTest.QTest.mouseRelease(viewport,QtCore.Qt.MouseButton.LeftButton,pos=end)
                self.assertTrue((np.ptp(view.viewRange(),axis=1)<.8*span).all())
                QtTest.QTest.mouseDClick(viewport,QtCore.Qt.MouseButton.LeftButton,pos=viewport.rect().center())
                self.assertIsNone(plot.fullscreen_controller.window)
                restored = np.array(view.viewRange())
                self.assertTrue((np.max(abs(restored-before),axis=1)<.06*span).all())
                parent = plot.parentWidget()
                button.click()
                APP.processEvents()
                self.assertIsNotNone(plot.fullscreen_controller.window)
                plot.fullscreen_controller.close()
                APP.processEvents()
                self.assertIs(plot.parentWidget(),parent)

    def test_export_fullscreen_and_restore(self):
        w=self.w;w.add_roi([3.60,3.68],group='Sacharoza');settle(w)
        controller=w.decay.fullscreen_controller;parent=w.decay.parentWidget()
        button=next(b for b in w.findChildren(QtWidgets.QPushButton) if b.text()=='Pełny ekran zaników')
        button.click();APP.processEvents();self.assertIsNotNone(controller.window)
        controller.close();APP.processEvents();self.assertEqual(w.decay.parentWidget(),parent)
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)/'export';w.export_to(folder);settle(w)
            for name in ['DOSY_results.xlsx','integrals_and_fits.csv','spectrum.png','spectrum.svg','decays.svg','Raport.md']:
                self.assertTrue((folder/name).exists(),name)
                self.assertGreater((folder/name).stat().st_size,50)

    def test_unverified_b_still_shows_measured_points(self):
        w=self.w
        w.raw=replace(w.raw,b=None)
        w.add_roi([3.60,3.68])
        settle(w)
        self.assertFalse(w.results[0]['fit']['ok'])
        scatter=next(item for item in w.decay.plotItem.items if isinstance(item,pg.ScatterPlotItem))
        self.assertEqual(len(scatter.points()),len(w.raw.fid))
        np.testing.assert_allclose(scatter.getData()[0],w.raw.gradients)
        self.assertIn('brak zweryfikowanego b',w.decay.getAxis('bottom').labelText)
        self.assertIn('b=brak',w.point_tooltip('S1',0))

    def test_dosy_alignment_modes_preserve_cached_fits_and_allow_manual_membership(self):
        w = self.w
        w.add_roi([3.60, 3.68]); w.add_roi([5.26, 5.33])
        settle(w)
        original = w.results
        d_original = [row['fit']['D_m2_s'] for row in original]
        with patch('ui.analyze', side_effect=AssertionError('Map modes must preserve independent fits')):
            w.show_map()
            window = w.map_window
            panel = window.alignment

            def wait_for_alignment():
                until = time.monotonic()+5
                while not panel.ready and not panel.error and time.monotonic() < until:
                    QtTest.QTest.qWait(15)
                self.assertTrue(panel.ready, panel.error)

            self.assertEqual(panel.mode.currentIndex(), 0)
            self.assertFalse(panel.mode.model().item(3).isEnabled())
            panel.tolerance.setValue(20.)
            panel.mode.setCurrentIndex(1)
            wait_for_alignment()
            overrides, _, _ = panel.map_data()
            self.assertEqual(set(overrides), {'S1','S2'})
            self.assertEqual(overrides['S1'], overrides['S2'])
            group = panel.tree.topLevelItem(0)
            group.setCheckState(0, QtCore.Qt.CheckState.Unchecked)
            self.assertEqual(len(window.images), 0)
            group.setCheckState(0, QtCore.Qt.CheckState.Checked)
            group.child(0).setCheckState(0, QtCore.Qt.CheckState.Unchecked)
            wait_for_alignment()
            self.assertFalse(panel.map_data()[0])
            outside = panel.tree.topLevelItem(0)
            excluded = next(outside.child(i) for i in range(outside.childCount())
                            if outside.child(i).checkState(0) == QtCore.Qt.CheckState.Unchecked)
            excluded.setCheckState(0, QtCore.Qt.CheckState.Checked)
            panel.mode.setCurrentIndex(2)
            wait_for_alignment()
            self.assertEqual(len(panel.groups[0].diagnostics), 2)
            self.assertIn('Wspólny fit D', window.note.text())
            panel.mode.setCurrentIndex(0)
            self.assertFalse(panel.map_data()[0])
            self.assertIs(w.results, original)
            self.assertEqual([row['fit']['D_m2_s'] for row in w.results], d_original)

    def test_diffusion_projection_links_d_and_delete_clears_all_guides(self):
        w = self.w
        w.add_roi([3.60, 3.68])
        w.add_roi([5.26, 5.33])
        settle(w)
        w.show_map()
        window = w.map_window
        self.assertFalse(window.diffusion_plot.isVisible())
        window.diffusion_button.click()
        APP.processEvents()
        self.assertTrue(window.diffusion_plot.isVisible())
        self.assertTrue(window.settings.value('dosy_map/diffusion_projection', False, type=bool))
        original, d = window.diffusion_trace.getData()
        self.assertGreater(np.max(original), 0.)
        self.assertEqual(len(original), len(d))
        before_x, before_d = window.diffusion_plot.getViewBox().viewRange()
        viewport = window.diffusion_plot.viewport()
        center = viewport.rect().center()
        for delta, factor in [(120, 1/1.2), (-120, 1.)]:
            event = QtGui.QWheelEvent(QtCore.QPointF(center), QtCore.QPointF(viewport.mapToGlobal(center)),
                                      QtCore.QPoint(), QtCore.QPoint(0, delta), QtCore.Qt.MouseButton.NoButton,
                                      QtCore.Qt.KeyboardModifier.NoModifier, QtCore.Qt.ScrollPhase.NoScrollPhase, False)
            APP.sendEvent(viewport, event)
            APP.processEvents()
            x_range, d_range = window.diffusion_plot.getViewBox().viewRange()
            self.assertAlmostEqual(np.ptp(x_range), factor*np.ptp(before_x), places=7)
            np.testing.assert_allclose(d_range, before_d)
        window.floor.setValue(25.)
        window.gain.setValue(10)
        np.testing.assert_array_equal(window.diffusion_trace.getData()[0], original)
        window.floor.setValue(0.)
        window.rois.item(1).setCheckState(QtCore.Qt.CheckState.Unchecked)
        remaining, remaining_d = window.diffusion_trace.getData()
        self.assertLess(remaining.sum()*np.diff(remaining_d).mean(), original.sum()*np.diff(d).mean())

        low, high = window.full_range[1]
        window.plot.setYRange(low, (low+high)/2, padding=0.)
        APP.processEvents()
        y = (3*low+high)/4
        window.show_coordinates(3.64, y)
        guide = window.diffusion_projection_guide
        self.assertTrue(guide.isVisible())
        self.assertEqual(guide.value(), window.diffusion_guide.value())

        def screen_y(plot):
            scene = plot.getViewBox().mapViewToScene(QtCore.QPointF(0., y))
            return plot.viewport().mapToGlobal(plot.mapFromScene(scene)).y()

        def screen_x(plot):
            scene = plot.getViewBox().mapViewToScene(QtCore.QPointF(3.64, 0.))
            return plot.viewport().mapToGlobal(plot.mapFromScene(scene)).x()

        self.assertLessEqual(abs(screen_y(window.diffusion_plot)-screen_y(window.plot)), 1)
        self.assertLessEqual(abs(screen_x(window.projection)-screen_x(window.plot)), 1)
        window.activateWindow()
        window.plot.setFocus()
        APP.processEvents()
        QtTest.QTest.keyClicks(window.plot, 'p')
        self.assertTrue(window.tools.button(1).isChecked())
        QtTest.QTest.keyClicks(window.plot, 'z')
        self.assertTrue(window.tools.button(0).isChecked())
        QtTest.QTest.keyClicks(window.plot, 'f')
        np.testing.assert_allclose(window.plot.getViewBox().viewRange()[1], window.full_range[1])
        with tempfile.TemporaryDirectory() as folder:
            marked, cleared = [Path(folder)/name for name in ('marked.png', 'cleared.png')]
            self.assertTrue(window.save_view(marked))
            QtTest.QTest.keyClick(window.plot.viewport(), QtCore.Qt.Key.Key_Delete)
            for item in (guide, window.diffusion_guide, window.ppm_guide, window.projection_guide, window.d_label):
                self.assertFalse(item.isVisible())
            self.assertTrue(window.save_view(cleared))
            first, second = [QtGui.QImage(str(path)).convertToFormat(QtGui.QImage.Format.Format_RGB32) for path in (marked, cleared)]
            self.assertEqual(first.size(), second.size())
            width = window.diffusion_plot.width()
            a, b = (np.frombuffer(image.constBits(), dtype=np.uint32).reshape(image.height(), -1) for image in (first, second))
            self.assertGreater(np.count_nonzero(a[:, :width] != b[:, :width]), 20)
        window.diffusion_button.click()
        self.assertFalse(window.diffusion_plot.isVisible())

    def test_adaptive_local_ilt_uses_cached_rows_and_updates_mode(self):
        from dosy_alignment_panel import ADAPTIVE

        w = self.w
        for bounds in ([3.60,3.68], [5.26,5.33], [4.66,4.74]):
            w.add_roi(bounds)
        settle(w)
        original = w.results
        with patch('ui.analyze', side_effect=AssertionError('Adaptive map must preserve independent fits')):
            w.show_map()
            window = w.map_window
            panel = window.alignment
            self.assertFalse(panel.adaptive_controls.isVisible())
            panel.tolerance.setValue(20.)
            panel.mode.setCurrentIndex(ADAPTIVE)
            self.assertTrue(panel.adaptive_controls.isVisible())

            def wait():
                end = time.monotonic()+8
                while not panel.ready and not panel.error and time.monotonic() < end:
                    QtTest.QTest.qWait(15)
                self.assertTrue(panel.ready, panel.error)

            wait()
            self.assertIn(('S1','S2'), [group.ids for group in panel.local_groups])
            group = panel.tree.topLevelItem(0)
            panel.tree.setCurrentItem(group)
            self.assertEqual(panel.highlighted, {'S1','S2'})
            self.assertTrue(any(image.opacity() < 1. for image in window.images))
            group.setCheckState(0, QtCore.Qt.CheckState.Unchecked)
            self.assertLess(len(window.images), 3)
            panel.balance.setValue(3.)
            panel.regularization.setText('2e-6')
            panel.recalculate()
            wait()
            self.assertIs(w.results, original)
            panel.reset_defaults()
            wait()
            self.assertEqual(panel.balance.value(), 10.)
            self.assertEqual(panel.regularization.text(), '1e-6')

    def test_dosy_map_navigation_coordinates_floor_and_export(self):
        w = self.w
        w.add_roi([3.60, 3.68], group='Sacharoza')
        w.add_roi([5.26, 5.33], group='Sacharoza')
        settle(w)
        with patch('ui.analyze', side_effect=AssertionError('Mapa nie może ponawiać fitu.')):
            w.map_button.click()
            APP.processEvents()
            window = w.map_window
            self.assertTrue(window.isVisible())
            self.assertEqual(len(window.images), 2)
            group = window.groups.topLevelItem(0)
            self.assertEqual(group.text(0), 'Sacharoza')
            window.rois.item(1).setCheckState(QtCore.Qt.CheckState.Unchecked)
            group.setCheckState(0, QtCore.Qt.CheckState.Unchecked)
            self.assertEqual(len(window.images), 0)
            group.setCheckState(0, QtCore.Qt.CheckState.Checked)
            self.assertEqual(len(window.images), 1)
            self.assertEqual(window.rois.item(1).checkState(), QtCore.Qt.CheckState.Unchecked)
            window.rois.item(1).setCheckState(QtCore.Qt.CheckState.Checked)
            self.assertEqual(len(window.images), 2)
            layer = w.dataset.options['reference_layer']
            np.testing.assert_array_equal(window.trace.getData()[1], w.dataset.spectrum[layer].real)
            plot, view = window.plot, window.map_view
            original = np.array(view.viewRange())
            def point(x, y):
                pos = original[:,0] + np.array([x, y])*np.ptp(original, axis=1)
                return plot.mapFromScene(view.mapViewToScene(QtCore.QPointF(*pos)))
            drag(plot.viewport(), point(.2,.25), point(.75,.75))
            self.assertLess(np.ptp(view.viewRange()[0]), .7*np.ptp(original[0]))
            np.testing.assert_allclose(window.reference_view.viewRange()[0], view.viewRange()[0], atol=.005)
            QtTest.QTest.mouseDClick(plot.viewport(), QtCore.Qt.MouseButton.LeftButton, pos=plot.viewport().rect().center())
            np.testing.assert_allclose(view.viewRange(), original, atol=.005)
            window.tools.button(1).click()
            drag(plot.viewport(), point(.4,.4), point(.58,.45))
            self.assertGreater(abs(view.viewRange()[0][0] - original[0,0]), .1)
            window.reset_view()
            window.tools.button(2).click()
            QtTest.QTest.mouseClick(plot.viewport(), QtCore.Qt.MouseButton.LeftButton, pos=point(.5,.5))
            self.assertTrue(window.ppm_guide.isVisible())
            self.assertAlmostEqual(window.ppm_guide.value(), original[0].mean(), delta=.02)
            self.assertAlmostEqual(window.diffusion_guide.value(), original[1].mean(), delta=.02)
            self.assertIn('m²/s', window.coordinate.text())
            self.assertTrue(window.d_label.isVisible())
            self.assertIn(f'{window.diffusion_guide.value()*window.D_SCALE:.3e}',window.d_label.toPlainText())
            self.assertEqual(window.d_label.pos(),QtCore.QPointF(window.ppm_guide.value(),window.diffusion_guide.value()))
            self.assertTrue(view.sceneBoundingRect().contains(window.d_label.sceneBoundingRect().center()))
            def wheel(widget):
                viewport = widget.viewport()
                pos = viewport.rect().center()
                APP.sendEvent(viewport, QtGui.QWheelEvent(
                    QtCore.QPointF(pos), QtCore.QPointF(viewport.mapToGlobal(pos)),
                    QtCore.QPoint(), QtCore.QPoint(0,120), QtCore.Qt.MouseButton.NoButton,
                    QtCore.Qt.KeyboardModifier.NoModifier, QtCore.Qt.ScrollPhase.ScrollUpdate, False))
            floor = window.floor.value()
            wheel(plot)
            self.assertGreater(window.floor.value(), floor)
            np.testing.assert_allclose(view.viewRange(), original)
            height = np.ptp(window.reference_view.viewRange()[1])
            wheel(window.projection)
            self.assertLess(np.ptp(window.reference_view.viewRange()[1]), height)
            image = window.images[0]
            peak = window.peak_intensity
            image.setImage(np.array([[-peak, 0., .02*peak, peak]], dtype=np.float32), autoLevels=False)
            window.floor.setValue(0)
            image.render()
            weak_alpha = image.qimage.pixelColor(2,0).alpha()
            self.assertEqual(image.qimage.pixelColor(3,0).name(), '#fb923c')
            window.gain.setValue(20)
            image.render()
            self.assertEqual(window.gain_label.text(), '10×')
            self.assertEqual(image.qimage.pixelColor(0,0).alpha(), 0)
            self.assertEqual(image.qimage.pixelColor(1,0).alpha(), 0)
            self.assertGreater(image.qimage.pixelColor(2,0).alpha(), weak_alpha)
            window.floor.setValue(3)
            image.render()
            self.assertEqual(image.qimage.pixelColor(2,0).alpha(), 0)
            with patch.object(QtWidgets.QColorDialog, 'getColor',
                              side_effect=[QtGui.QColor('#22c55e'), QtGui.QColor('#a78bfa')]):
                window.map_color_action.trigger()
                window.trace_color_action.trigger()
            image.render()
            self.assertEqual(image.qimage.pixelColor(3,0).name(), '#22c55e')
            self.assertEqual(window.trace.opts['pen'].color().name(), '#a78bfa')
            window.grid_button.click()
            for chart in (window.plot, window.projection):
                self.assertFalse(chart.getAxis('bottom').grid)
                self.assertFalse(chart.getAxis('left').grid)
            self.assertEqual(w.settings.value('dosy_map/color'), '#22c55e')
            self.assertEqual(w.settings.value('dosy_map/trace_color'), '#a78bfa')
            self.assertFalse(w.settings.value('dosy_map/grid', type=bool))
            window.rois.item(1).setCheckState(QtCore.Qt.CheckState.Unchecked)
            self.assertEqual(len(window.images), 1)
            with tempfile.TemporaryDirectory() as directory:
                png = Path(directory)/'map.png'
                self.assertTrue(window.save_view(png))
                saved = QtGui.QImage(str(png)).convertToFormat(QtGui.QImage.Format.Format_RGB888)
                self.assertFalse(saved.isNull())
                window.d_label.hide()
                without_label=Path(directory)/'without_label.png'
                self.assertTrue(window.save_view(without_label))
                window.d_label.show()
                plain=QtGui.QImage(str(without_label)).convertToFormat(QtGui.QImage.Format.Format_RGB888)
                self.assertGreater(np.count_nonzero(np.frombuffer(saved.constBits(),dtype=np.uint8)
                                                    !=np.frombuffer(plain.constBits(),dtype=np.uint8)),30)
                rows = np.frombuffer(saved.constBits(), dtype=np.uint8).reshape(saved.height(), saved.bytesPerLine())
                pixels = rows[:, :3*saved.width()].reshape(saved.height(), saved.width(), 3)
                self.assertGreater(np.mean(np.all(pixels == 255, axis=2)), .80)
                self.assertTrue(np.any(np.all(pixels == 0, axis=2)))
                for panel in (pixels[:window.projection.height()], pixels[window.projection.height()+7:]):
                    self.assertTrue(np.any(np.all(panel == [20, 43, 92], axis=2)))
                self.assertEqual(window.trace.opts['pen'].color().name(), '#a78bfa')
                np.testing.assert_array_equal(window.images[0].lut[-1,:3], [34,197,94])
                np.testing.assert_allclose(view.viewRange(), original)

    def test_calibration_rescales_cached_integrals_and_round_trips_project(self):
        w = self.w
        w.add_roi([3.60,3.68])
        w.add_roi([5.26,5.33])
        settle(w)
        measured = [result['result']['corrected_integral'].copy() for result in w.results]
        revision = w.revision
        fits = w.results
        w.roi_table.item(0, INTEGRAL_COLUMN).setText('1,00')
        self.assertEqual(w.roi_table.item(0, INTEGRAL_COLUMN).text(), '1,00')
        expected = f'{measured[1][0]/measured[0][0]:.2f}'.replace('.', ',')
        self.assertEqual(w.roi_table.item(1, INTEGRAL_COLUMN).text(), expected)
        self.assertEqual(w.integral_labels['S2'].textItem.toPlainText(), expected)
        self.assertEqual(w.signal_labels['S2'].textItem.toPlainText(), 'S2')
        self.assertGreater(w.signal_labels['S2'].pos().y(), w.integral_labels['S2'].pos().y())
        self.assertEqual(w.revision, revision)
        self.assertIs(w.results, fits)
        w.roi_table.item(0, INTEGRAL_COLUMN).setText('0')
        self.assertEqual(w.roi_table.item(0, INTEGRAL_COLUMN).text(), '1,00')
        w.layer.setCurrentIndex(3)
        expected = f'{measured[1][3]/measured[0][0]:.2f}'.replace('.', ',')
        self.assertEqual(w.roi_table.item(1, INTEGRAL_COLUMN).text(), expected)
        self.assertEqual(w.integral_labels['S2'].textItem.toPlainText(), expected)
        path = str(Path(self.settings_path).with_name('project.proj'))
        with patch('ui.choose_save_file', return_value=path):
            w.save_project()
        w.open_path(path, project=True)
        settle(w)
        with patch('ui.choose_save_file', return_value=None) as dialog:
            w.save_project()
            self.assertEqual(dialog.call_args.args[3], f'{self.source.stem}.proj')
        self.assertEqual(w.roi_table.item(0, INTEGRAL_COLUMN).text(), '1,00')
        self.assertEqual(w.integral_labels['S1'].textItem.toPlainText(), '1,00')
        self.assertEqual(w.signal_labels['S1'].textItem.toPlainText(), 'S1')
        for result, original in zip(w.results, measured):
            np.testing.assert_array_equal(result['result']['corrected_integral'], original)

    def test_decay_rectangle_preview_and_double_click_resets_views(self):
        w = self.w
        w.add_roi([3.60,3.68])
        settle(w)
        w.spectrum.setXRange(3.6,3.7,padding=0)
        viewport = w.spectrum.viewport()
        QtTest.QTest.mouseDClick(viewport,QtCore.Qt.MouseButton.LeftButton,pos=viewport.rect().center())
        self.assertGreater(np.ptp(w.view.viewRange()[0]),1.)
        self.assertIsNone(w.spectrum.fullscreen_controller.window)
        for plot in (w.decay, w.residual):
            with self.subTest(plot=plot):
                w.tabs.setCurrentWidget(plot.parentWidget())
                w.reset_decay_views()
                APP.processEvents()
                view = plot.getViewBox()
                before = np.array(view.viewRange())
                start = plot.mapFromScene(view.mapViewToScene(QtCore.QPointF(*(before[:,0]+.2*np.ptp(before,axis=1)))))
                end = plot.mapFromScene(view.mapViewToScene(QtCore.QPointF(*(before[:,0]+.7*np.ptp(before,axis=1)))))
                viewport = plot.viewport()
                QtTest.QTest.mousePress(viewport,QtCore.Qt.MouseButton.LeftButton,pos=start)
                try:
                    QtTest.QTest.mouseMove(viewport,end,70)
                    APP.processEvents()
                    self.assertTrue(view.rbScaleBox.isVisible())
                    np.testing.assert_allclose(view.viewRange(),before)
                finally:
                    QtTest.QTest.mouseRelease(viewport,QtCore.Qt.MouseButton.LeftButton,pos=end)
                self.assertFalse(view.rbScaleBox.isVisible())
                self.assertTrue((np.ptp(view.viewRange(),axis=1)<.8*np.ptp(before,axis=1)).all())
                QtTest.QTest.mouseDClick(viewport,QtCore.Qt.MouseButton.LeftButton,pos=viewport.rect().center())
                self.assertIsNone(plot.fullscreen_controller.window)
                reset = np.array(view.viewRange())
                # Pixel-sized markers can shift auto-range padding slightly after a zoom.
                self.assertTrue((np.max(abs(reset-before),axis=1)<.03*np.ptp(before,axis=1)).all())

    def test_wheel_changes_display_height_only_and_fullscreen_button(self):
        w=self.w;view=w.view;viewport=w.spectrum.viewport();pos=viewport.rect().center()
        original=w.dataset.spectrum.copy();before=np.array(view.viewRange())
        event=QtGui.QWheelEvent(QtCore.QPointF(pos),QtCore.QPointF(viewport.mapToGlobal(pos)),
            QtCore.QPoint(),QtCore.QPoint(0,120),QtCore.Qt.MouseButton.NoButton,
            QtCore.Qt.KeyboardModifier.NoModifier,QtCore.Qt.ScrollPhase.NoScrollPhase,False)
        with patch('core.fit_decay',side_effect=AssertionError('Display change refitted data')):
            APP.sendEvent(viewport,event);APP.processEvents()
            after=np.array(view.viewRange())
            np.testing.assert_allclose(after[0],before[0])
            self.assertAlmostEqual(np.ptp(after[1])/np.ptp(before[1]),1/1.2,places=7)
            np.testing.assert_array_equal(w.dataset.spectrum,original)
            parent=w.spectrum.parentWidget();w.fullscreen_button.click();APP.processEvents()
            self.assertIsNotNone(w.spectrum.fullscreen_controller.window)
            w.spectrum.fullscreen_controller.close();APP.processEvents()
            self.assertEqual(w.spectrum.parentWidget(),parent)

    def test_default_drag_zooms_and_ROI_is_explicit_toggle(self):
        w=self.w;view=w.view;viewport=w.spectrum.viewport()
        self.assertEqual(view.mode,'zoom');self.assertFalse(w.roi_button.isChecked())
        def point(x):return w.spectrum.mapFromScene(view.mapViewToScene(QtCore.QPointF(x,0)))
        drag(viewport,point(3.50),point(4.1));APP.processEvents()
        self.assertEqual(w.entries,[]);self.assertLess(np.ptp(view.viewRange()[0]),1)
        w.roi_button.click();self.assertEqual(view.mode,'roi')
        drag(viewport,point(3.60),point(3.68));settle(w)
        self.assertEqual(len(w.entries),1);bounds=list(w.entries[0]['bounds'])
        w.roi_button.click();self.assertEqual(view.mode,'zoom')
        drag(viewport,point(3.61),point(3.66));APP.processEvents()
        self.assertEqual(len(w.entries),1);self.assertEqual(w.entries[0]['bounds'],bounds)
        self.assertLess(np.ptp(view.viewRange()[0]),.1)

    def test_spectrum_shortcuts_select_tools_reset_range_and_preserve_text_input(self):
        w = self.w
        w.activateWindow()
        w.spectrum.setFocus()
        APP.processEvents()
        for key, mode in (('i', 'roi'), ('i', 'roi'), ('p', 'pan'), ('z', 'zoom')):
            QtTest.QTest.keyClicks(w.spectrum, key)
            self.assertEqual(w.view.mode, mode)
            self.assertTrue(w.mode_buttons[mode].isChecked())
            self.assertTrue(w.mode_buttons[mode].text().endswith(f'[{key.upper()}]'))

        buttons = w.findChildren(QtWidgets.QPushButton)
        full_range = next(button for button in buttons if button.text() == 'Cały zakres [F]')
        narrow = dict(xRange=(3.6, 3.7), yRange=(0., 1.), padding=0.)
        w.spectrum.setRange(**narrow)
        full_range.click()
        full = w.view.viewRange()
        w.spectrum.setRange(**narrow)
        QtTest.QTest.keyClicks(w.spectrum, 'f')
        np.testing.assert_allclose(w.view.viewRange(), full)
        self.assertGreater(np.ptp(full[0]), 1.)

        w.add_roi([3.60, 3.68])
        settle(w)
        w.roi_table.setFocus()
        QtTest.QTest.keyClicks(w.roi_table, 'p')
        self.assertEqual(w.view.mode, 'pan')
        w.roi_table.editItem(w.roi_table.item(0, GROUP_COLUMN))
        editor = w.roi_table.findChild(QtWidgets.QComboBox).lineEdit()
        editor.setFocus()
        QtTest.QTest.keyClicks(editor, 'izpf')
        self.assertEqual(editor.text(), 'izpf')
        self.assertEqual(w.view.mode, 'pan')
        QtTest.QTest.keyClick(editor, QtCore.Qt.Key.Key_Return)
        settle(w)
        self.assertEqual(w.entries[0]['group'], 'izpf')

        w.language_button.click()
        self.assertEqual(w.roi_button.text(), 'Mouse integration [I]')
        self.assertEqual(w.mode_buttons['pan'].text(), 'Pan [P]')
        self.assertEqual(full_range.text(), 'Full range [F]')
        w.language_button.click()
        self.assertEqual(w.mode_buttons['pan'].text(), 'Przesuwaj [P]')

    def test_selection_feedback_visible_on_press_and_tracks_drag_before_commit(self):
        w=self.w;view=w.view;viewport=w.spectrum.viewport()
        for mode in ('roi','baseline','zoom'):
            with self.subTest(mode=mode):
                w.set_mode(mode);w.spectrum.setXRange(3.4,4.0,padding=0);APP.processEvents()
                def point(x):return w.spectrum.mapFromScene(view.mapViewToScene(QtCore.QPointF(x,0)))
                start,first,end=point(3.60),point(3.65),point(3.75)
                committed=QtTest.QSignalSpy(view.selected)
                before=np.array(view.viewRange());raw=w.dataset.spectrum.copy()
                roi_count=len(w.entries);base_count=len(w.processing.regions)
                QtTest.QTest.mousePress(viewport,QtCore.Qt.MouseButton.LeftButton,pos=start)
                try:
                    APP.processEvents()
                    self.assertIsNotNone(view.selection_preview)
                    self.assertTrue(view.selection_preview.isVisible())
                    anchor=view.selection_preview.lines[0].value()
                    self.assertAlmostEqual(anchor,3.60,delta=.002)
                    with patch('core.fit_decay',side_effect=AssertionError('Unfinished selection triggered fit')):
                        QtTest.QTest.mouseMove(viewport,first,25);APP.processEvents()
                        first_edge=view.selection_preview.lines[1].value()
                        QtTest.QTest.mouseMove(viewport,end,25);APP.processEvents()
                        self.assertEqual(view.selection_preview.lines[0].value(),anchor)
                        self.assertGreater(view.selection_preview.lines[1].value(),first_edge)
                        self.assertTrue(view.selection_preview.lines[1].isVisible())
                        for line,label in zip(view.selection_preview.lines,view.selection_labels):
                            self.assertEqual(label.textItem.toPlainText(),f'{line.value():.3f} ppm')
                        self.assertGreater(view.selection_preview.brush.color().alpha(),0)
                        self.assertLess(view.selection_preview.brush.color().alpha(),100)
                        self.assertEqual(committed.count(),0)
                        self.assertEqual(len(w.entries),roi_count)
                        self.assertEqual(len(w.processing.regions),base_count)
                        np.testing.assert_allclose(view.viewRange(),before)
                        np.testing.assert_array_equal(w.dataset.spectrum,raw)
                finally:
                    QtTest.QTest.mouseRelease(viewport,QtCore.Qt.MouseButton.LeftButton,pos=end)
                settle(w)
                self.assertFalse(view.selection_preview.isVisible())
                self.assertEqual(committed.count(),0 if mode=='zoom' else 1)
                if mode=='roi':self.assertEqual(len(w.entries),roi_count+1)
                if mode=='baseline':self.assertEqual(len(w.processing.regions),base_count+1)
                if mode=='zoom':self.assertLess(np.ptp(view.viewRange()[0]),.2)

    def test_click_without_drag_clears_marker_and_phase_gestures_do_not_select(self):
        w=self.w;viewport=w.spectrum.viewport();pos=viewport.rect().center()
        QtTest.QTest.mousePress(viewport,QtCore.Qt.MouseButton.LeftButton,pos=pos);APP.processEvents()
        self.assertTrue(w.view.selection_preview.isVisible())
        QtTest.QTest.mouseRelease(viewport,QtCore.Qt.MouseButton.LeftButton,pos=pos);APP.processEvents()
        self.assertFalse(w.view.selection_preview.isVisible());self.assertFalse(w.entries)
        w.phase_button.click();APP.processEvents();d=w.phase_window;viewport=d.plot.viewport()
        pos=QtCore.QPoint(round(viewport.width()*.3),round(viewport.height()*.5))
        QtTest.QTest.mousePress(viewport,QtCore.Qt.MouseButton.LeftButton,pos=pos);APP.processEvents()
        self.assertIsNone(d.view.selection_preview)
        QtTest.QTest.mouseRelease(viewport,QtCore.Qt.MouseButton.LeftButton,pos=pos);APP.processEvents()
        d.zoom_button.click()
        pos=pos+QtCore.QPoint(60,0)  # Avoid the pivot that the preceding click placed.
        QtTest.QTest.mousePress(viewport,QtCore.Qt.MouseButton.LeftButton,pos=pos);APP.processEvents()
        self.assertTrue(d.view.selection_preview.isVisible())
        QtTest.QTest.mouseRelease(viewport,QtCore.Qt.MouseButton.LeftButton,pos=pos);APP.processEvents()
        self.assertFalse(d.view.selection_preview.isVisible());d.reject();APP.processEvents()

    def test_mouse_phase_preview_all_layers_pivot_and_explicit_apply(self):
        w=self.w;raw=w.raw.fid.copy();before=w.dataset.phased.copy();options=dict(w.dataset.options)
        w.phase_button.click();APP.processEvents();d=w.phase_window
        self.assertIsNotNone(d);viewport=d.plot.viewport()
        start=QtCore.QPoint(round(viewport.width()*.3),round(viewport.height()*.55))
        finish=start+QtCore.QPoint(0,-50)
        with patch('core.fit_decay',side_effect=AssertionError('Phase preview refitted data')):
            ph0=d.controls['phase0_deg'].value();ph1=d.controls['phase1_deg'].value()
            drag(viewport,start,finish)
            self.assertAlmostEqual(d.controls['phase0_deg'].value()-ph0,25,places=3)
            drag(viewport,start,finish,QtCore.Qt.MouseButton.RightButton)
            self.assertAlmostEqual(d.controls['phase1_deg'].value()-ph1,500,places=3)
            np.testing.assert_array_equal(w.dataset.phased,before)
            old=d.preview.copy();d.set_pivot(4.7)
            np.testing.assert_allclose(d.preview,old,rtol=1e-7,atol=1e-7)
            d.layer.setCurrentIndex(8);APP.processEvents()
            expected=process_fids(w.raw,d.options())
            np.testing.assert_allclose(d.preview,expected.phased,atol=1e-9)
            phase=d.controls['phase0_deg'].value();d.zoom_button.click()
            a=d.plot.mapFromScene(d.view.mapViewToScene(QtCore.QPointF(3.5,0)))
            b=d.plot.mapFromScene(d.view.mapViewToScene(QtCore.QPointF(4.1,0)))
            drag(viewport,a,b)
            self.assertLess(np.ptp(d.view.viewRange()[0]),1)
            self.assertEqual(d.controls['phase0_deg'].value(),phase)
            d.pan_button.click()
            self.assertFalse(d.zoom_button.isChecked())
            d.plot.setYRange(100,300,padding=0)
            prior_range = np.array(d.view.viewRange())
            prior_options = d.options()
            drag(viewport,start,finish)
            moved = np.array(d.view.viewRange())
            np.testing.assert_allclose(np.ptp(moved,axis=1),np.ptp(prior_range,axis=1))
            self.assertLess(moved[1,0],prior_range[1,0])
            QtTest.QTest.mouseClick(viewport,QtCore.Qt.MouseButton.LeftButton,pos=start)
            self.assertEqual(d.options(),prior_options)
            d.pan_button.click()
            self.assertTrue(d.view.phasing)
        applied=d.options();d.apply_button.click();settle(w)
        for key in ('phase0_deg','phase1_deg','pivot_ppm'):self.assertAlmostEqual(w.dataset.options[key],applied[key],places=6)
        self.assertNotEqual(w.dataset.options['phase0_deg'],options['phase0_deg'])
        np.testing.assert_allclose(w.dataset.phased,expected.phased,atol=1e-9)
        np.testing.assert_array_equal(w.raw.fid,raw)

    def test_phase_cancel_keeps_project_and_manual_delay_is_shared(self):
        w=self.w;before=w.dataset.spectrum.copy();options=dict(w.dataset.options)
        w.phase_button.click();APP.processEvents();d=w.phase_window
        ph1=d.controls['phase1_deg'].value();d.delay.setValue(20);d.add_delay()
        self.assertEqual(d.controls['phase1_deg'].value(),ph1+7200)
        d.toggle_fullscreen();APP.processEvents();self.assertTrue(d.isFullScreen())
        d.toggle_fullscreen();APP.processEvents();self.assertFalse(d.isFullScreen())
        d.reject();APP.processEvents()
        np.testing.assert_array_equal(w.dataset.spectrum,before);self.assertEqual(w.dataset.options,options)


if __name__=='__main__':unittest.main()
