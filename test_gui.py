"""Actual Qt widgets, worker updates, groups, point masks and PNG/SVG export."""
import os
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
from pathlib import Path
from unittest.mock import patch
import tempfile,time,unittest
import numpy as np
from PySide6 import QtCore,QtGui,QtWidgets,QtTest
import pyqtgraph as pg
from ui import Window, GROUP_COLUMN, INTEGRAL_COLUMN
from spectrum_settings import SpectrumAppearance, SpectrumSettingsDialog
from core import default_options
from test_dosy_lab import STE,synthetic
from examples.generate_demo import demo_bytes
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

    def test_baseline_preview_apply_and_undo(self):
        w=self.w;w.add_roi([3.60,3.68],group='Sacharoza');settle(w);original=w.dataset.spectrum.copy()
        w.processing.add_region([2.65,2.95]);w.processing.add_region([5.65,6.1]);settle(w)
        self.assertIsNotNone(w.dataset.diagnostics['baseline']);np.testing.assert_array_equal(w.dataset.spectrum,original)
        w.processing.controls['baseline_enabled'].setChecked(True);settle(w)
        self.assertFalse(np.array_equal(w.dataset.spectrum,original))
        w.processing.undo_baseline();settle(w);np.testing.assert_array_equal(w.dataset.spectrum,original)

    def test_statistics_groups_selection_and_diagnostics_never_refit(self):
        w=self.w;w.add_roi([3.60,3.68],name='H1',group='A');w.add_roi([5.26,5.33],name='H2',group='B');settle(w)
        w.show_statistics();s=w.statistics_window;APP.processEvents()
        self.assertEqual(s.summary['N'],0)
        s.group_combo.setCurrentIndex(s.group_combo.findData('A'));self.assertEqual(s.summary['N'],1)
        with patch('core.fit_decay',side_effect=AssertionError('Unexpected refit')):
            s.select_roi('S1');s.table.item(0,0).setCheckState(QtCore.Qt.CheckState.Unchecked)
            self.assertEqual(s.summary['N'],0);s.table.item(0,0).setCheckState(QtCore.Qt.CheckState.Checked)
            self.assertEqual(s.summary['N'],1);s.show_point(2)
        w.processing.controls['lb_Hz'].setValue(.4);settle(w)
        self.assertEqual(s.current['statistics_group'],'A');self.assertEqual(s.summary['N'],1)
        self.assertEqual(s.current['processing']['lb_Hz'],.4)

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

    def test_dosy_map_navigation_coordinates_floor_and_export(self):
        w = self.w
        w.add_roi([3.60, 3.68])
        w.add_roi([5.26, 5.33])
        settle(w)
        with patch('ui.analyze', side_effect=AssertionError('Mapa nie może ponawiać fitu.')):
            w.map_button.click()
            APP.processEvents()
            window = w.map_window
            self.assertTrue(window.isVisible())
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
            window.rois.item(1).setCheckState(QtCore.Qt.CheckState.Unchecked)
            self.assertEqual(len(window.images), 1)
            with tempfile.TemporaryDirectory() as directory:
                png = Path(directory)/'map.png'
                self.assertTrue(window.save_view(png))
                self.assertFalse(QtGui.QImage(str(png)).isNull())
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
        path = str(Path(self.settings_path).with_name('project.json'))
        with patch.object(QtWidgets.QFileDialog, 'getSaveFileName', return_value=(path, '')):
            w.save_project()
        w.open_path(path, project=True)
        settle(w)
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
