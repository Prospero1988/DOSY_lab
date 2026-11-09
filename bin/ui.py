"""DOSY Lab desktop UI: large 1D spectrum, worker processing, cached statistics."""
from copy import deepcopy
from datetime import datetime
from pathlib import Path
import argparse,html,json,sys,uuid
import numpy as np
from PySide6 import QtCore,QtGui,QtWidgets
import pyqtgraph as pg
import pyqtgraph.exporters
from core import process_fids,initial_dataset,analyze
from jeol import read_jdf,manual_b
from processing_panel import ProcessingPanel,BDialog
from statistics_window import StatisticsWindow
from dosy_map import DosyMapWindow
from storage import DATA, PROJECTS, PROJECT_FILTER, choose_save_file, dialog_directory
from plot_export import save_light_plot
from d_statistics import build_snapshot,number
from fullscreen_plot import enable_fullscreen
from exports import export_bundle,save_project,read_project
from spectrum_view import SpectrumBox
from phase_window import PhaseWindow
from spectrum_settings import SpectrumAppearance, SpectrumSettingsDialog
from experiment_info import ExperimentInfo

STYLE='''QMainWindow,QWidget {background:#101827;color:#dbe4ef;font-size:12px;}
QPushButton {background:#26364e;border:1px solid #40516d;padding:7px;border-radius:5px;}
QPushButton:hover {background:#344966;} QPushButton:checked {background:#155e75;}
QPushButton:disabled {color:#66758b;} QComboBox,QDoubleSpinBox,QSpinBox,QLineEdit {padding:4px;background:#26364e;}
QTableWidget,QListWidget {gridline-color:#334155;} QHeaderView::section {background:#26364e;padding:4px;}
QTabBar::tab {padding:6px;background:#26364e;} QTabBar::tab:selected {background:#40516d;}
QDockWidget::title {background:#26364e;padding:7px;}'''
COLORS=['#67e8f9','#fb923c','#4ade80','#f472b6','#facc15','#60a5fa','#c084fc','#f87171']
GROUP_COLUMN = 1
LOWER_PPM_COLUMN = 2
UPPER_PPM_COLUMN = 3
INTEGRAL_COLUMN = 4


class GroupDelegate(QtWidgets.QStyledItemDelegate):
    """Edit a group as text or reuse a label from the current project."""
    def __init__(self, window):
        super().__init__(window.roi_table)
        self.window = window

    def createEditor(self, parent, option, index):
        editor = QtWidgets.QComboBox(parent)
        editor.setEditable(True)
        editor.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.NoInsert)
        editor.addItems(sorted({e['group'] for e in self.window.entries if e['group']}, key=str.casefold))
        editor.activated.connect(lambda _: self.commit_selection(editor))
        return editor

    def setEditorData(self, editor, index):
        editor.setCurrentText(index.data())
        editor.lineEdit().selectAll()

    def setModelData(self, editor, model, index):
        model.setData(index, editor.currentText().strip())

    def commit_selection(self, editor):
        self.commitData.emit(editor)
        self.closeEditor.emit(editor)


class Job(QtCore.QThread):
    completed=QtCore.Signal(object)
    failed=QtCore.Signal(str)
    def __init__(self,fn,parent):super().__init__(parent);self.fn=fn
    def run(self):
        try:self.completed.emit(self.fn())
        except Exception as exc:self.failed.emit(str(exc))


class Window(QtWidgets.QMainWindow):
    def __init__(self, settings=None):
        super().__init__();self.setWindowTitle('DOSY Lab — klasyczne 1D DOSY');self.resize(1560,1000);self.setStyleSheet(STYLE)
        self.settings = settings if settings is not None else QtCore.QSettings('DOSY Lab', 'DOSY Lab')
        self.appearance = SpectrumAppearance.load(self.settings)
        self.raw=None;self.dataset=None;self.entries=[];self.results=[];self.active_id=None
        self.integral_scale = 1.
        self.integral_labels = {}
        self.signal_labels = {}
        self.info_window = None
        self.map_window = None
        self.worker=None;self.job_callback=None;self.pending=False;self.revision=0;self.job_revision=0
        self.valid=False;self._updating=False;self.close_pending=False;self.statistics_window=None;self.phase_window=None
        self.experiment_id=uuid.uuid4().hex;self.region_items={};self.base_items=[];self.next_id=1
        self.timer=QtCore.QTimer(self);self.timer.setSingleShot(True);self.timer.timeout.connect(self.start_calculation)
        root=QtWidgets.QWidget();self.setCentralWidget(root);layout=QtWidgets.QVBoxLayout(root)
        top=QtWidgets.QHBoxLayout();layout.addLayout(top)
        self.open_button = self.button('Otwórz DOSY', self.choose_file, top, 'Otwórz pomiar DOSY z pliku JEOL JDF.')
        self.button('Wczytaj projekt…', self.choose_project, top, 'Przywróć zapisaną analizę wraz z ustawieniami, integracjami i kalibracją całek. Projekt odwołuje się do pliku JDF.')
        self.button('Zapisz projekt', self.save_project, top, 'Zapisz ustawienia analizy, integracje i kalibrację całek w pliku projektu JSON.')
        self.info_button = self.button('Info', self.show_info, top, 'Pokaż parametry eksperymentu zapisane w JDF i aktualne ustawienia widma.')
        top.addStretch()
        self.button('Lista b', self.edit_b, top, 'Sprawdź lub zmień wartości b przypisane warstwom gradientowym.')
        self.button('Statystyka D', self.show_statistics, top, 'Porównaj współczynniki dyfuzji w obrębie grup związków.')
        self.export_button = self.button('Eksport Statystyki D', self.choose_export, top, 'Zapisz pakiet wyników, statystyki, tabele i wykresy.')
        self.export_button.setEnabled(False)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        row.addWidget(QtWidgets.QLabel('Warstwa:'));self.layer=QtWidgets.QComboBox();row.addWidget(self.layer,1)
        self.layer.setToolTip('Wybierz warstwę gradientową wyświetlaną na widmie i w kolumnie Całka.')
        self.layer.currentIndexChanged.connect(self.refresh_displayed_layer)
        self.button('Pokaż referencję',lambda:self.layer.setCurrentIndex(int(self.processing.options()['reference_layer'])),row,
                    'Pokaż warstwę referencyjną ustawioną w panelu Processing.')
        self.overlay=QtWidgets.QCheckBox('Nałóż warstwy:');row.addWidget(self.overlay);self.overlay.toggled.connect(self.refresh_displayed_layer)
        self.overlay.setToolTip('Porównaj kilka warstw na jednym widmie.')
        self.overlay_layers=QtWidgets.QLineEdit('1,8,16');self.overlay_layers.setMaximumWidth(110);row.addWidget(self.overlay_layers)
        self.overlay_layers.setToolTip('Wpisz numery warstw oddzielone przecinkami. Pierwsza warstwa wyznacza całki pokazane przy widmie.')
        self.overlay_layers.editingFinished.connect(self.refresh_displayed_layer)
        self.imaginary=QtWidgets.QCheckBox('Pokaż Im');row.addWidget(self.imaginary);self.imaginary.toggled.connect(self.render_spectrum)
        self.imaginary.setToolTip('Pokaż część urojoną widma obok części rzeczywistej.')
        self.coordinate=QtWidgets.QLabel('ppm —');row.addWidget(self.coordinate)
        spectrum_tools = QtWidgets.QHBoxLayout()
        layout.addLayout(spectrum_tools)
        self.mode_buttons = {}
        self.roi_button = self.tool_button('roi','Integracja myszą',spectrum_tools,
                                           'Przeciągnij po widmie, aby zaznaczyć granice nowej całki.')
        self.button('Dodaj całkę ręcznie', self.add_numeric_roi, spectrum_tools, 'Wpisz dolną i górną granicę całkowania w ppm.')
        self.button('Usuń całkę', self.delete_roi, spectrum_tools, 'Usuń całkę wybraną w tabeli.')
        for key,label,hint in (
            ('modify','Modyfikuj całkę','Przeciągnij środek obszaru, aby go przesunąć, lub krawędź, aby zmienić granicę.'),
            ('split','Podziel całkę','Przeciągnij jedną z krawędzi do środka obszaru i puść w miejscu podziału.'),
            ('zoom','Zoom','Przeciągnij, aby powiększyć zakres ppm.'),
            ('pan','Przesuwaj','Przeciągnij, aby przesunąć widok widma.'),
            ('pivot','Ustaw pivot','Kliknij widmo, aby ustawić punkt obrotu korekcji fazy.'),
            ('baseline','Zaznacz bazę','Zaznacz fragmenty bez pików, z których program wyznaczy linię bazową. Żółta linia to podgląd.'),
        ):
            self.tool_button(key,label,spectrum_tools,hint)
        self.clear_baseline_button = self.button('Wyczyść bazę',self.clear_baseline,spectrum_tools,
                                                 'Usuń wszystkie niebieskie obszary bazy i wyłącz jej odejmowanie.')
        self.clear_baseline_button.setEnabled(False)
        self.map_button = self.button('Mapa DOSY 2D', self.show_map, spectrum_tools,
                                      'Pokaż mapę ppm × D z całek i aktualnych dopasowań.')
        self.map_button.setEnabled(False)
        self.mode_buttons['zoom'].setChecked(True)
        spectrum_tools.addStretch()
        self.button('Processing', lambda:self.processing.setVisible(not self.processing.isVisible()), spectrum_tools,
                    'Pokaż ustawienia przetwarzania widma i korekcji bazy.')
        self.phase_button = self.button('Korekcja fazy', self.show_phase, spectrum_tools, 'Dostrój wspólną fazę wszystkich warstw za pomocą myszy.')
        self.button('Widok', self.edit_spectrum_settings, spectrum_tools, 'Ustaw kolor i grubość linii widma oraz widoczność siatki.')
        self.button('Cały zakres', lambda:self.spectrum.autoRange(), spectrum_tools, 'Przywróć cały zakres widma. Działa też dwuklik.')
        self.fullscreen_button = self.button('Pełny ekran', lambda:self.spectrum.fullscreen_controller.open(), spectrum_tools,
                                             'Otwórz widmo na pełnym ekranie. Esc wraca do głównego okna.')
        self.save_spectrum_button = self.button('Zapisz PNG', self.choose_spectrum_png, spectrum_tools,
                                                'Zapisz aktualny widok widma 1D w jasnej kolorystyce.')
        self.save_spectrum_button.setEnabled(False)
        vertical=QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical);layout.addWidget(vertical,1)
        self.view=SpectrumBox();self.spectrum=pg.PlotWidget(viewBox=self.view);vertical.addWidget(self.spectrum)
        self.view.bind_selection_feedback(self.spectrum)
        self.spectrum.setLabel('bottom','¹H — przesunięcie chemiczne',units='ppm');self.spectrum.setLabel('left','Re / Im',units='j.u.')
        self.spectrum.invertX(True)
        self.appearance.apply_grid(self.spectrum)
        self.view.selected.connect(self.range_selected);self.view.clicked.connect(self.pivot_clicked)
        self.view.split_requested.connect(self.split_roi)
        self.view.sigYRangeChanged.connect(self.position_integral_labels)
        self.spectrum.scene().sigMouseMoved.connect(self.mouse_moved)
        bottom=QtWidgets.QSplitter();vertical.addWidget(bottom)
        left=QtWidgets.QWidget();ll=QtWidgets.QVBoxLayout(left);bottom.addWidget(left)
        bar=QtWidgets.QHBoxLayout();ll.addLayout(bar)
        self.button('Dopasuj / odśwież',self.request_calculation,bar,'Przelicz widma, całki i dopasowania D z bieżących ustawień.')
        self.roi_table=QtWidgets.QTableWidget(0,9)
        self.roi_table.setHorizontalHeaderLabels(['ID','Grupa / związek','Od ppm','Do ppm','Całka','D / 10⁻¹⁰','SE(D) / 10⁻¹⁰','R²','QC'])
        self.roi_table.setItemDelegateForColumn(GROUP_COLUMN, GroupDelegate(self))
        self.roi_table.horizontalHeaderItem(6).setToolTip(
            'Błąd standardowy oszacowania D z dopasowania zaniku, w 10⁻¹⁰ m²/s.')
        self.roi_table.horizontalHeaderItem(INTEGRAL_COLUMN).setToolTip(
            'Całki pierwszej widocznej warstwy. Wpisz wartość odniesienia, aby proporcjonalnie przeskalować wszystkie całki.')
        self.roi_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.roi_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.roi_table.itemChanged.connect(self.roi_edited);self.roi_table.currentCellChanged.connect(self.roi_selected);ll.addWidget(self.roi_table)
        self.fit_info=QtWidgets.QLabel('Zaznacz sygnał. Grupa jest potrzebna do statystyki zbiorczej.');self.fit_info.setWordWrap(True);ll.addWidget(self.fit_info)
        self.tabs=QtWidgets.QTabWidget();bottom.addWidget(self.tabs)
        decay_widget=QtWidgets.QWidget();dl=QtWidgets.QVBoxLayout(decay_widget)
        row=QtWidgets.QHBoxLayout();dl.addLayout(row)
        self.all_decays=QtWidgets.QCheckBox('Wszystkie ROI');self.all_decays.setChecked(True);row.addWidget(self.all_decays)
        self.all_decays.setToolTip('Pokaż krzywe zaniku wszystkich integracji. Po wyłączeniu widać tylko wybraną.')
        self.normalized=QtWidgets.QCheckBox('Podgląd I / pierwsza całka');self.normalized.setChecked(True);row.addWidget(self.normalized)
        self.normalized.setToolTip('Wyświetl każdy zanik w jednostkach jego pierwszej całki.')
        self.all_decays.toggled.connect(self.render_diagnostics);self.normalized.toggled.connect(self.render_diagnostics)
        row.addStretch()
        self.button('Pełny ekran zaników', lambda:self.decay.fullscreen_controller.open(), row, 'Otwórz wykres zaników na pełnym ekranie.')
        self.button('Pełny ekran reszt', lambda:self.residual.fullscreen_controller.open(), row, 'Otwórz wykres reszt na pełnym ekranie.')
        self.decay=pg.PlotWidget();self.decay.addLegend();self.decay.showGrid(x=True,y=True,alpha=.2);dl.addWidget(self.decay,3)
        self.residual=pg.PlotWidget();self.residual.showGrid(x=True,y=True,alpha=.2);self.residual.setXLink(self.decay);dl.addWidget(self.residual,1)
        for plot in (self.decay,self.residual):
            view = plot.getViewBox()
            view.setMouseMode(pg.ViewBox.RectMode)
            view.mouse_help = 'Przeciągnij prostokąt, aby powiększyć obszar. Dwuklik przywraca cały zakres.'
        self.decay.setLabel('bottom','b / 10⁹ [s/m²]');self.residual.setLabel('bottom','b / 10⁹ [s/m²]')
        self.tabs.addTab(decay_widget,'Zaniki i reszty')
        points=QtWidgets.QWidget();pl=QtWidgets.QVBoxLayout(points)
        self.weighted=QtWidgets.QCheckBox('Fit ważony: użyj podanych σ(I) aktywnych punktów');pl.addWidget(self.weighted)
        self.weighted.setToolTip('Użyj niepewności całek wpisanych w tabeli przy dopasowaniu D.')
        self.weighted.toggled.connect(self.weighting_changed)
        self.point_table=QtWidgets.QTableWidget(0,8)
        self.point_table.setHorizontalHeaderLabels(['Fit','Warstwa','g [T/m]','b [s/m²]','Całka','Przewidywanie','Reszta','σ(I) użytkownika'])
        self.point_table.itemChanged.connect(self.point_changed);self.point_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        pl.addWidget(self.point_table);self.tabs.addTab(points,'Całki / punkty fitu')
        self.roi_overlay=pg.PlotWidget();self.roi_overlay.setLabel('bottom','¹H [ppm]');self.roi_overlay.invertX(True)
        self.appearance.apply_grid(self.roi_overlay)
        self.tabs.addTab(self.roi_overlay,'Widma w aktywnym ROI')
        vertical.setSizes([560,330]);bottom.setSizes([680,750])
        self.processing=ProcessingPanel(self);self.addDockWidget(QtCore.Qt.DockWidgetArea.LeftDockWidgetArea,self.processing)
        self.processing.changed.connect(self.request_calculation);self.processing.mode_requested.connect(self.set_mode)
        self.statusBar().showMessage('Przeciągnięcie powiększa zakres. Rolka zmienia wysokość widma. Dwuklik przywraca cały zakres.')
        context = lambda:html.escape(self.dataset.name) if self.dataset else ''
        enable_fullscreen(self.spectrum,'Widmo 1D',context,reset_view=self.spectrum.autoRange)
        for plot,title in [(self.decay,'Zaniki ROI'),(self.residual,'Reszty')]:
            enable_fullscreen(plot,title,context,reset_view=self.reset_decay_views)
        enable_fullscreen(self.roi_overlay,'Nakładanie widm w ROI',context)

    @staticmethod
    def button(text,fn,layout,tooltip):
        b=QtWidgets.QPushButton(text);b.setToolTip(tooltip);b.clicked.connect(fn);layout.addWidget(b);return b

    def tool_button(self,key,text,layout,tooltip):
        button = self.button(text,lambda checked:self.set_mode(key if checked else 'zoom'),layout,tooltip)
        button.setCheckable(True)
        self.mode_buttons[key] = button
        return button
    @property
    def active(self):return next((e for e in self.entries if e['id']==self.active_id),None)
    @property
    def active_result(self):return next((e for e in self.results if e['id']==self.active_id),None)

    def choose_file(self):
        try:
            path,_=QtWidgets.QFileDialog.getOpenFileName(self,'Otwórz klasyczny JEOL DOSY',dialog_directory(DATA),'JEOL (*.jdf)')
            if path:self.open_path(path)
        except OSError as exc:
            QtWidgets.QMessageBox.warning(self,'Otwórz DOSY',str(exc))

    def open_path(self,path,project=False):
        if self.worker is not None:return
        def work():
            if project:
                raw,o,entries,scale=read_project(path,with_integral_scale=True)
                ds=process_fids(raw,o,[e['bounds'] for e in entries])
            else:
                raw=read_jdf(path);entries=[];ds=initial_dataset(raw);scale=1.
            return raw,ds,entries,analyze(ds,entries),scale
        self.start_job(work,self.loaded)

    def start_job(self,fn,callback):
        self.valid=False;self.export_button.setEnabled(False);self.open_button.setEnabled(False)
        self.invalidate_map()
        self.statusBar().showMessage('Obliczanie w tle. Widok odświeży się po zakończeniu.')
        self.job_callback=callback;self.worker=Job(fn,self)
        self.worker.completed.connect(self.job_success);self.worker.failed.connect(self.job_error);self.worker.finished.connect(self.job_finished)
        self.worker.start()

    @QtCore.Slot(object)
    def job_success(self,value):self.job_callback(value)
    @QtCore.Slot(str)
    def job_error(self,error):self.statusBar().showMessage('Nie zastosowano zmian. '+error)
    @QtCore.Slot()
    def job_finished(self):
        old=self.worker;self.worker=None;old.deleteLater();self.open_button.setEnabled(True)
        if self.close_pending:self.close();return
        if self.pending:self.pending=False;self.timer.start(1)

    def loaded(self,payload):
        self.raw,self.dataset,self.entries,self.results,self.integral_scale=payload
        self.experiment_id=uuid.uuid4().hex
        self.active_id=self.entries[0]['id'] if self.entries else None
        self.next_id=max([int(e['id'][1:]) for e in self.entries if e['id'].startswith('S') and e['id'][1:].isdigit()]+[0])+1
        self.processing.set_options(self.dataset.options,len(self.raw.fid))
        self.layer.blockSignals(True);self.layer.clear()
        for i,g in enumerate(self.raw.gradients):
            b='—' if self.raw.b is None else f'{self.raw.b[i]:.4e}'
            self.layer.addItem(f'{i+1}: g={g:.5g} T/m · b={b} s/m²',i)
        self.layer.blockSignals(False);self.layer.setCurrentIndex(self.dataset.options['reference_layer'])
        self.overlay_layers.setText(','.join(str(i+1) for i in np.unique([0,len(self.raw.fid)//2,len(self.raw.fid)-1])))
        self.set_mode('zoom')
        self.accepted();self.spectrum.autoRange()

    def request_calculation(self,*args):
        if self.raw is None or self._updating:return
        self.revision+=1;self.valid=False;self.export_button.setEnabled(False);self.timer.start(180)
        self.invalidate_map()

    def invalidate_map(self):
        self.map_button.setEnabled(False)
        if self.map_window is not None:self.map_window.invalidate()

    def start_calculation(self):
        if self.raw is None:return
        if self.worker is not None:self.pending=True;return
        raw=self.raw;options=self.processing.options();entries=deepcopy(self.entries);rev=self.revision
        def work():
            ds=process_fids(raw,options,[e['bounds'] for e in entries]);return rev,ds,analyze(ds,entries)
        self.start_job(work,self.calculated)

    def calculated(self,payload):
        rev,ds,results=payload
        if rev!=self.revision:self.pending=True;return
        previous={e['id']:e.get('fit',{}).get('D_m2_s') for e in self.results}
        for e in results:e['previous_D']=previous.get(e['id'])
        self.dataset=ds;self.results=results;self.accepted()

    def accepted(self):
        self.valid=True;self.export_button.setEnabled(True)
        self.save_spectrum_button.setEnabled(True)
        self.map_button.setEnabled(True)
        meta=self.raw.metadata
        delay=self.dataset.options['filter_delay_samples']
        self.layer.blockSignals(True)
        for i,g in enumerate(self.raw.gradients):
            b='—' if self.raw.b is None else f'{self.raw.b[i]:.4e}'
            self.layer.setItemText(i,f'{i+1}: g={number(g)} T/m · b={b} s/m²')
        self.layer.blockSignals(False)
        self.render_spectrum();self.render_roi_table();self.render_diagnostics()
        if self.info_window is not None and self.info_window.isVisible():
            self.info_window.set_dataset(self.dataset,self.integral_scale)
        if self.statistics_window is not None:self.statistics_window.update_live(self.snapshot())
        if self.map_window is not None and self.map_window.isVisible():
            self.map_window.update_data(self.experiment_id,self.dataset,self.results)
        self.statusBar().showMessage('; '.join(self.dataset.diagnostics['warnings']) or 'Gotowe. Zmiana wyboru widma lub ROI nie zmienia danych pomiarowych.')
        if meta.get('digital_filter') and not delay and abs(self.dataset.options['phase1_deg'])<1:
            self.statusBar().showMessage('Brak zastosowanej korekcji opóźnienia filtra. Sprawdź fazę; falowanie może pochodzić z filtra cyfrowego.')

    def selected_layers(self):
        if self.raw is None:return []
        try:
            if self.overlay.isChecked():
                values=[int(x.strip())-1 for x in self.overlay_layers.text().split(',')]
                if not values or any(i<0 or i>=len(self.raw.fid) for i in values):raise ValueError()
                return list(dict.fromkeys(values))
        except ValueError:self.statusBar().showMessage('Podaj numery istniejących warstw rozdzielone przecinkami.')
        return [max(0,self.layer.currentIndex())]

    def render_spectrum(self,*args):
        if self.dataset is None:return
        self.integral_labels = {}
        self.signal_labels = {}
        self.spectrum.clear();self.region_items={};self.base_items=[]
        self.view.split_regions = self.region_items
        self.clear_baseline_button.setEnabled(bool(self.processing.regions))
        d=self.dataset
        for k,i in enumerate(self.selected_layers()):
            color = COLORS[k % len(COLORS)] if k else self.appearance.color
            self.spectrum.plot(d.ppm,d.spectrum[i].real,pen=self.appearance.pen(color))
            if self.imaginary.isChecked():self.spectrum.plot(d.ppm,d.spectrum[i].imag,pen=self.appearance.pen(color,QtCore.Qt.PenStyle.DotLine))
            if self.processing.regions and d.options['baseline_preview'] and d.diagnostics['baseline']:
                self.spectrum.plot(d.ppm,d.baseline[i],pen=pg.mkPen('#fbbf24',width=2,style=QtCore.Qt.PenStyle.DashLine))
                if d.options['baseline_enabled']:self.spectrum.plot(d.ppm,d.phased[i].real,pen=self.appearance.pen('#64748b',QtCore.Qt.PenStyle.DotLine))
        for e in self.entries:
            movable = self.view.mode=='modify'
            color = QtGui.QColor(e['color'])
            color.setAlpha(24)
            item=pg.LinearRegionItem(e['bounds'],movable=movable,brush=pg.mkBrush(color),pen=pg.mkPen(e['color']))
            buttons = QtCore.Qt.MouseButton.LeftButton if movable else QtCore.Qt.MouseButton.NoButton
            item.setAcceptedMouseButtons(buttons)
            for line in item.lines:line.setAcceptedMouseButtons(buttons)
            item.setZValue(5);self.spectrum.addItem(item);self.region_items[e['id']]=item
            item.sigRegionChangeFinished.connect(lambda region,sid=e['id']:self.region_changed(sid,region.getRegion()))
            for labels,text in ((self.signal_labels,e['id']),(self.integral_labels,'')):
                label = pg.TextItem(text, anchor=(.5,0), color=e['color'], fill=pg.mkBrush(16,24,39,220))
                label.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
                label.textItem.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
                label.setZValue(6)
                self.spectrum.addItem(label,ignoreBounds=True)
                labels[e['id']] = label
            item.sigRegionChanged.connect(self.position_integral_labels)
        for bounds in self.processing.regions:
            item=pg.LinearRegionItem(bounds,movable=False,brush=pg.mkBrush(60,110,200,35),pen=pg.mkPen('#3b82f6'))
            item.setZValue(-2);self.spectrum.addItem(item);self.base_items.append(item)
        self.pivot_line=pg.InfiniteLine(self.processing.controls['pivot_ppm'].value(),movable=self.view.mode=='pivot',pen=pg.mkPen('#c084fc',style=QtCore.Qt.PenStyle.DashLine))
        self.pivot_line.setAcceptedMouseButtons(QtCore.Qt.MouseButton.LeftButton if self.view.mode=='pivot' else QtCore.Qt.MouseButton.NoButton)
        self.pivot_line.sigPositionChangeFinished.connect(lambda line:self.pivot_clicked(line.value()))
        self.spectrum.addItem(self.pivot_line)
        self.update_integral_labels()

    def integral_value(self, entry, result):
        if result is None or list(result['bounds']) != list(entry['bounds']):
            return None
        integrals = (result.get('result') or {}).get('corrected_integral')
        layers = self.selected_layers()
        if integrals is None or not layers:
            return None
        return float(integrals[layers[0]])

    def format_integral(self, value):
        if value is None:
            return '—'
        scaled = value * self.integral_scale
        return f'{scaled:.2f}'.replace('.', ',') if np.isfinite(scaled) else '—'

    def position_integral_labels(self, *args):
        low, high = self.view.viewRange()[1]
        layers = self.selected_layers()
        if not self.integral_labels or not layers:
            return
        span = high-low
        d = self.dataset
        for sid,label in self.integral_labels.items():
            left,right = self.region_items[sid].getRegion()
            center = (left+right)/2
            self.signal_labels[sid].setPos(center,high-.04*span)
            start,end = np.searchsorted(d.ppm,[left,right])
            trace = d.spectrum[layers,start:end]
            floor = min(0., float(trace.real.min())) if trace.size else 0.
            if trace.size and self.imaginary.isChecked():
                floor = min(floor,float(trace.imag.min()))
            y = max(low+.06*span,min(high-.06*span,floor-.025*span))
            label.setPos(center,y)

    def update_integral_labels(self):
        results = {result['id']:result for result in self.results}
        for entry in self.entries:
            sid = entry['id']
            label = self.integral_labels.get(sid)
            if label is None:
                continue
            value = self.format_integral(self.integral_value(entry,results.get(sid)))
            label.setText(value)
            label.setVisible(value != '—')
        self.position_integral_labels()

    def refresh_displayed_layer(self,*args):
        self.render_spectrum()
        self.render_roi_table()
        self.render_diagnostics()

    def reset_decay_views(self):
        self.decay.autoRange(padding=.05)
        self.residual.autoRange(padding=.05)
        self.residual.getViewBox().linkedXChanged()

    def set_mode(self,key):
        for mode,button in self.mode_buttons.items():
            button.setChecked(mode==key)
        if self.view.mode==key:
            return
        self.view.hide_selection_preview()
        self.view.mode=key
        self.render_spectrum()

    def clear_baseline(self):
        self.processing.clear_regions()
        if self.view.mode=='baseline':
            self.set_mode('zoom')
        else:
            self.render_spectrum()
    def pivot_clicked(self,value):self.processing.controls['pivot_ppm'].setValue(value)
    def mouse_moved(self,pos):
        if self.spectrum.sceneBoundingRect().contains(pos):self.coordinate.setText(f'{self.view.mapSceneToView(pos).x():.5f} ppm')
    def range_selected(self,a,b):
        if self.dataset is None:return
        if self.view.mode=='baseline':self.processing.add_region([a,b])
        else:self.add_roi([a,b])
    def add_numeric_roi(self):
        if self.dataset is None:return
        text,ok=QtWidgets.QInputDialog.getText(self,'Dodaj całkę ręcznie','Wpisz dwie granice całkowania w ppm, rozdzielone średnikiem.')
        if ok:
            try:self.add_roi([float(v.strip().replace(',', '.')) for v in text.split(';')])
            except ValueError as e:self.statusBar().showMessage(str(e),5000)

    def new_roi(self,bounds,name=None,group=''):
        if len(bounds)!=2 or not np.isfinite(bounds).all() or min(bounds)>=max(bounds):raise ValueError('Nieprawidłowe granice ROI.')
        sid=f'S{self.next_id}';self.next_id+=1
        return dict(id=sid,name=name or sid,group=group,bounds=sorted(bounds),color=COLORS[(self.next_id-2)%len(COLORS)],
                    fit_mask=[True]*len(self.raw.fid),sigma=[None]*len(self.raw.fid),weighted=False)

    def add_roi(self,bounds,name=None,group=''):
        if self.raw is None:return
        entry = self.new_roi(bounds,name,group)
        self.entries.append(entry)
        self.active_id=entry['id'];self.render_spectrum();self.render_roi_table();self.request_calculation()

    def split_roi(self,sid,cut):
        index = next((i for i,entry in enumerate(self.entries) if entry['id']==sid),None)
        if index is None:
            return
        entry = self.entries[index]
        left,right = entry['bounds']
        if not left<cut<right:
            return
        ppm = self.dataset.ppm
        ranges = ([left,cut],[cut,right])
        for lo,hi in ranges:
            points = np.searchsorted(ppm,hi,side='right')-np.searchsorted(ppm,lo)
            if points<2:
                self.statusBar().showMessage('Każda część całki musi obejmować co najmniej dwa punkty widma.',5000)
                return
        children = [self.new_roi(bounds,f"{entry['name']} ({part})",entry['group'])
                    for part,bounds in enumerate(ranges,1)]
        for child in children:
            child['fit_mask'] = list(entry['fit_mask'])
        self.entries[index:index+1] = children
        self.active_id = children[0]['id']
        self.render_spectrum()
        self.render_roi_table()
        self.request_calculation()
    def region_changed(self,sid,bounds):
        if self._updating:return
        e=next(e for e in self.entries if e['id']==sid)
        lo,hi=sorted(map(float,bounds))
        ppm=self.dataset.ppm
        points=np.searchsorted(ppm,hi,side='right')-np.searchsorted(ppm,lo)
        if lo<ppm[0] or hi>ppm[-1] or points<2:
            region=self.region_items[sid]
            with QtCore.QSignalBlocker(region):
                region.setRegion(e['bounds'])
            self.position_integral_labels()
            self.statusBar().showMessage('Granice całki muszą obejmować co najmniej dwa punkty i mieścić się w widmie.',5000)
            return
        if e['bounds']==[lo,hi]:
            return
        e['bounds']=[lo,hi]
        self.active_id=sid
        self.render_roi_table();self.update_integral_labels();self.request_calculation()
    def delete_roi(self):
        if self.active_id:
            self.entries=[e for e in self.entries if e['id']!=self.active_id]
            self.active_id=self.entries[0]['id'] if self.entries else None
            self.render_spectrum();self.render_roi_table();self.request_calculation()

    def render_roi_table(self):
        self._updating=True;self.roi_table.setRowCount(len(self.entries));fits={e['id']:e for e in self.results}
        for i,e in enumerate(self.entries):
            fit=fits.get(e['id'],{}).get('fit') or {}
            integral=self.integral_value(e,fits.get(e['id']))
            values=[e['id'],e['group'],*[f'{x:.6f}' for x in e['bounds']],
                    self.format_integral(integral),
                    f"{fit['D_m2_s']/1e-10:.5f}" if fit.get('ok') else '—',
                    f"{fit['D_standard_error_m2_s']/1e-10:.5f}" if fit.get('D_standard_error_m2_s') is not None else '—',
                    f"{fit['R2']:.6f}" if fit.get('R2') is not None else '—',
                    ('UWAGA' if fit.get('warnings') else 'OK') if fit.get('ok') else 'BRAK FITU']
            for j,value in enumerate(values):
                item=QtWidgets.QTableWidgetItem(value)
                editable = j in (GROUP_COLUMN,LOWER_PPM_COLUMN,UPPER_PPM_COLUMN) or (j==INTEGRAL_COLUMN and integral is not None)
                if not editable:item.setFlags(item.flags()&~QtCore.Qt.ItemFlag.ItemIsEditable)
                if j==0 and e['name']!=e['id']:item.setToolTip(e['name'])
                if j==GROUP_COLUMN:
                    item.setToolTip('Dwuklik otwiera edycję. Wpisz nową grupę lub wybierz istniejącą z listy.')
                if j==INTEGRAL_COLUMN:
                    item.setToolTip('Wpisz wartość odniesienia. Wszystkie całki przy widmie i w tej kolumnie zmienią się w tej samej proporcji.')
                self.roi_table.setItem(i,j,item)
            if e['id']==self.active_id:self.roi_table.setCurrentCell(i,0)
        self._updating=False
    def roi_edited(self,item):
        if self._updating:return
        e=self.entries[item.row()]
        try:
            if item.column()==INTEGRAL_COLUMN:
                if not self.valid or self.worker is not None:
                    raise ValueError('Poczekaj na aktualne całki przed kalibracją.')
                result=next((r for r in self.results if r['id']==e['id']),None)
                source=self.integral_value(e,result)
                target=float(item.text().strip().replace(',','.'))
                scale=target/source if source else 0.
                if not np.isfinite(scale) or scale<=0:
                    raise ValueError('Wpisz skończoną, niezerową wartość o znaku zmierzonej całki.')
                self.integral_scale=scale
                self.render_roi_table()
                self.update_integral_labels()
                if self.info_window is not None and self.info_window.isVisible():
                    self.info_window.set_dataset(self.dataset,self.integral_scale)
                return
            if item.column()==GROUP_COLUMN:e['group']=item.text().strip()
            elif item.column() in (LOWER_PPM_COLUMN,UPPER_PPM_COLUMN):
                e['bounds'][item.column()-LOWER_PPM_COLUMN]=float(item.text());e['bounds']=sorted(e['bounds'])
            self.render_spectrum();self.request_calculation()
        except Exception as exc:self.statusBar().showMessage(str(exc),5000);self.render_roi_table()
    def roi_selected(self,row,col,oldrow,oldcol):
        if not self._updating and 0<=row<len(self.entries):self.active_id=self.entries[row]['id'];self.render_diagnostics()

    def render_diagnostics(self,*args):
        self.decay.clear();self.residual.clear();self.render_roi_overlay();self._updating=True;self.point_table.setRowCount(0)
        d=self.dataset
        if d is None:self._updating=False;return
        normalized=self.normalized.isChecked();self.decay.setLabel('left','I / pierwsza całka' if normalized else 'Całka [j.u.·ppm]')
        self.residual.setLabel('left','Reszta / Iref' if normalized else 'Reszta')
        if d.b is not None:
            x=d.b/1e9;x_label='b / 10⁹ [s/m²]'
        elif np.isfinite(d.gradients).all():
            x=d.gradients;x_label='g [T/m] · brak zweryfikowanego b'
        else:
            x=np.arange(1,len(d.raw.fid)+1);x_label='Warstwa akwizycji · brak zweryfikowanego b'
        for plot in (self.decay,self.residual):plot.setLabel('bottom',x_label)
        self.residual.addItem(pg.InfiniteLine(0,angle=0,pen=pg.mkPen('#64748b')))
        for entry in self.results:
            if not self.all_decays.isChecked() and entry['id']!=self.active_id:continue
            result=entry.get('result') or {};fit=entry.get('fit') or {};y=result.get('corrected_integral')
            if y is None:continue
            factor=float(y[0]) if normalized and y[0]!=0 else 1.
            mask=np.asarray(fit.get('fit_mask',entry['fit_mask']),bool);color=entry['color']
            for state,symbol in ((True,'o'),(False,'x')):
                spots=[dict(pos=(x[i],y[i]/factor),data=(entry['id'],i)) for i in range(len(y)) if mask[i]==state]
                scatter=pg.ScatterPlotItem(spots=spots,symbol=symbol,size=8,brush=pg.mkBrush(color),pen=pg.mkPen(color),hoverable=True,
                    tip=lambda x,y,data:self.point_tooltip(*data),useCache=False)
                scatter.sigClicked.connect(self.decay_point_clicked);self.decay.addItem(scatter)
            if fit.get('ok'):self.decay.plot(fit['b_grid']/1e9,fit['curve']/factor,pen=pg.mkPen(color,width=1.5),name=entry['name'])
        active=self.active_result
        if active:
            e=self.active;fit=active.get('fit') or {};result=active.get('result') or {};y=result.get('corrected_integral')
            self.weighted.setChecked(e.get('weighted',False));self.point_table.setRowCount(len(d.raw.fid))
            if y is not None:
                pred=fit.get('prediction',[None]*len(y));res=fit.get('residuals',[None]*len(y))
                for i in range(len(y)):
                    check=QtWidgets.QTableWidgetItem();check.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled|QtCore.Qt.ItemFlag.ItemIsUserCheckable|QtCore.Qt.ItemFlag.ItemIsSelectable)
                    check.setCheckState(QtCore.Qt.CheckState.Checked if e['fit_mask'][i] else QtCore.Qt.CheckState.Unchecked);self.point_table.setItem(i,0,check)
                    values=[str(i+1),number(d.gradients[i]),number(d.b[i]) if d.b is not None else '—',number(y[i]),number(pred[i]),number(res[i]),number(e['sigma'][i]) if e['sigma'][i] is not None else '']
                    for j,v in enumerate(values,1):
                        item=QtWidgets.QTableWidgetItem(v)
                        if j!=7:item.setFlags(item.flags()&~QtCore.Qt.ItemFlag.ItemIsEditable)
                        self.point_table.setItem(i,j,item)
                if fit.get('ok') and d.b is not None:
                    factor=y[0] if normalized and y[0]!=0 else 1
                    for state,sym in ((True,'o'),(False,'x')):
                        mask=np.asarray(e['fit_mask'])==state;self.residual.plot(d.b[mask]/1e9,np.asarray(res)[mask]/factor,pen=None,symbol=sym,symbolBrush=e['color'],symbolSize=6)
            info=f"{e['id']} · {e['name']} · grupa: {e['group'] or 'nieprzypisana'}\n"
            if fit.get('ok'):
                previous=active.get('previous_D');change=f' · poprzednio {previous/1e-10:.4f}' if previous is not None and previous!=fit['D_m2_s'] else ''
                info+=f"D={fit['D_m2_s']/1e-10:.5f} ×10⁻¹⁰ m²/s{change}\nSE={number(fit.get('D_standard_error_m2_s'))} m²/s · CI95={fit.get('D_CI95_m2_s')} · R²={fit.get('R2')}\n"
                info+=' '.join(fit.get('warnings',[]))
            else:
                info+=fit.get('message','Brak poprawnego fitu.')
                if d.b is None:
                    info+='\n'+d.metadata['b_source']['description']
            self.fit_info.setText(info)
        self._updating=False

    def render_roi_overlay(self, reset_range=True):
        self.roi_overlay.clear()
        if self.dataset is None or self.active_result is None:
            return
        d = self.dataset
        lo, hi = sorted(self.active['bounds'])
        mask = (d.ppm >= lo) & (d.ppm <= hi)
        for k, i in enumerate(self.selected_layers()):
            color = COLORS[k % len(COLORS)] if k else self.appearance.color
            self.roi_overlay.plot(d.ppm[mask], d.spectrum[i, mask].real, pen=self.appearance.pen(color))
        if reset_range:
            self.roi_overlay.setXRange(lo, hi, padding=.02)

    def point_tooltip(self,sid,index):
        e=next((e for e in self.results if e['id']==sid),None)
        if e is None:return ''
        f=e.get('fit') or {};y=e.get('result',{}).get('corrected_integral',[])
        b=number(self.raw.b[index])+' s/m²' if self.raw.b is not None else 'brak'
        return f"{sid}, warstwa {index+1}\ng={number(self.raw.gradients[index])} T/m · b={b}\nI={y[index]:.6g} · przewidywanie={number(f.get('prediction',[None]*len(y))[index])}\nreszta={number(f.get('residuals',[None]*len(y))[index])}"
    def decay_point_clicked(self,item,points,event):
        if not points:return
        sid,index=points[0].data();self.active_id=sid;self.render_roi_table();self.render_diagnostics();self.point_table.setCurrentCell(index,0);self.statusBar().showMessage(self.point_tooltip(sid,index))
    def point_changed(self,item):
        if self._updating or self.active is None:return
        e=self.active;i=item.row()
        try:
            if item.column()==0:e['fit_mask'][i]=item.checkState()==QtCore.Qt.CheckState.Checked
            elif item.column()==7:e['sigma'][i]=float(item.text()) if item.text().strip() else None
            self.request_calculation()
        except ValueError:self.statusBar().showMessage('σ(I) musi być liczbą dodatnią albo pustym polem.',5000)
    def weighting_changed(self,checked):
        if self._updating or self.active is None:return
        self.active['weighted']=checked;self.request_calculation()

    def edit_b(self):
        if self.raw is None:return
        dialog=BDialog(self.raw,self)
        if dialog.exec()==QtWidgets.QDialog.DialogCode.Accepted:
            self.raw=manual_b(self.raw,dialog.values,dialog.reason.text());self.request_calculation()
    def snapshot(self):return build_snapshot(self.dataset,self.results,self.experiment_id,self.dataset.name)

    def show_info(self):
        if self.info_window is None:
            self.info_window = ExperimentInfo(self)
        self.info_window.set_dataset(self.dataset,self.integral_scale)
        self.info_window.show()
        self.info_window.raise_()
    def show_statistics(self):
        if self.dataset is None:return
        if self.statistics_window is None:self.statistics_window=StatisticsWindow(self)
        self.statistics_window.update_live(self.snapshot());self.statistics_window.show();self.statistics_window.raise_()
    def show_map(self):
        if not self.valid:return
        if self.map_window is None:self.map_window=DosyMapWindow(self)
        self.map_window.show()
        self.map_window.update_data(self.experiment_id,self.dataset,self.results)
        self.map_window.raise_()
    def edit_spectrum_settings(self):
        dialog = SpectrumSettingsDialog(self.appearance, self)
        if dialog.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return
        self.appearance = dialog.appearance()
        self.appearance.save(self.settings)
        self.appearance.apply_grid(self.spectrum)
        self.appearance.apply_grid(self.roi_overlay)
        self.render_spectrum()
        self.render_roi_overlay(reset_range=False)

    def show_phase(self):
        if self.dataset is None:return
        if not self.valid or self.worker is not None:
            self.statusBar().showMessage('Poczekaj na zakończenie aktualnego processingu.');return
        if self.phase_window is not None:self.phase_window.raise_();return
        self.phase_window=PhaseWindow(self.dataset,self,self.appearance)
        self.phase_window.applied.connect(self.apply_phase)
        self.phase_window.finished.connect(self.phase_closed)
        self.phase_window.show()
    def phase_closed(self,*args):
        dialog=self.phase_window;self.phase_window=None
        if dialog is not None:dialog.deleteLater()
    def apply_phase(self,values):
        options=self.processing.options();options.update(values)
        self.processing.set_options(options,len(self.raw.fid));self.request_calculation()
    def save_project(self):
        if self.dataset is None or not self.valid:return
        try:
            path=choose_save_file(self,'Zapisz projekt',PROJECTS,'DOSY.proj',PROJECT_FILTER)
            if path:save_project(path,self.dataset,self.entries,integral_scale=self.integral_scale)
        except (OSError,ValueError) as exc:
            QtWidgets.QMessageBox.warning(self,'Zapisz projekt',str(exc))
    def choose_project(self):
        try:
            path,_=QtWidgets.QFileDialog.getOpenFileName(self,'Wczytaj projekt',dialog_directory(PROJECTS),PROJECT_FILTER)
            if path:self.open_path(path,project=True)
        except OSError as exc:
            QtWidgets.QMessageBox.warning(self,'Wczytaj projekt',str(exc))

    def choose_spectrum_png(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Zapisz widmo 1D', 'widmo_1D.png', 'PNG (*.png)')
        if path:
            if not path.lower().endswith('.png'):
                path += '.png'
            try:
                save_light_plot(self.spectrum, path)
            except Exception as exc:
                QtWidgets.QMessageBox.warning(self, 'Zapis widma 1D', str(exc))
    def choose_export(self):
        if not self.valid:return
        parent=QtWidgets.QFileDialog.getExistingDirectory(self,'Folder eksportu')
        if parent:self.export_to(Path(parent)/datetime.now().strftime('DOSY_Lab_%Y%m%d_%H%M%S_%f'))
    def export_to(self,path):
        if not self.valid or self.worker is not None:raise ValueError('Poczekaj na aktualne obliczenia.')
        d=self.dataset;entries=deepcopy(self.results);choices={};revision=self.revision
        scale=self.integral_scale
        if self.statistics_window is not None:choices=self.statistics_window.choices.get(self.experiment_id,{})
        def done(folder):
            for name,plot in [('spectrum',self.spectrum),('decays',self.decay),('residuals',self.residual)]:
                save = save_light_plot if plot is self.spectrum else self.export_plot
                save(plot,Path(folder)/(name+'.png'));save(plot,Path(folder)/(name+'.svg'))
            self.valid=revision==self.revision;self.export_button.setEnabled(self.valid)
            self.statusBar().showMessage('Zapisano eksport w '+str(folder))
        self.start_job(lambda:export_bundle(path,d,entries,choices,integral_scale=scale),done)
    @staticmethod
    def export_plot(plot,path):
        path=Path(path)
        if path.suffix.lower()=='.svg':
            from statistics_svg import export_svg
            export_svg(plot.plotItem,path)
        else:
            ex=pg.exporters.ImageExporter(plot.plotItem);ex.parameters()['width']=1800;ex.export(str(path))
    def closeEvent(self,event):
        self.timer.stop()
        if self.worker is not None and self.worker.isRunning():self.close_pending=True;event.ignore();return
        if self.statistics_window is not None:self.statistics_window.close()
        if self.map_window is not None:self.map_window.close()
        if self.info_window is not None:self.info_window.close()
        super().closeEvent(event)


def main():
    parser=argparse.ArgumentParser(description='DOSY Lab — klasyczne 1D DOSY JEOL')
    parser.add_argument('path',nargs='?',type=Path,help='Plik pomiarowy .jdf lub projekt .proj')
    args=parser.parse_args()
    suffix=args.path.suffix.lower() if args.path else None
    if suffix is not None and suffix not in ('.jdf','.proj'):
        parser.error('Wybierz plik .jdf albo projekt .proj.')
    app=QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    app.setApplicationName('DOSY Lab');app.setFont(QtGui.QFont('Segoe UI',10))
    pg.setConfigOptions(background='#101827',foreground='#dbe4ef',antialias=True)
    window=Window();window.showMaximized()
    if args.path:window.open_path(args.path,project=suffix=='.proj')
    sys.exit(app.exec())
