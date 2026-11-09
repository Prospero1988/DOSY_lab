"""Non-modal PySide6 window for descriptive D statistics and cached fit diagnostics."""
from i18n import tr

from copy import deepcopy
from datetime import datetime
import html
import math
from pathlib import Path
import uuid

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets
import pyqtgraph as pg
from fullscreen_plot import enable_fullscreen
from plot_export import save_light_plot
from storage import source_save_name

from d_statistics import (SCIENTIFIC_NOTE, diagnostic_points, export_statistics,
                          finite, import_statistics, number, summarize, groups, included_records)


class StatisticsWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent, QtCore.Qt.WindowType.Window)
        self.setWindowTitle(tr('DOSY Lab · Statystyka D')); self.resize(1460, 980)
        if parent is not None: self.setStyleSheet(parent.styleSheet())
        self.experiments = {}; self.choices = {}; self.weighted_choices = {}
        self.live_id = None; self.current_id = None; self.active_roi = None
        self.records_by_id = {}; self.row_by_id = {}; self.summary = {}; self.point_rows = []
        self.scale = 1e-10; self._updating = False
        root = QtWidgets.QWidget(); self.setCentralWidget(root)
        layout = QtWidgets.QVBoxLayout(root)
        row = QtWidgets.QHBoxLayout(); layout.addLayout(row)
        row.addWidget(QtWidgets.QLabel(tr('Eksperyment:')))
        self.experiment_combo = QtWidgets.QComboBox(); self.experiment_combo.setMinimumWidth(420)
        row.addWidget(self.experiment_combo, 1); self.experiment_combo.currentIndexChanged.connect(self.experiment_changed)
        self.pin_button = self.button(tr('Zachowaj kopię do porównania'), self.pin_current, row)
        self.button(tr('Wczytaj zapis statystyki'), self.choose_import, row)
        self.source_label = QtWidgets.QLabel(tr('Brak eksperymentu.')); self.source_label.setWordWrap(True); layout.addWidget(self.source_label)
        note = QtWidgets.QLabel(tr(SCIENTIFIC_NOTE)); note.setWordWrap(True)
        note.setStyleSheet('color:#fbbf24;background:#292415;padding:7px;border-radius:5px;'); layout.addWidget(note)
        row = QtWidgets.QHBoxLayout(); layout.addLayout(row)
        self.weighted = QtWidgets.QCheckBox(tr('Średnia ważona wg SE(D)'))
        self.weighted.setToolTip(
            tr('Każde D dostaje wagę 1/SE(D)², więc mniejsze SE daje większy wpływ na wynik. '
            'Włącz, gdy SE poszczególnych dopasowań są wiarygodne i porównywalne. '
            'Bieżące SE z reszt pomijają kalibrację g i korelację między ROI.'))
        self.weighted.toggled.connect(self.weighting_changed); row.addWidget(self.weighted)
        self.two_sd = QtWidgets.QCheckBox(tr('Linie ±2 SD')); self.two_sd.toggled.connect(self.refresh_statistics); row.addWidget(self.two_sd)
        self.show_residuals = QtWidgets.QCheckBox(tr('Pokaż reszty')); self.show_residuals.setChecked(True)
        self.show_residuals.toggled.connect(lambda checked: self.residual_plot.setVisible(checked)); row.addWidget(self.show_residuals)
        row.addStretch()
        self.summary_label = QtWidgets.QLabel('N = 0'); self.summary_label.setWordWrap(True)
        self.summary_label.setStyleSheet('font-size:14px;padding:5px;color:#67e8f9;')
        summary_row = QtWidgets.QHBoxLayout()
        layout.addLayout(summary_row)
        summary_row.addWidget(self.summary_label,1)
        summary_row.addSpacing(18)
        summary_row.addWidget(QtWidgets.QLabel(tr('Związek')))
        self.group_combo = QtWidgets.QComboBox()
        self.group_combo.setToolTip(tr('Wybierz grupę sygnałów, dla której liczone jest podsumowanie D.'))
        self.group_combo.currentIndexChanged.connect(self.group_changed)
        summary_row.addWidget(self.group_combo)
        summary_row.addStretch(1)
        self.weighted_label = QtWidgets.QLabel(''); self.weighted_label.setWordWrap(True); layout.addWidget(self.weighted_label)
        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical); layout.addWidget(splitter, 1)
        upper = QtWidgets.QSplitter(); splitter.addWidget(upper)
        self.table = QtWidgets.QTableWidget(0, 8)
        self.table.setHorizontalHeaderLabels([tr('Uwzględnij'), 'ID', tr('Grupa / związek'), 'ppm', 'D [m²/s]', 'SE(D) [m²/s]', 'R²', 'QC'])
        self.table.horizontalHeaderItem(5).setToolTip(tr('Błąd standardowy oszacowania D z dopasowania zaniku.'))
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemChanged.connect(self.inclusion_changed); self.table.currentCellChanged.connect(self.table_selection_changed)
        upper.addWidget(self.table)
        overview = QtWidgets.QWidget()
        overview_layout = QtWidgets.QVBoxLayout(overview)
        overview_layout.setContentsMargins(0,0,0,0)
        plot_tools = QtWidgets.QHBoxLayout()
        plot_tools.addStretch()
        overview_layout.addLayout(plot_tools)
        self.d_fullscreen = self.button(tr('Pełny ekran D'), lambda:self.d_plot.fullscreen_controller.open(), plot_tools)
        self.d_plot = pg.PlotWidget(); self.d_plot.showGrid(x=True, y=True, alpha=.2)
        self.d_plot.getAxis('left').enableAutoSIPrefix(False)
        overview_layout.addWidget(self.d_plot,1)
        upper.addWidget(overview); upper.setSizes([760, 660])
        self.highlight = None; self.d_positions = {}
        self.d_grid = self.grid_button(plot_tools, self.d_plot)
        self.tabs = QtWidgets.QTabWidget(); splitter.addWidget(self.tabs)
        detail = QtWidgets.QWidget(); detail_layout = QtWidgets.QVBoxLayout(detail)
        self.roi_info = QtWidgets.QLabel(tr('Wybierz wiersz lub punkt D, aby zobaczyć zapisany zanik.'))
        self.roi_info.setTextFormat(QtCore.Qt.TextFormat.PlainText); self.roi_info.setWordWrap(True); detail_layout.addWidget(self.roi_info)
        self.point_info = QtWidgets.QLabel(tr('Kliknij lub najedź na punkt eksperymentalny.')); self.point_info.setWordWrap(True)
        self.point_info.setTextFormat(QtCore.Qt.TextFormat.PlainText); detail_layout.addWidget(self.point_info)
        plot_tools = QtWidgets.QHBoxLayout()
        plot_tools.addStretch()
        detail_layout.addLayout(plot_tools)
        self.decay_fullscreen = self.button(tr('Pełny ekran zaniku'), lambda:self.decay_plot.fullscreen_controller.open(), plot_tools)
        self.residual_fullscreen = self.button(tr('Pełny ekran reszt'), lambda:self.residual_plot.fullscreen_controller.open(), plot_tools)
        self.show_residuals.toggled.connect(self.residual_fullscreen.setEnabled)
        self.decay_plot = pg.PlotWidget(); self.decay_plot.showGrid(x=True, y=True, alpha=.2)
        self.decay_plot.setLabel('left', tr('Całka'), units=tr('j.u.·ppm')); self.decay_plot.setLabel('bottom', 'b / 10⁹ [s/m²]')
        detail_layout.addWidget(self.decay_plot, 3)
        self.residual_plot = pg.PlotWidget(); self.residual_plot.showGrid(x=True, y=True, alpha=.2)
        self.residual_plot.setLabel('left', tr('I − przewidywanie')); self.residual_plot.setLabel('bottom', 'b / 10⁹ [s/m²]')
        self.residual_plot.setXLink(self.decay_plot); detail_layout.addWidget(self.residual_plot, 1)
        self.detail_grid = self.grid_button(plot_tools, self.decay_plot, self.residual_plot)
        self.tabs.addTab(detail, tr('Zanik wybranego ROI'))
        self.comparison = QtWidgets.QTableWidget(0, 7)
        self.comparison.setHorizontalHeaderLabels([tr('Eksperyment (oddzielne grupy)'), 'N', tr('Średnia D [m²/s]'), tr('Mediana D [m²/s]'), 'SD [m²/s]', 'CV [%]', tr('SEM opisowe [m²/s]')])
        self.comparison.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.comparison.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.comparison.horizontalHeader().setStretchLastSection(True)
        self.comparison.cellDoubleClicked.connect(self.comparison_selected)
        self.tabs.addTab(self.comparison, tr('Porównanie eksperymentów — bez łączenia ROI'))
        splitter.setSizes([320, 390])
        row = QtWidgets.QHBoxLayout(); layout.addLayout(row)
        self.csv_button = self.button(tr('Eksport CSV + zapis statystyki'), self.choose_csv_export, row)
        row.addStretch(); row.addWidget(QtWidgets.QLabel(tr('Eksportowany wykres:')))
        self.export_kind = QtWidgets.QComboBox(); self.export_kind.addItems([tr('Zbiorczy D'), tr('Zanik ROI'), tr('Reszty ROI')]); row.addWidget(self.export_kind)
        self.png_button = self.button('PNG', lambda: self.choose_plot_export('png'), row)
        self.svg_button = self.button('SVG', lambda: self.choose_plot_export('svg'), row)
        self.statusBar().showMessage(tr('Zaznaczenie do statystyki nie zmienia ROI, punktów ani dopasowań w głównym oknie.'))
        self.refresh_statistics()
        overview_context = lambda:html.escape(self.summary_label.text()).replace('\n','<br>')
        detail_context = lambda:html.escape(self.roi_info.text()+'\n'+self.point_info.text()).replace('\n','<br>')
        for plot,title,context,reset in (
            (self.d_plot,tr('Statystyka D'),overview_context,self.reset_overview),
            (self.decay_plot,tr('Zanik wybranego ROI'),detail_context,self.reset_diagnostics),
            (self.residual_plot,tr('Reszty wybranego ROI'),detail_context,self.reset_diagnostics),
        ):
            view = plot.getViewBox()
            view.setMouseMode(pg.ViewBox.RectMode)
            view.mouse_help = tr('Przeciągnij prostokąt, aby powiększyć obszar. Dwuklik przywraca cały zakres.')
            enable_fullscreen(plot,title,context,reset_view=reset)

    def reset_overview(self):
        self.d_plot.autoRange(padding=.05)
        count = len(self.d_positions)
        self.d_plot.setXRange(-.6,max(.6,count-.4),padding=0)

    def reset_diagnostics(self):
        self.decay_plot.autoRange(padding=.05)
        self.residual_plot.autoRange(padding=.05)
        self.residual_plot.getViewBox().linkedXChanged()

    @staticmethod
    def button(text, callback, row):
        button = QtWidgets.QPushButton(text); button.clicked.connect(callback); row.addWidget(button); return button

    @staticmethod
    def grid_button(row, *plots):
        button = QtWidgets.QPushButton(tr('Wyłącz siatkę'))
        button.setCheckable(True)
        def toggle(hidden):
            for plot in plots:
                plot.showGrid(x=not hidden, y=not hidden, alpha=.2)
            button.setText(tr('Włącz siatkę') if hidden else tr('Wyłącz siatkę'))
        button.toggled.connect(toggle)
        row.addWidget(button)
        return button

    @property
    def current(self): return self.experiments.get(self.current_id)

    def update_live(self, snapshot):
        if snapshot is None: return
        follow = self.current_id is None or self.current_id == self.live_id
        key = snapshot['id']; group=self.experiments.get(key,{}).get('statistics_group')
        self.experiments[key] = deepcopy(snapshot)
        if group is not None:self.experiments[key]['statistics_group']=group
        self.choices.setdefault(key, {}); self.weighted_choices.setdefault(key, False)
        self.live_id = key
        if follow: self.current_id = key
        self.rebuild_experiment_combo(); self.refresh_current()

    def rebuild_experiment_combo(self):
        self.experiment_combo.blockSignals(True); self.experiment_combo.clear()
        for key, exp in self.experiments.items():
            self.experiment_combo.addItem((tr('Bieżący · ') if key == self.live_id else '')+exp['label'], key)
        index = self.experiment_combo.findData(self.current_id)
        self.experiment_combo.setCurrentIndex(index); self.experiment_combo.blockSignals(False)

    def experiment_changed(self, index):
        if index < 0: return
        self.current_id = self.experiment_combo.itemData(index); self.active_roi = None; self.refresh_current()

    def refresh_current(self):
        exp = self.current
        if exp is None: return
        self.group_combo.blockSignals(True);self.group_combo.clear()
        available=groups(exp)
        self.group_combo.addItem(tr('Wybierz grupę / przypisz ROI'),'')
        for group in available:self.group_combo.addItem(group,group)
        selected=exp.get('statistics_group',available[0] if len(available)==1 else '')
        self.group_combo.setCurrentIndex(max(0,self.group_combo.findData(selected)))
        exp['statistics_group']=self.group_combo.currentData() or ''
        self.group_combo.blockSignals(False)
        self._updating = True
        self.records_by_id = {r['id']: r for r in exp['records']}; self.row_by_id = {}
        self.weighted.blockSignals(True); self.weighted.setChecked(self.weighted_choices.get(self.current_id, False)); self.weighted.blockSignals(False)
        self.table.setRowCount(len(exp['records']))
        for row, record in enumerate(exp['records']):
            sid = record['id']; self.row_by_id[sid] = row
            check = QtWidgets.QTableWidgetItem(); check.setData(QtCore.Qt.ItemDataRole.UserRole, sid)
            check.setFlags(QtCore.Qt.ItemFlag.ItemIsSelectable | QtCore.Qt.ItemFlag.ItemIsUserCheckable | (QtCore.Qt.ItemFlag.ItemIsEnabled if record['eligible'] else QtCore.Qt.ItemFlag.NoItemFlags))
            check.setCheckState(QtCore.Qt.CheckState.Checked if record['eligible'] and self.choices[self.current_id].get(sid, True) else QtCore.Qt.CheckState.Unchecked)
            check.setToolTip(tr('Udział w statystyce; nie usuwa sygnału.') if record['eligible'] else tr(record['qc_detail']))
            self.table.setItem(row, 0, check)
            vals = [sid, record.get('group',''), self.ppm(record['ppm']), number(record['D']), number(record['SE']),
                    f"{record['R2']:.5f}" if record['R2'] is not None else '—', tr(record['qc'])]
            for col, value in enumerate(vals, 1):
                item = QtWidgets.QTableWidgetItem(value); item.setToolTip(tr(record['qc_detail']))
                if col == 1 and record['name'] != sid:
                    item.setToolTip(f"{record['name']}\n{tr(record['qc_detail'])}")
                self.table.setItem(row, col, item)
        self.source_label.setText(tr(f"{exp['source_name']} · {exp.get('experiment',{}).get('sequence','DOSY')} · SHA256: {exp.get('source_sha256') or 'brak'} · oddzielna seria wyników"))
        self._updating = False
        if self.active_roi not in self.records_by_id: self.active_roi = next(iter(self.records_by_id), None)
        self.refresh_statistics()
        self.select_roi(self.active_roi)

    @staticmethod
    def ppm(v): return f'{v:.2f}' if finite(v) is not None else '—'

    def inclusion_changed(self, item):
        if self._updating or item.column() != 0 or not self.current: return
        sid = item.data(QtCore.Qt.ItemDataRole.UserRole)
        self.choices[self.current_id][sid] = item.checkState() == QtCore.Qt.CheckState.Checked
        self.refresh_statistics()

    def weighting_changed(self, checked):
        if self.current_id is not None: self.weighted_choices[self.current_id] = checked
        self.refresh_statistics()

    def refresh_statistics(self, *args):
        exp = self.current
        self.d_plot.clear(); self.d_positions = {}; self.highlight = None
        if exp is None:
            self.summary_label.setText(tr('N = 0 · Brak wyników. Otwórz dane i zaznacz ROI w głównym oknie.'))
            self.csv_button.setEnabled(False); return
        self.csv_button.setEnabled(True)
        records = included_records(exp, self.choices[self.current_id])
        self.d_positions = {r['id']: i for i, r in enumerate(records)}
        self.summary = summarize(exp, self.choices[self.current_id], self.weighted.isChecked())
        s = self.summary
        cv = f"{s['CV_percent']:.3f}%" if s['CV_percent'] is not None else '—'
        self.summary_label.setText(tr(f"N = {s['N']}    Średnia: {number(s['mean'])}    Mediana: {number(s['median'])}    SD (ddof=1): {number(s['SD'])} m²/s\n"
                                   f"CV: {cv}    SEM opisowe: {number(s['SEM'])}    Min: {number(s['minimum'])}    Max: {number(s['maximum'])} m²/s"))
        self.weighted_label.setText((tr(f"Średnia ważona: {number(s['weighted_mean'])} m²/s. ") if self.weighted.isChecked() else '')+tr(s['weighted_status']))
        values = [abs(r['D']) for r in records]
        values += [r['SE'] for r in records if finite(r['SE']) is not None and r['SE'] >= 0]
        largest = max(values, default=0)
        exponent = max(-300, min(300, math.floor(math.log10(largest)))) if largest else -10
        self.scale = 10.**exponent
        self.d_plot.setLabel('left', f'D / 10^{exponent} [m²/s]'); self.d_plot.setLabel('bottom', tr('ID sygnału'))
        legend=tr('● uwzględnione · słupki ±SE(D) · linia średniej · pasmo ±1 SD')
        if not records: legend=tr('Brak zaznaczonych sygnałów w wybranej grupie.')
        if s['SD'] is None: legend+=tr('<br>SD i SEM wymagają co najmniej dwóch sygnałów')
        if self.two_sd.isChecked(): legend+=tr('<br>Przerywane: ±2 SD')
        if s['weighted_mean'] is not None: legend+=tr(' · linia kropkowana: średnia ważona')
        self.d_plot.setTitle(html.escape(exp['label'])+'<br><span style="font-size:9pt">'+legend+'</span>')
        # PlotItem defaults to a fixed 30 px title row, even for a multiline legend.
        title = self.d_plot.plotItem.titleLabel
        height = max(30,title.item.boundingRect().height())
        title.setMaximumHeight(height)
        self.d_plot.plotItem.layout.setRowFixedHeight(0,height)
        self.d_plot.getAxis('bottom').setTicks([self.d_ticks()])
        if s['mean'] is not None:
            mean = s['mean']/self.scale
            if s['SD'] is not None:
                sd = s['SD']/self.scale
                band = pg.LinearRegionItem((mean-sd, mean+sd), orientation=pg.LinearRegionItem.Horizontal,
                                           movable=False, brush=pg.mkBrush(103,232,249,28), pen=pg.mkPen(None))
                band.setZValue(-20); band.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton); self.d_plot.addItem(band)
                if self.two_sd.isChecked():
                    for y in (mean-2*sd, mean+2*sd): self.d_plot.addItem(pg.InfiniteLine(y, angle=0, pen=pg.mkPen('#64748b', style=QtCore.Qt.PenStyle.DashLine)))
            self.d_plot.addItem(pg.InfiniteLine(mean, angle=0, pen=pg.mkPen('#67e8f9', width=1.5)))
        if s['weighted_mean'] is not None:
            self.d_plot.addItem(pg.InfiniteLine(s['weighted_mean']/self.scale, angle=0, pen=pg.mkPen('#c4b5fd', style=QtCore.Qt.PenStyle.DotLine)))
        spots = []; ex=[]; ey=[]; es=[]
        for i, r in enumerate(records):
            spots.append({'pos': (i, r['D']/self.scale), 'data': r['id'], 'symbol': 'o',
                          'brush': pg.mkBrush(r['color']), 'pen': pg.mkPen(r['color']), 'size': 10})
            if r['SE'] is not None and r['SE'] >= 0:
                ex.append(i); ey.append(r['D']/self.scale); es.append(r['SE']/self.scale)
        if ex:
            errors = pg.ErrorBarItem(x=np.array(ex), y=np.array(ey), top=np.array(es), bottom=np.array(es), beam=.2, pen=pg.mkPen('#94a3b8'))
            self.d_plot.addItem(errors)
        self.scatter = pg.ScatterPlotItem(spots=spots, hoverable=True, tip=self.d_tooltip, useCache=False)
        self.scatter.sigClicked.connect(self.scatter_clicked); self.d_plot.addItem(self.scatter)
        self.highlight = pg.ScatterPlotItem(size=18, pen=pg.mkPen('#ffffff', width=2), brush=pg.mkBrush(None), useCache=False)
        self.d_plot.addItem(self.highlight)
        self.d_plot.setXRange(-.6, max(.6, len(records)-.4), padding=0)
        self.d_plot.enableAutoRange(axis=pg.ViewBox.YAxis, enable=True)
        self.draw_highlight(); self.refresh_comparison()

    def d_ticks(self):
        return [(x, sid) for sid, x in self.d_positions.items()]

    def d_tooltip(self, x, y, data):
        r = self.records_by_id.get(data)
        return f"{data} · D={number(r['D'])} m²/s · SE={number(r['SE'])}\n{tr(r['qc'])}: {tr(r['qc_detail'])}" if r else str(data)

    def scatter_clicked(self, item, points, event):
        if len(points): self.select_roi(points[0].data())

    def table_selection_changed(self, row, col, old_row, old_col):
        if not self._updating and row >= 0:
            item = self.table.item(row, 1)
            if item: self.select_roi(item.text())

    def select_roi(self, sid):
        self.active_roi = sid
        row = self.row_by_id.get(sid)
        if row is not None and self.table.currentRow() != row:
            self._updating = True; self.table.setCurrentCell(row, 1); self._updating = False
        self.draw_highlight(); self.refresh_diagnostic()

    def draw_highlight(self):
        if self.highlight is None: return
        r = self.records_by_id.get(self.active_roi)
        if r is not None and r['id'] in self.d_positions:
            self.highlight.setData([self.d_positions[r['id']]], [r['D']/self.scale])
        else: self.highlight.setData([], [])

    def refresh_diagnostic(self):
        self.decay_plot.clear(); self.residual_plot.clear(); self.point_rows = []
        self.residual_plot.addItem(pg.InfiniteLine(0, angle=0, pen=pg.mkPen('#64748b')))
        self.point_info.setText(tr('Kliknij lub najedź na punkt: g, b, całka, przewidywanie, reszta i S/N, jeśli dostępne.'))
        record = self.records_by_id.get(self.active_roi)
        if record is None:
            self.roi_info.setText(tr('Brak sygnałów. Wybierz ROI w głównym programie.')); return
        fit, result = record['fit'], record['result']
        self.point_rows, warning = diagnostic_points(self.current, record)
        bg = tr('Wielomianowa korekcja bazy') if result.get('baseline_applied') else tr('Bez korekcji bazy')
        model = fit.get('model', tr('brak modelu'))
        r2 = f"{record['R2']:.5f}" if record['R2'] is not None else '—'
        info = (f"{record['id']} · {record['name']} · {self.ppm(record['ppm'])} ppm · grupa: {record.get('group') or tr('nieprzypisana')}\n"
                f"D {number(record['D'])} m²/s · SE(D) {number(record['SE'])} m²/s · R² {r2} · {model}\n{bg}")
        if record['qc'] != 'OK': info += '\n'+tr(record['qc_detail'])
        if self.point_rows:
            used=sum(p['fit_included'] is True for p in self.point_rows)
            excluded=sum(p['fit_included'] is False for p in self.point_rows)
            unknown=sum(p['fit_included'] is None for p in self.point_rows)
            info+=tr(f'\nPunkty: {used} w dopasowaniu · {excluded} poza fitem · {unknown} o nieustalonym statusie.')
        if warning: info += '\n'+tr(warning)
        self.roi_info.setText(tr(info))
        self.decay_plot.setTitle(html.escape(record['id'])+tr(' · rzeczywiste I(b) i zapisany model · × poza fitem · □ status nieznany'))
        self.measurement_scatters = []
        for status, symbol, color in [(True,'o',record['color']),(False,'x','#fb923c'),(None,'s','#94a3b8')]:
            spots = [{'pos': (r['b']/1e9, r['integral']), 'data': r['index']} for r in self.point_rows
                     if r['fit_included'] is status and r['b'] is not None and r['integral'] is not None]
            if not spots: continue
            scatter = pg.ScatterPlotItem(spots=spots, symbol=symbol, size=9, pen=pg.mkPen(color), brush=pg.mkBrush(color), hoverable=True, tip=self.measurement_tooltip, useCache=False)
            scatter.sigClicked.connect(self.measurement_clicked); scatter.sigHovered.connect(self.measurement_hovered)
            self.decay_plot.addItem(scatter); self.measurement_scatters.append(scatter)
        if fit.get('ok'):
            try:
                bx = np.asarray(fit.get('b_grid'), dtype=float); curve = np.asarray(fit.get('curve'), dtype=float)
                if bx.ndim == 1 and bx.shape == curve.shape and len(bx):
                    valid = np.isfinite(bx)&np.isfinite(curve)
                    self.decay_plot.plot(bx[valid]/1e9, curve[valid], pen=pg.mkPen(record['color'], width=2))
                else:
                    predicted = [r for r in self.point_rows if r['b'] is not None and r['prediction'] is not None]
                    self.decay_plot.plot([r['b']/1e9 for r in predicted], [r['prediction'] for r in predicted], pen=pg.mkPen(record['color']))
            except (ValueError, TypeError): pass
        for status, symbol, color in [(True,'o',record['color']),(False,'x','#fb923c'),(None,'s','#94a3b8')]:
            rows = [r for r in self.point_rows if r['fit_included'] is status and r['b'] is not None and r['residual'] is not None]
            self.residual_plot.plot([r['b']/1e9 for r in rows], [r['residual'] for r in rows], pen=None, symbol=symbol, symbolBrush=color, symbolPen=None, symbolSize=7, useCache=False)
        self.measurement_highlight = pg.ScatterPlotItem(size=15, pen=pg.mkPen('#ffffff',width=2), brush=pg.mkBrush(None), useCache=False)
        self.decay_plot.addItem(self.measurement_highlight)

    def measurement_tooltip(self, x, y, data): return self.point_description(data)

    def point_description(self, index):
        r = next((r for r in self.point_rows if r['index'] == index), None)
        if r is None: return tr('Brak danych punktu.')
        state = tr('w dopasowaniu') if r['fit_included'] else (tr('poza fitem') if r['fit_included'] is False else tr('udział w fit nieustalony'))
        return (tr(f"Punkt {index+1} · {state} · g={number(r['g'])} T/m · b={number(r['b'])} s/m²\n"
                f"Całka={number(r['integral'])} · przewidywanie={number(r['prediction'])} · reszta={number(r['residual'])} · S/N={number(r['snr']) if r['snr'] is not None else tr('nie oszacowano')}"))

    def measurement_clicked(self, item, points, event):
        if len(points): self.show_point(points[0].data())

    def measurement_hovered(self, item, points, event):
        if len(points): self.point_info.setText(self.point_description(points[0].data()))

    def show_point(self, index):
        self.point_info.setText(self.point_description(index))
        r = next((r for r in self.point_rows if r['index'] == index), None)
        if r and r['b'] is not None and r['integral'] is not None:
            self.measurement_highlight.setData([r['b']/1e9], [r['integral']])

    def group_changed(self,index):
        if self.current is not None:
            self.current['statistics_group']=self.group_combo.currentData() or ''
            self.refresh_statistics()

    def refresh_comparison(self):
        rows=[(key,e,g) for key,e in self.experiments.items() for g in groups(e)]
        self.comparison.setRowCount(len(rows))
        for i,(key,e,g) in enumerate(rows):
            s=summarize(e,self.choices[key],group=g)
            values=[e['label']+' · '+g,str(s['N']),number(s['mean']),number(s['median']),number(s['SD']),
                    number(s['CV_percent']),number(s['SEM'])]
            for j,value in enumerate(values):
                item=QtWidgets.QTableWidgetItem(value);item.setData(QtCore.Qt.ItemDataRole.UserRole,(key,g));self.comparison.setItem(i,j,item)

    def comparison_selected(self,row,column):
        item=self.comparison.item(row,column)
        if item:
            key,group=item.data(QtCore.Qt.ItemDataRole.UserRole)
            self.experiment_combo.setCurrentIndex(self.experiment_combo.findData(key))
            self.group_combo.setCurrentIndex(self.group_combo.findData(group))

    def pin_current(self):
        if self.current is None: return
        key = uuid.uuid4().hex; exp = deepcopy(self.current); exp['id'] = key; exp['label'] += tr(f' · kopia {len(self.experiments)+1}')
        self.experiments[key] = exp; self.choices[key] = dict(self.choices[self.current_id]); self.weighted_choices[key] = self.weighted.isChecked()
        self.rebuild_experiment_combo(); self.refresh_comparison()
        self.statusBar().showMessage(tr('Zachowano niezależną kopię. Kolejne fity w głównym oknie jej nie zmienią.'))

    def load_snapshot(self, path):
        exp, choices = import_statistics(path); key = exp['id']
        self.experiments[key] = exp; self.choices[key] = choices; self.weighted_choices[key] = False
        self.current_id = key; self.active_roi = None; self.rebuild_experiment_combo(); self.refresh_current()

    def choose_import(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, tr('Wczytaj statistics_snapshot.json z DOSY Lab'), '', tr('Zapis statystyki (*.json)'))
        if path:
            try: self.load_snapshot(path)
            except Exception as e: QtWidgets.QMessageBox.warning(self, tr('Odczyt statystyki'), str(e))

    def export_csv(self, directory):
        if self.current is None: raise ValueError(tr('Brak eksperymentu.'))
        return export_statistics(directory, self.current, self.choices[self.current_id], self.weighted.isChecked(),
                                 [(exp, self.choices[key]) for key, exp in self.experiments.items()], self.active_roi)

    def choose_csv_export(self):
        parent = QtWidgets.QFileDialog.getExistingDirectory(self, tr('Wybierz folder na eksport Statystyki D'))
        if parent:
            try:
                folder = Path(parent)/datetime.now().strftime('statystyka_D_%Y%m%d_%H%M%S_%f')
                self.export_csv(folder); self.statusBar().showMessage(tr(f'Zapisano CSV i kopię wyników: {folder}'))
            except Exception as e: QtWidgets.QMessageBox.warning(self, tr('Eksport statystyki'), str(e))

    def export_plot(self, path, kind=None):
        kind = self.export_kind.currentIndex() if kind is None else kind
        plot = [self.d_plot, self.decay_plot, self.residual_plot][kind]
        save_light_plot(plot, path, bottom_ticks=self.d_ticks() if kind == 0 else None)

    def choose_plot_export(self, suffix):
        source_name = (self.current or {}).get('source_name')
        filename = source_save_name(source_name, '.' + suffix, 'statystyka_D')
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, tr('Zapisz wykres'), filename, tr(f'Wykres (*.{suffix})'))
        if path:
            try:
                target = Path(path)
                if target.suffix.lower() not in ('.png', '.svg'): target = target.with_name(target.name + '.' + suffix)
                self.export_plot(target); self.statusBar().showMessage(tr(f'Zapisano wykres: {target}'))
            except Exception as e: QtWidgets.QMessageBox.warning(self, tr('Eksport wykresu'), str(e))
