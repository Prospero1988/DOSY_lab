"""Interactive joint ILT for selected, already fitted DOSY integrals."""
from threading import Event

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from global_ilt import GlobalParameters, solve_global_ilt
from i18n import tr
from ilt import D_UNIT, GRID_SIZES, ILTParameters, solve_ilt
from ilt_help import GLOBAL_ILT_HELP, ILT_WINDOW_FLAGS, add_ilt_info
from ilt_tools import ILTPlotTools, SolveJob
from ilt_regularization import SweepInput, has_systematic_warning, supplied_noise
from ilt_regularization_window import RegularizationControls


COLORS = ('#67e8f9', '#fb923c', '#a7f3d0', '#f9a8d4', '#c4b5fd', '#fde047', '#93c5fd', '#fca5a5')
PRINT_COLORS = ('#155e75', '#9a3412', '#166534', '#9d174d', '#5b21b6', '#854d0e', '#1e40af', '#991b1b')


class GlobalILTWindow(QtWidgets.QDialog, ILTPlotTools):
    solution_changed = QtCore.Signal(object)

    def __init__(self, parent, source, b, entries):
        super().__init__(parent, ILT_WINDOW_FLAGS)
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)
        self.setStyleSheet(parent.styleSheet())
        self.setWindowTitle(tr('DOSY Lab · Global ILT'))
        self.resize(1190, 870)
        self.b = np.array(b, dtype=float, copy=True)
        self.y = np.column_stack([e['result']['corrected_integral'] for e in entries])
        self.fits = tuple(e['fit'] for e in entries)
        self.mask = np.column_stack([fit['fit_mask'] for fit in self.fits]).astype(bool)
        self.ids = [e['id'] for e in entries]
        self.ppm = [(e['bounds'][0]+e['bounds'][1])/2 for e in entries]
        self.context = source + ' · ' + ', '.join(self.ids)
        self.sigma = np.ones_like(self.y)
        self.has_sigma = False
        self.noise_sigma = supplied_noise(entries)
        self.systematic = has_systematic_warning(entries)
        for column, entry in enumerate(entries):
            sigma = np.asarray(entry['sigma'], dtype=float)
            active = sigma[self.mask[:, column]]
            if np.isfinite(active).all() and np.all(active > 0):
                self.sigma[:, column] = sigma
                self.has_sigma = True
        d = [fit['D_m2_s'] for fit in self.fits if fit['D_m2_s'] > 0]
        center = 1/np.ptp(self.b)
        low, high = np.clip([min(d, default=center), max(d, default=center)], 1e-15, 1e-5)
        self.defaults = GlobalParameters(ILTParameters(low/20, high*20, regularization=1e-6, order=0), sparsity=1e-4)
        self.job, self.cancel, self.result, self.independent = None, None, None, ()
        self.revision = 0
        self.closed = False
        self.timer = QtCore.QTimer(self, singleShot=True, interval=180)
        self.timer.timeout.connect(self.start_solve)
        self.export_buttons = []

        root = QtWidgets.QVBoxLayout(self)
        root.addWidget(QtWidgets.QLabel(self.context, textFormat=QtCore.Qt.TextFormat.PlainText))
        note = QtWidgets.QLabel(tr('Migawka bieżącego dopasowania. Otwórz analizę ponownie po zmianie danych.'))
        root.addWidget(note)
        self.build_controls(root)
        tabs = QtWidgets.QTabWidget()
        root.addWidget(tabs, 1)
        distributions = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        tabs.addTab(distributions, tr('Rozkłady'))
        self.aggregate_plot = self.create_plot(distributions, 'A · Rozkład wspólny', 'D / 10⁻¹⁰ [m²/s]', 'P · G / ΣG')
        self.roi_plot = self.create_plot(distributions, 'B · Rozkłady ROI', 'D / 10⁻¹⁰ [m²/s]', '')
        self.roi_plot.addLegend()
        self.roi_plot.export_options = {'series_colors': {sid: PRINT_COLORS[i % len(PRINT_COLORS)] for i,sid in enumerate(self.ids)}}
        distributions.setStretchFactor(0, 1)
        distributions.setStretchFactor(1, 1)

        reconstruction = QtWidgets.QWidget()
        tabs.addTab(reconstruction, tr('Rekonstrukcja'))
        column = QtWidgets.QVBoxLayout(reconstruction)
        self.roi = QtWidgets.QComboBox()
        for sid, ppm in zip(self.ids, self.ppm):
            self.roi.addItem(f'{sid} · {ppm:.4f} ppm')
        self.roi.currentIndexChanged.connect(self.draw_reconstruction)
        column.addWidget(self.roi)
        curves = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        column.addWidget(curves, 1)
        self.decay = self.create_plot(curves, 'C · Zanik', 'b / 10⁹ [s/m²]', 'I / Iₛ')
        self.residual = self.create_plot(curves, 'C · Reszty', 'b / 10⁹ [s/m²]', 'ΔI / Iₛ')
        self.decay.addLegend()
        self.residual.addLegend()
        curves.setStretchFactor(0, 1)
        curves.setStretchFactor(1, 1)

        contribution = QtWidgets.QWidget()
        tabs.addTab(contribution, tr('Udziały'))
        panel = QtWidgets.QVBoxLayout(contribution)
        self.absolute = QtWidgets.QCheckBox(tr('Udziały w oryginalnych jednostkach'))
        self.absolute.toggled.connect(self.draw_heatmap)
        panel.addWidget(self.absolute)
        self.heatmap = self.create_plot(panel, 'D · Udziały składowych', 'ROI', 'D / 10⁻¹⁰ [m²/s]')
        self.heatmap.invertY(True)
        self.heatmap.showGrid(x=False, y=False)
        self.heatmap.format_coordinate = self.heatmap_coordinate
        self.table = QtWidgets.QTableWidget(0, 9)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setHorizontalHeaderLabels(['ROI', 'ppm', 'D mono [m²/s]', tr('D dominujące [m²/s]'),
            tr('Udziały wg rosnącego D [%]'), 'RMS mono / Iₛ', 'RMS Global / Iₛ', 'RMS ILT / Iₛ', 'Status'])
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        tabs.addTab(self.table, tr('Wyniki'))
        self.plots = (self.aggregate_plot, self.roi_plot, self.decay, self.residual, self.heatmap)

        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        root.addWidget(self.summary)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        root.addWidget(self.status)
        root.addWidget(QtWidgets.QLabel(tr('Wspólne maksimum D jest wskazówką. Nie identyfikuje cząsteczki ani stężenia.')))
        self.reset_parameters()

    def build_controls(self, root):
        tuning = QtWidgets.QHBoxLayout()
        root.addLayout(tuning)
        self.penalties = []
        for symbol, lower, upper in [('λ', -120, 20), ('μ', -120, 20)]:
            tuning.addWidget(QtWidgets.QLabel(symbol))
            slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
            slider.setRange(lower, upper)
            edit = QtWidgets.QLineEdit()
            edit.setMaximumWidth(95)
            slider.valueChanged.connect(lambda value, field=edit: self.slider_changed(value, field))
            edit.editingFinished.connect(lambda field=edit, control=slider: self.penalty_changed(field, control))
            self.penalties.append((slider, edit))
            tuning.addWidget(slider, 1)
            tuning.addWidget(edit)
        self.minimum, self.maximum = QtWidgets.QDoubleSpinBox(), QtWidgets.QDoubleSpinBox()
        for label, spin in [('D min', self.minimum), ('D max', self.maximum)]:
            spin.setDecimals(8)
            spin.setRange(1e-7, 1e7)
            spin.setSuffix(' ×10⁻¹⁰ m²/s')
            spin.setKeyboardTracking(False)
            spin.valueChanged.connect(self.queue_solve)
            tuning.addWidget(QtWidgets.QLabel(label))
            tuning.addWidget(spin)
        self.points = QtWidgets.QComboBox()
        for size in GRID_SIZES:
            self.points.addItem(str(size), size)
        self.spacing = QtWidgets.QComboBox()
        self.spacing.addItem(tr('Logarytmiczna'), 'log')
        self.spacing.addItem(tr('Liniowa'), 'linear')
        for label, combo in [('Punkty siatki', self.points), ('Siatka D', self.spacing)]:
            tuning.addWidget(QtWidgets.QLabel(tr(label)))
            tuning.addWidget(combo)
            combo.currentIndexChanged.connect(self.queue_solve)

        options = QtWidgets.QHBoxLayout()
        root.addLayout(options)
        self.iterations = QtWidgets.QSpinBox()
        self.iterations.setRange(1, 100000)
        self.iterations.setKeyboardTracking(False)
        self.iterations.valueChanged.connect(self.queue_solve)
        self.tolerance = QtWidgets.QLineEdit()
        self.tolerance.setMaximumWidth(85)
        self.tolerance.editingFinished.connect(self.queue_solve)
        for label, widget in [('Limit iteracji', self.iterations), ('Tolerancja', self.tolerance)]:
            options.addWidget(QtWidgets.QLabel(tr(label)))
            options.addWidget(widget)
        self.weighted = QtWidgets.QCheckBox(tr('Wagi z SE całek'))
        self.weighted.setEnabled(self.has_sigma)
        self.compare = QtWidgets.QCheckBox(tr('Porównaj niezależne ILT'))
        for button in (self.weighted, self.compare):
            options.addWidget(button)
            button.toggled.connect(self.queue_solve)
        for label, callback in [('Przywróć ustawienia', self.reset_parameters), ('Przelicz ILT', self.recalculate)]:
            button = QtWidgets.QPushButton(tr(label))
            button.clicked.connect(callback)
            options.addWidget(button)
        options.addStretch()
        add_ilt_info(options, self, 'Instrukcja Global ILT', GLOBAL_ILT_HELP)
        combined = QtWidgets.QPushButton(tr('Zapisz wszystkie wykresy PNG'))
        combined.clicked.connect(lambda: self.export_png(self.plots, 'Global_ILT'))
        self.export_buttons.append(combined)
        options.addWidget(combined)
        self.auto_lambda = RegularizationControls(self,options,self.lambda_snapshot,lambda: self.revision,
                                                  self.apply_lambda,self.set_lambda)

        view = QtWidgets.QHBoxLayout()
        root.addLayout(view)
        self.log_axis = QtWidgets.QCheckBox(tr('Logarytmiczna oś D'))
        self.log_axis.setChecked(True)
        self.log_axis.toggled.connect(self.draw_distributions)
        self.common_amplitude = QtWidgets.QCheckBox(tr('Wspólna skala amplitudy ROI'))
        self.common_amplitude.toggled.connect(self.draw_distributions)
        view.addWidget(self.log_axis)
        view.addWidget(self.common_amplitude)
        view.addStretch()
        visibility = QtWidgets.QToolButton(text=tr('Pokaż rozkłady ROI'))
        visibility.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        menu = QtWidgets.QMenu(visibility)
        self.visible_rois = []
        for sid in self.ids:
            action = menu.addAction(sid)
            action.setCheckable(True)
            action.setChecked(True)
            action.toggled.connect(self.draw_distributions)
            self.visible_rois.append(action)
        visibility.setMenu(menu)
        view.addWidget(visibility)

    @staticmethod
    def read_number(widget):
        try:
            return float(widget.text().strip().replace(',', '.'))
        except ValueError as exc:
            raise ValueError('Wpisz prawidłową wartość liczbową.') from exc

    def slider_changed(self, value, edit):
        edit.setText(f'{10.**(value/10):.10g}')
        self.queue_solve()

    def penalty_changed(self, edit, slider):
        try:
            value = self.read_number(edit)
            with QtCore.QSignalBlocker(slider):
                if np.isfinite(value) and value > 0:
                    slider.setValue(round(10*np.log10(value)))
                elif value == 0:
                    slider.setValue(slider.minimum())
        except ValueError:
            pass
        self.queue_solve()

    def reset_parameters(self):
        p = self.defaults
        controls = [self.minimum, self.maximum, self.points, self.spacing, self.iterations, self.weighted, self.compare]
        blockers = [QtCore.QSignalBlocker(control) for control in controls + [s for s,_ in self.penalties]]
        for (slider, edit), value in zip(self.penalties, (p.grid.regularization, p.sparsity)):
            slider.setValue(round(10*np.log10(value)))
            edit.setText(str(value))
        self.minimum.setValue(p.grid.d_min/D_UNIT)
        self.maximum.setValue(p.grid.d_max/D_UNIT)
        self.points.setCurrentIndex(self.points.findData(p.grid.points))
        self.spacing.setCurrentIndex(self.spacing.findData(p.grid.spacing))
        self.iterations.setValue(p.max_iterations)
        self.tolerance.setText(str(p.tolerance))
        self.weighted.setChecked(any(fit['weighted'] for fit in self.fits))
        self.compare.setChecked(False)
        for blocker in blockers:
            blocker.unblock()
        self.queue_solve()

    def recalculate(self):
        self.queue_solve()
        self.timer.stop()
        self.start_solve()

    def queue_solve(self, *_):
        self.auto_lambda.manual_override()
        self.revision += 1
        self.solution_changed.emit(None)
        if self.cancel is not None:
            self.cancel.set()
        for button in self.export_buttons:
            button.setEnabled(False)
        self.set_status(tr('Obliczanie Global ILT…'), stale=True)
        self.timer.start()

    def set_status(self, text, *, stale=False):
        if stale and self.result is not None:
            text += ' ' + tr('Wykresy pokazują poprzednie rozwiązanie.')
        self.status.setText(text)

    def parameters(self):
        regularization, sparsity = [self.read_number(edit) for _,edit in self.penalties]
        grid = ILTParameters(self.minimum.value()*D_UNIT,self.maximum.value()*D_UNIT,
                             self.points.currentData(),self.spacing.currentData(),regularization,0)
        if not 1e-12 <= regularization <= 100:
            raise ValueError('λ musi mieścić się w zakresie 10⁻¹²–100.')
        return GlobalParameters(grid,sparsity,self.iterations.value(),self.read_number(self.tolerance))

    def lambda_snapshot(self):
        source = SweepInput(self.b,self.y,self.mask,self.parameters(),self.sigma if self.weighted.isChecked() else None,
                            self.noise_sigma,systematic=self.systematic)
        ids = tuple(self.ids)
        return self.revision,source.current,lambda: {ids:source},self.context

    def set_lambda(self,value):
        slider,edit = self.penalties[0]
        edit.setText(f'{value:.12g}')
        self.penalty_changed(edit,slider)

    def apply_lambda(self,recommendations):
        previous = self.parameters().grid.regularization
        self.set_lambda(recommendations[tuple(self.ids)].recommendation)
        return previous

    def start_solve(self):
        if self.closed or self.job is not None:
            return
        try:
            params = self.parameters()
        except ValueError as exc:
            self.set_status(tr(str(exc)), stale=True)
            return
        b, y, mask = self.b, self.y, self.mask
        sigma = self.sigma if self.weighted.isChecked() else None
        compare = self.compare.isChecked()
        self.cancel = cancel = Event()

        def compute(p):
            result = solve_global_ilt(b, y, p, include=mask, sigma=sigma, cancelled=cancel.is_set)
            independent = ()
            if compare:
                fits = []
                for index, (signal, active) in enumerate(zip(y.T, mask.T)):
                    if cancel.is_set():
                        break
                    try:
                        fits.append(solve_ilt(b, signal, p.grid, include=active, sigma=None if sigma is None else sigma[:,index]))
                    except (ValueError, RuntimeError, np.linalg.LinAlgError):
                        fits.append(None)
                independent = tuple(fits)
            return result, independent

        self.job = SolveJob(self.revision, params, compute)
        self.job.signals.finished.connect(self.solved)
        QtCore.QThreadPool.globalInstance().start(self.job)

    @QtCore.Slot(int, object, object, str)
    def solved(self, revision, parameters, payload, error):
        self.job, self.cancel = None, None
        if self.closed:
            return
        if revision != self.revision:
            if not self.timer.isActive():
                self.start_solve()
            return
        if error:
            self.set_status(tr(error), stale=True)
            return
        self.result, self.independent = payload
        self.solved_parameters = parameters
        self.draw_distributions()
        self.draw_reconstruction()
        self.draw_heatmap()
        self.fill_table()
        r = self.result
        d = ', '.join(f'{r.diffusion[i]:.4e} m²/s' for i in r.components) or '—'
        self.summary.setText(f"{tr('Maksima Global ILT')}  {d}\nλ = {parameters.grid.regularization:.5g} · μ = {parameters.sparsity:.5g} · "
            f"{tr('Funkcja celu')} = {r.objective:.5g} · {tr('Reszta stacjonarności')} = {r.stationarity:.3g} · {tr('Iteracje')} = {r.iterations}")
        status = 'Solver zakończony' if r.converged else 'Osiągnięto limit iteracji. Wynik przybliżony.'
        self.set_status(' '.join(tr(line) for line in (status, *r.warnings)))
        if r.converged:
            self.solution_changed.emit(r)
        for button in self.export_buttons:
            button.setEnabled(True)

    def draw_distributions(self, *_):
        r = self.result
        if r is None:
            return
        a, o = self.aggregate_plot, self.roi_plot
        log = self.log_axis.isChecked()
        for plot in (a, o):
            plot.clear()
            plot.setLogMode(x=log, y=False)
        d = r.diffusion/D_UNIT
        count = len(d)
        g = r.aggregate
        g = g*count/g.sum() if g.sum() else g
        a.plot(d, g, pen=pg.mkPen('#fb923c', width=2))
        peaks = list(r.components)
        a.plot(d[peaks], g[peaks], pen=None, symbol='o', symbolBrush='#fb923c', symbolSize=8)
        weights = r.weights*r.scales
        common = self.common_amplitude.isChecked()
        totals = weights.sum(axis=0)
        denominator = totals.max() if common else totals
        densities = np.divide(count*weights, denominator, out=np.zeros_like(weights), where=denominator > 0)
        o.setLabel('left', tr('Gęstość na wspólnej skali') if common else tr('Gęstość o polu 1'))
        styles = (QtCore.Qt.PenStyle.SolidLine, QtCore.Qt.PenStyle.DashLine, QtCore.Qt.PenStyle.DotLine)
        for index, (sid, curve, visible) in enumerate(zip(self.ids, densities.T, self.visible_rois)):
            if visible.isChecked():
                style = styles[index // len(COLORS) % len(styles)]
                o.plot(d, curve, pen=pg.mkPen(COLORS[index % len(COLORS)], width=2, style=style), name=sid)
        for peak in peaks:
            o.addItem(pg.InfiniteLine(np.log10(d[peak]) if log else d[peak],
                      pen=pg.mkPen('#64748b', style=QtCore.Qt.PenStyle.DotLine)), ignoreBounds=True)
        for plot in (a, o):
            plot.autoRange()

    def draw_reconstruction(self, *_):
        r = self.result
        if r is None:
            return
        i = self.roi.currentIndex()
        fit, active = self.fits[i], self.mask[:,i]
        scale, observed = r.scales[i], self.y[:,i]
        b, grid = self.b/1e9, fit['b_grid']
        models = [(tr('Monoeksponenta'), fit['curve'], fit['prediction'], '#67e8f9', QtCore.Qt.PenStyle.SolidLine, 'o'),
                  ('Global ILT', r.predict(grid)[:,i], r.predict(self.b)[:,i], '#fb923c', QtCore.Qt.PenStyle.DashLine, 't')]
        independent = self.independent[i] if self.independent else None
        if independent is not None:
            models.append((tr('Niezależne ILT'), independent.predict(grid), independent.predict(self.b),
                           '#a7f3d0', QtCore.Qt.PenStyle.DotLine, 's'))
        self.decay.clear()
        self.residual.clear()
        self.residual.addItem(pg.InfiniteLine(angle=0, pen=pg.mkPen('#64748b')))
        for label, curve, prediction, color, style, symbol in models:
            self.decay.plot(grid/1e9, curve/scale, pen=pg.mkPen(color, width=2, style=style), name=label)
            residual = (observed-prediction)/scale
            self.residual.plot(b[active], residual[active], pen=None, symbol=symbol, symbolBrush=color, name=label)
            if (~active).any():
                self.residual.plot(b[~active], residual[~active], pen=None, symbol='x', symbolPen=color, symbolBrush=color)
        for mask, symbol, label in [(active,'o','Użyte punkty'), (~active,'x','Poza dopasowaniem')]:
            if mask.any():
                self.decay.plot(b[mask], observed[mask]/scale, pen=None, symbol=symbol, symbolBrush='#e2e8f0', name=tr(label))
        self.decay.autoRange()
        self.residual.autoRange()

    def heatmap_values(self):
        return self.result.contributions if self.absolute.isChecked() else 100*self.result.fractions

    def heatmap_coordinate(self, x, y):
        col, row = round(x), round(y)
        if self.result is None or not (0 <= col < len(self.ids) and 0 <= row < len(self.result.components)):
            return ''
        d = self.result.diffusion[self.result.components[row]]
        unit = tr('j.u.·ppm') if self.absolute.isChecked() else '%'
        return f'{self.ids[col]} · D = {d:.4e} m²/s · {self.heatmap_values()[row,col]:.4g} {unit}'

    def draw_heatmap(self, *_):
        r = self.result
        if r is None:
            return
        plot = self.heatmap
        plot.clear()
        values = self.heatmap_values()
        unit = tr('j.u.·ppm') if self.absolute.isChecked() else '%'
        high = max(float(values.max(initial=0.)), 1e-30) if self.absolute.isChecked() else 100.
        plot.setTitle(tr('D · Udziały składowych') + f' [0–{high:.4g} {unit}]')
        ticks = {'bottom_ticks': list(enumerate(self.ids)),
                 'left_ticks': [(i, f'{r.diffusion[peak]/D_UNIT:.3g}') for i,peak in enumerate(r.components)]}
        plot.export_options = ticks
        for axis in ('left', 'bottom'):
            plot.getAxis(axis).setTicks([ticks[axis+'_ticks']])
        if not len(values):
            return
        lut = pg.ColorMap([0., 1.], ['#eff6ff', '#60a5fa']).getLookupTable()
        image = pg.ImageItem(values, axisOrder='row-major', levels=(0., high), lut=lut, autoLevels=False)
        image.setRect(QtCore.QRectF(-.5, -.5, len(self.ids), len(values)))
        plot.addItem(image)
        if values.size <= 100:
            for (row, col), value in np.ndenumerate(values):
                label = pg.TextItem(f'{value:.3g}', color='#101827', anchor=(.5,.5))
                label.setPos(col, row)
                plot.addItem(label, ignoreBounds=True)
        plot.autoRange()

    def fill_table(self):
        r = self.result
        self.table.setRowCount(len(self.ids))
        for i, (sid, ppm, fit, scale, rms) in enumerate(zip(self.ids, self.ppm, self.fits, r.scales, r.rms)):
            fractions = r.fractions[:,i]
            dominant = f'{r.diffusion[r.components[np.argmax(fractions)]]:.4e}' if fractions.any() else '—'
            independent = self.independent[i] if self.independent else None
            values = [sid, f'{ppm:.4f}', f"{fit['D_m2_s']:.4e}", dominant,
                      ' / '.join(f'{100*x:.1f}' for x in fractions) if fractions.any() else '—',
                      f"{fit['RMSE']/scale:.4g}", f'{rms/scale:.4g}', f'{independent.rms/scale:.4g}' if independent else '—',
                      tr('Zbieżny') if r.converged else tr('Wynik przybliżony')]
            for col, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                if col == 4:
                    item.setToolTip(tr('Udziały z sumy wag między minimami rozkładu wspólnego.'))
                self.table.setItem(i, col, item)

    def done(self, result):
        self.auto_lambda.dismiss()
        self.closed = True
        self.timer.stop()
        if self.cancel is not None:
            self.cancel.set()
        super().done(result)
