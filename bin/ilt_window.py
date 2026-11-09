"""On-demand ILT diagnostics for a snapshot of one DOSY integral."""
import html
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from d_statistics import number
from fullscreen_plot import enable_fullscreen
from i18n import tr
from ilt import D_UNIT, GRID_SIZES, ILTParameters, solve_ilt
from plot_export import save_light_plot


class SolveSignals(QtCore.QObject):
    finished = QtCore.Signal(int, object, object, str)


class SolveJob(QtCore.QRunnable):
    def __init__(self, revision, parameters, b, intensity, mask, sigma):
        super().__init__()
        self.revision, self.parameters = revision, parameters
        self.data = b, intensity, mask, sigma
        self.signals = SolveSignals()

    def run(self):
        b, y, mask, sigma = self.data
        result, error = None, ''
        try:
            result = solve_ilt(b, y, self.parameters, include=mask, sigma=sigma)
        except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as exc:
            error = str(exc)
        self.signals.finished.emit(self.revision, self.parameters, result, error)


class ILTWindow(QtWidgets.QDialog):
    def __init__(self, parent, source, b, entry):
        super().__init__(parent, QtCore.Qt.WindowType.Window | QtCore.Qt.WindowType.WindowMaximizeButtonHint)
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)
        self.setWindowTitle(tr('DOSY Lab · Eksperymentalne ILT'))
        self.setStyleSheet(parent.styleSheet())
        self.resize(1150, 920)
        self.b = np.array(b, dtype=float, copy=True)
        self.y = np.array(entry['result']['corrected_integral'], copy=True)
        self.fit = entry['fit']  # the accepted fit is a snapshot, read only
        self.mask = np.array(self.fit['fit_mask'], dtype=bool, copy=True)
        self.sigma = np.array(entry['sigma'], dtype=float) if self.fit['weighted'] else None
        self.scale = float(np.max(np.abs(self.y[self.mask])))
        center = self.fit['D_m2_s'] or 1 / np.ptp(self.b[self.mask])
        center = float(np.clip(center, 1e-15, 1e-5))
        self.defaults = ILTParameters(center / 20, center * 20)
        self.result, self.solved_parameters, self.job = None, None, None
        self.revision = 0
        self.closed = False
        self.timer = QtCore.QTimer(self, singleShot=True, interval=160)
        self.timer.timeout.connect(self.start_solve)

        root = QtWidgets.QVBoxLayout(self)
        lo, hi = entry['bounds']
        self.context = f"{source} · {entry['id']} · {(lo+hi)/2:.4f} ppm [{lo:.4f}–{hi:.4f}] · {entry['group']}"
        root.addWidget(QtWidgets.QLabel(self.context, textFormat=QtCore.Qt.TextFormat.PlainText))
        root.addWidget(QtWidgets.QLabel(tr('Migawka bieżącego dopasowania. Otwórz analizę ponownie po zmianie danych.')))
        self.create_controls(root)
        warning = QtWidgets.QLabel(tr('Maksima zależą od danych, siatki i λ. Samo maksimum nie potwierdza osobnej populacji, a jego szerokość nie jest niepewnością D.'))
        warning.setWordWrap(True)
        root.addWidget(warning)
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        root.addWidget(self.summary)

        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        root.addWidget(splitter, 1)
        self.attenuation = self.create_plot(splitter, 'A · Zanik sygnału', 'b / 10⁹ [s/m²]', 'I / Iₛ')
        self.distribution = self.create_plot(splitter, 'B · Rozkład ILT', 'D / 10⁻¹⁰ [m²/s]', '')
        self.residual = self.create_plot(splitter, 'C · Reszty', 'b / 10⁹ [s/m²]', 'ΔI / Iₛ')
        self.attenuation.addLegend()
        self.residual.addLegend()
        self.plot_list = (self.attenuation, self.distribution, self.residual)
        for index, plot in enumerate(self.plot_list):
            splitter.setStretchFactor(index, 1)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        root.addWidget(self.status)

        self.render_comparison()
        self.reset_parameters()

    def create_controls(self, root):
        row = QtWidgets.QHBoxLayout()
        root.addLayout(row)
        row.addWidget(QtWidgets.QLabel('λ'))
        self.slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.slider.setRange(-120, 20)
        self.slider.setMinimumWidth(100)
        self.slider.valueChanged.connect(self.slider_changed)
        row.addWidget(self.slider, 1)
        self.lambda_edit = QtWidgets.QLineEdit()
        self.lambda_edit.setMaximumWidth(90)
        self.lambda_edit.editingFinished.connect(self.lambda_changed)
        row.addWidget(self.lambda_edit)
        self.minimum, self.maximum = QtWidgets.QDoubleSpinBox(), QtWidgets.QDoubleSpinBox()
        for label, spin in (('D min', self.minimum), ('D max', self.maximum)):
            row.addWidget(QtWidgets.QLabel(label))
            spin.setDecimals(8)
            spin.setRange(1e-7, 1e7)
            spin.setSuffix(' ×10⁻¹⁰ m²/s')
            spin.setKeyboardTracking(False)
            spin.valueChanged.connect(self.queue_solve)
            row.addWidget(spin)
        for label, attr, items in (
            ('Punkty siatki', 'points', [(str(n), n) for n in GRID_SIZES]),
            ('Siatka D', 'spacing', [(tr('Logarytmiczna'), 'log'), (tr('Liniowa'), 'linear')]),
            ('Kara', 'order', [(tr('Druga pochodna'), 2), (tr('Rząd zerowy'), 0)]),
        ):
            row.addWidget(QtWidgets.QLabel(tr(label)))
            combo = QtWidgets.QComboBox()
            for text, value in items:
                combo.addItem(text, value)
            combo.currentIndexChanged.connect(self.queue_solve)
            row.addWidget(combo)
            setattr(self, attr, combo)
        actions = QtWidgets.QHBoxLayout()
        root.addLayout(actions)
        self.log_axis = QtWidgets.QCheckBox(tr('Logarytmiczna oś D'))
        self.log_axis.toggled.connect(self.render_distribution)
        actions.addWidget(self.log_axis)
        for label, callback in [('Przywróć ustawienia', self.reset_parameters), ('Przelicz ILT', self.recalculate)]:
            button = QtWidgets.QPushButton(tr(label))
            button.clicked.connect(callback)
            actions.addWidget(button)
        actions.addStretch()
        self.combined_export = QtWidgets.QPushButton(tr('Zapisz wszystkie wykresy PNG'))
        self.combined_export.clicked.connect(lambda: self.export_png(self.plot_list, 'ILT'))
        actions.addWidget(self.combined_export)
        self.export_buttons = [self.combined_export]

    def create_plot(self, splitter, title, bottom, left):
        panel = QtWidgets.QWidget()
        splitter.addWidget(panel)
        layout = QtWidgets.QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        bar = QtWidgets.QHBoxLayout()
        layout.addLayout(bar)
        coordinate = QtWidgets.QLabel()
        bar.addWidget(coordinate, 1)
        plot = pg.PlotWidget()
        layout.addWidget(plot, 1)
        plot.setTitle(tr(title))
        plot.setLabel('bottom', bottom)
        plot.setLabel('left', left)
        plot.showGrid(x=True, y=True, alpha=.2)
        view = plot.getViewBox()
        view.setMouseMode(pg.ViewBox.RectMode)
        enable_fullscreen(plot, tr(title), lambda: html.escape(self.context), reset_view=plot.autoRange)
        for label, callback in [('Cały zakres', plot.autoRange), ('Pełny ekran', plot.fullscreen_controller.open)]:
            button = QtWidgets.QPushButton(tr(label))
            button.clicked.connect(callback)
            bar.addWidget(button)
        pan = QtWidgets.QPushButton(tr('Przesuwaj'))
        pan.setCheckable(True)
        pan.toggled.connect(lambda checked: view.setMouseMode(pg.ViewBox.PanMode if checked else pg.ViewBox.RectMode))
        bar.addWidget(pan)
        save = QtWidgets.QPushButton(tr('Zapisz PNG'))
        save.clicked.connect(lambda: self.export_png((plot,), 'ILT_' + title[0]))
        bar.addWidget(save)
        self.export_buttons.append(save)

        def read_coordinate(pos):
            if view.sceneBoundingRect().contains(pos):
                point = view.mapSceneToView(pos)
                x = 10.**point.x() if plot.getAxis('bottom').logMode and abs(point.x()) < 300 else point.x()
                coordinate.setText(f'x = {x:.5g}   y = {point.y():.5g}')

        plot.scene().sigMouseMoved.connect(read_coordinate)
        return plot

    def reset_parameters(self):
        p = self.defaults
        with (QtCore.QSignalBlocker(self.slider), QtCore.QSignalBlocker(self.minimum),
              QtCore.QSignalBlocker(self.maximum), QtCore.QSignalBlocker(self.points),
              QtCore.QSignalBlocker(self.spacing), QtCore.QSignalBlocker(self.order)):
            self.slider.setValue(round(10*np.log10(p.regularization)))
            self.lambda_edit.setText(str(p.regularization))
            self.minimum.setValue(p.d_min / D_UNIT)
            self.maximum.setValue(p.d_max / D_UNIT)
            for combo, value in ((self.points, p.points), (self.spacing, p.spacing), (self.order, p.order)):
                combo.setCurrentIndex(combo.findData(value))
        self.log_axis.setChecked(True)
        self.queue_solve()

    def slider_changed(self, step):
        self.lambda_edit.setText(f'{10.**(step/10):.10g}')
        self.queue_solve()

    def lambda_value(self):
        try:
            value = float(self.lambda_edit.text().strip().replace(',', '.'))
        except ValueError:
            value = np.nan
        if not 1e-12 <= value <= 100:
            raise ValueError('λ musi mieścić się w zakresie 10⁻¹²–100.')
        return value

    def lambda_changed(self):
        try:
            with QtCore.QSignalBlocker(self.slider):
                self.slider.setValue(round(10*np.log10(self.lambda_value())))
        except ValueError:
            pass
        self.queue_solve()

    def recalculate(self):
        self.queue_solve()
        self.timer.stop()
        self.start_solve()

    def queue_solve(self, *_):
        self.revision += 1
        for button in self.export_buttons:
            button.setEnabled(False)
        self.show_error('Oczekiwanie na rozwiązanie ILT…')
        self.timer.start()

    def start_solve(self):
        if self.closed or self.job is not None:
            return
        try:
            p = ILTParameters(self.minimum.value()*D_UNIT, self.maximum.value()*D_UNIT,
                              self.points.currentData(), self.spacing.currentData(), self.lambda_value(), self.order.currentData())
        except ValueError as exc:
            self.show_error(str(exc))
            return
        self.job = SolveJob(self.revision, p, self.b, self.y, self.mask, self.sigma)
        self.job.signals.finished.connect(self.solved)
        QtCore.QThreadPool.globalInstance().start(self.job)

    @QtCore.Slot(int, object, object, str)
    def solved(self, revision, parameters, result, error):
        self.job = None
        if self.closed:
            return
        if revision != self.revision:
            if not self.timer.isActive():
                self.start_solve()
            return
        if error:
            self.show_error(error)
            return
        self.result, self.solved_parameters = result, parameters
        self.render_comparison()
        self.render_distribution()
        self.status.setText(' '.join(map(tr, result.warnings)))
        for button in self.export_buttons:
            button.setEnabled(True)

    def show_error(self, error):
        suffix = ' ' + tr('Wykresy pokazują poprzednie rozwiązanie.') if self.result is not None else ''
        self.status.setText(tr(error) + suffix)

    def render_comparison(self):
        a, r = self.attenuation, self.residual
        for plot in (a, r):
            plot.clear()
        b = self.b/1e9
        grid = self.fit['b_grid']
        mono = '#67e8f9'
        ilt = '#fb923c'
        a.plot(grid/1e9, self.fit['curve']/self.scale, pen=pg.mkPen(mono, width=2), name=tr('Monoeksponenta'))
        predicted = self.result.predict(self.b) if self.result is not None else None
        if predicted is not None:
            a.plot(grid/1e9, self.result.predict(grid)/self.scale,
                   pen=pg.mkPen(ilt, width=2, style=QtCore.Qt.PenStyle.DashLine), name='ILT')
        r.addItem(pg.InfiniteLine(angle=0, pen=pg.mkPen('#64748b')))
        for mask, symbol, label in ((self.mask, 'o', tr('Użyte punkty')), (~self.mask, 'x', tr('Poza dopasowaniem'))):
            if mask.any():
                a.plot(b[mask], self.y[mask]/self.scale, pen=None, symbol=symbol, symbolSize=8, symbolBrush='#e2e8f0', name=label)
        for residual, color, symbol, omitted, label in (
            (self.fit['residuals'], mono, 'o', 'x', tr('Monoeksponenta')),
            (self.y-predicted if predicted is not None else None, ilt, 't', '+', 'ILT'),
        ):
            if residual is None:
                continue
            r.plot(b[self.mask], residual[self.mask]/self.scale, pen=None, symbol=symbol, symbolSize=8, symbolBrush=color, name=label)
            if (~self.mask).any():
                r.plot(b[~self.mask], residual[~self.mask]/self.scale, pen=None, symbol=omitted, symbolSize=8,
                       symbolBrush=color, name=f"{label} · {tr('poza fitem')}")
        self.refresh_summary()

    def render_distribution(self, *_):
        plot = self.distribution
        plot.clear()
        log_axis = self.log_axis.isChecked()
        plot.setLogMode(x=log_axis, y=False)
        if self.result is None:
            return
        result, p = self.result, self.solved_parameters
        x = result.diffusion/D_UNIT
        y = result.density
        plot.setLabel('left', tr('Gęstość względem u') + (' [u ∝ log(D)]' if p.spacing == 'log' else ' [u ∝ D]'))
        plot.plot(x, y, pen=pg.mkPen('#fb923c', width=2))
        peaks = list(result.peaks)
        plot.plot(x[peaks], y[peaks], pen=None, symbol='o', symbolSize=7, symbolBrush='#fb923c')
        d = self.fit['D_m2_s']/D_UNIT
        if d > 0 and p.d_min <= self.fit['D_m2_s'] <= p.d_max:
            line = pg.InfiniteLine(np.log10(d) if log_axis else d, pen=pg.mkPen('#67e8f9', style=QtCore.Qt.PenStyle.DashLine),
                                   label=tr('D mono'), labelOpts={'position': .92})
            plot.addItem(line, ignoreBounds=True)
        plot.autoRange()

    def refresh_summary(self):
        f = self.fit
        text = (f"D mono = {number(f['D_m2_s'])} m²/s · SE = {number(f['D_standard_error_m2_s'])} m²/s · "
                f"RMS mono = {f['RMSE']/self.scale:.4g} Iₛ · Iₛ = {self.scale:.5g} " + tr('j.u.·ppm'))
        if self.result is not None:
            r, p = self.result, self.solved_parameters
            maxima = ', '.join(f'{r.diffusion[i]:.4e} m²/s ({r.weights[i]/r.weights.max():.0%})' for i in r.peaks) or '—'
            text += f"\n{tr('Maksima ILT i względne wysokości')}   {maxima}"
            text += (f"\nRMS ILT = {r.rms/self.scale:.4g} Iₛ · "
                     f"λ = {p.regularization:.12g} · D = [{p.d_min:.4e}, {p.d_max:.4e}] m²/s · N = {p.points} · "
                     f"{tr('Rząd efektywny jądra')} = {r.effective_rank} · ‖Lf‖ = {r.penalty_norm:.4g} · {tr('Solver zakończony')}")
            if self.sigma is not None:
                sigma = self.sigma[self.mask]
                mono = np.sqrt(np.average((f['residuals'][self.mask]/self.scale)**2, weights=(sigma.min()/sigma)**2))
                text += f"\n{tr('Ważone RMS mono / ILT')} = {mono:.4g} / {r.weighted_rms:.4g} Iₛ"
            if not p.d_min <= f['D_m2_s'] <= p.d_max:
                text += ' · ' + tr('D mono poza zakresem siatki')
        self.summary.setText(text)

    def export_png(self, plots, name):
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, tr('Zapisz PNG'), f'{name}.png', 'PNG (*.png)')
        if not filename:
            return
        path = Path(filename if filename.lower().endswith('.png') else filename+'.png')
        try:
            if len(plots) == 1:
                save_light_plot(plots[0], path, legend=True, dpi=300)
            else:
                self.save_combined(plots, path)
        except (OSError, ValueError) as exc:
            QtWidgets.QMessageBox.warning(self, tr('Zapisz PNG'), str(exc))

    @staticmethod
    def save_combined(plots, path):
        with TemporaryDirectory() as tmp:
            panels = []
            for index, plot in enumerate(plots):
                panel = Path(tmp)/f'{index}.png'
                save_light_plot(plot, panel, legend=True)
                image = QtGui.QImage(str(panel))
                if image.isNull():
                    raise OSError(tr('Nie udało się zapisać pliku PNG.'))
                panels.append(image)
            margin = 32
            combined = QtGui.QImage(max(image.width() for image in panels)+2*margin,
                                   sum(image.height()+margin for image in panels)+margin,
                                   QtGui.QImage.Format.Format_RGB32)
            combined.fill(QtGui.QColor('white'))
            painter = QtGui.QPainter(combined)
            y = margin
            for image in panels:
                painter.drawImage(margin, y, image)
                y += image.height()+margin
            painter.end()
            combined.setDotsPerMeterX(11811)
            combined.setDotsPerMeterY(11811)
            if not combined.save(str(path), 'PNG'):
                raise OSError(tr('Nie udało się zapisać pliku PNG.'))

    def done(self, result):
        self.closed = True
        self.timer.stop()
        super().done(result)
