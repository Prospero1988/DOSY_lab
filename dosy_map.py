"""ROI-based DOSY projection from the reference spectrum and cached fits."""
from dataclasses import dataclass

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from d_statistics import finite
from spectrum_view import SpectrumBox


@dataclass(frozen=True)
class Band:
    ppm: np.ndarray
    intensity: np.ndarray
    diffusion: float


def fitted_diffusion(result):
    fit = result.get('fit') or {}
    diffusion = finite(fit.get('D_m2_s'))
    return diffusion if fit.get('ok') and diffusion is not None and diffusion > 0 else None


def roi_bands(ppm, reference, results, excluded=()):
    """Use each FFT bin at most once. Mask ambiguous bins in overlapping ROIs."""
    intervals = []
    coverage = np.zeros(len(ppm) + 1, dtype=int)
    for result in results:
        diffusion = fitted_diffusion(result)
        if diffusion is None or result['id'] in excluded:
            continue
        start, stop = np.searchsorted(ppm, result['bounds'])
        if start < stop:
            intervals.append((start, stop, diffusion))
            coverage[start] += 1
            coverage[stop] -= 1
    shared = np.cumsum(coverage[:-1]) > 1
    bands = [Band(ppm[start:stop], np.where(shared[start:stop], 0., reference[start:stop]), diffusion)
             for start, stop, diffusion in intervals if not shared[start:stop].all()]
    return bands, bool(shared.any())


class MapViewBox(pg.ViewBox):
    point_clicked = QtCore.Signal(float, float)
    floor_stepped = QtCore.Signal(float)

    def __init__(self):
        super().__init__(enableMenu=False)
        self.pick_point = False

    def wheelEvent(self, event, axis=None):
        self.floor_stepped.emit(event.delta() / 120)
        event.accept()

    def mouseDragEvent(self, event, axis=None):
        if self.pick_point and event.button() == QtCore.Qt.MouseButton.LeftButton:
            event.accept()
        else:
            super().mouseDragEvent(event, axis)

    def mouseClickEvent(self, event):
        if self.pick_point and event.button() == QtCore.Qt.MouseButton.LeftButton and not event.double():
            point = self.mapToView(event.pos())
            self.point_clicked.emit(point.x(), point.y())
            event.accept()
        else:
            super().mouseClickEvent(event)


class DosyMapWindow(QtWidgets.QMainWindow):
    D_SCALE = 1e-10

    def __init__(self, parent):
        super().__init__(parent, QtCore.Qt.WindowType.Window)
        self.setWindowTitle('DOSY Lab · Mapa DOSY 2D')
        self.resize(1320, 880)
        self.setStyleSheet(parent.styleSheet())
        self.dataset = None
        self.experiment_id = None
        self.results = []
        self.excluded = set()
        self.images = []
        self.stale = True
        self.peak_intensity = 0.
        self.full_range = None
        self.settings = parent.settings
        self.map_color = self.saved_color('dosy_map/color', '#fb923c')
        self.trace_color = self.saved_color('dosy_map/trace_color', '#67e8f9')

        root = QtWidgets.QWidget()
        self.setCentralWidget(root)
        layout = QtWidgets.QVBoxLayout(root)
        toolbar = QtWidgets.QHBoxLayout()
        layout.addLayout(toolbar)
        self.tools = QtWidgets.QButtonGroup(self)
        for index, (name, tip) in enumerate((
            ('Zoom', 'Przeciągnij prostokąt na mapie lub zakres ppm nad widmem 1D.'),
            ('Przesuwaj', 'Przeciągnij, aby przesunąć widok.'),
            ('Punkt / prowadnice', 'Kliknij mapę, aby odczytać ppm i D w miejscu przecięcia prowadnic.'),
        )):
            button = QtWidgets.QPushButton(name)
            button.setToolTip(tip)
            button.setCheckable(True)
            self.tools.addButton(button, index)
            toolbar.addWidget(button)
        self.tools.button(0).setChecked(True)
        self.tools.idClicked.connect(self.set_tool)
        reset = QtWidgets.QPushButton('Cały zakres')
        reset.clicked.connect(self.reset_view)
        toolbar.addWidget(reset)
        toolbar.addStretch()
        toolbar.addWidget(QtWidgets.QLabel('Wzmocnienie 2D'))
        self.gain = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.gain.setRange(0, 40)
        self.gain.setFixedWidth(100)
        self.gain.setToolTip('Wzmocnij jasność mapy od 1× do 100×. Próg podłogi pozostaje niezależny.')
        self.gain.valueChanged.connect(self.update_intensity)
        toolbar.addWidget(self.gain)
        self.gain_label = QtWidgets.QLabel('1×')
        self.gain_label.setMinimumWidth(39)
        toolbar.addWidget(self.gain_label)
        self.floor = self.spin_box(toolbar, 'Podłoga', 0, 99, 0,
                                  'Próg w % najsilniejszej dodatniej intensywności przed wzmocnieniem. Rolka zmienia próg.')
        self.floor.valueChanged.connect(self.update_intensity)
        self.width = self.spin_box(toolbar, 'Szerokość pasma D', .2, 15, 3,
                                  'Szerokość FWHM w % maksymalnego D widocznych ROI. Ustawienie graficzne.')
        self.width.valueChanged.connect(self.render_map)
        view_button = QtWidgets.QToolButton()
        view_button.setText('Widok')
        view_button.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        appearance_menu = QtWidgets.QMenu(view_button)
        self.map_color_action = appearance_menu.addAction('Kolor mapy 2D…')
        self.trace_color_action = appearance_menu.addAction('Kolor widma 1D…')
        view_button.setMenu(appearance_menu)
        toolbar.addWidget(view_button)
        self.grid_button = QtWidgets.QPushButton('Siatka')
        self.grid_button.setCheckable(True)
        self.grid_button.setChecked(self.settings.value('dosy_map/grid', True, type=bool))
        self.grid_button.setToolTip('Pokaż lub ukryj siatkę na mapie i widmie 1D.')
        toolbar.addWidget(self.grid_button)
        self.save_button = QtWidgets.QPushButton('Zapisz widok PNG')
        self.save_button.setEnabled(False)
        self.save_button.clicked.connect(self.choose_save)
        toolbar.addWidget(self.save_button)

        content = QtWidgets.QSplitter()
        layout.addWidget(content, 1)
        self.graphics = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        content.addWidget(self.graphics)
        self.reference_view = SpectrumBox()
        self.projection = pg.PlotWidget(viewBox=self.reference_view)
        self.reference_view.bind_selection_feedback(self.projection)
        self.graphics.addWidget(self.projection)
        self.map_view = MapViewBox()
        self.plot = pg.PlotWidget(viewBox=self.map_view)
        self.graphics.addWidget(self.plot)
        self.graphics.setSizes([210, 540])
        self.graphics.setStretchFactor(1, 1)
        self.projection.setLabel('left', '1D · Re')
        self.projection.getAxis('bottom').setStyle(showValues=False)
        self.plot.setLabel('left', 'D / 10⁻¹⁰ [m²/s]')
        self.plot.setLabel('bottom', '¹H — przesunięcie chemiczne', units='ppm')
        self.plot.getAxis('bottom').enableAutoSIPrefix(False)
        for plot in (self.projection, self.plot):
            plot.invertX(True)
            plot.getAxis('left').setWidth(110)
            plot.viewport().installEventFilter(self)
        self.map_view.setLimits(yMin=0)
        self.projection.setXLink(self.plot)
        self.reference_view.sigRangeChangedManually.connect(self.sync_from_projection)
        self.map_view.floor_stepped.connect(lambda steps: self.floor.setValue(self.floor.value() + steps))
        self.map_view.point_clicked.connect(self.show_coordinates)
        pen = pg.mkPen('#facc15', style=QtCore.Qt.PenStyle.DashLine)
        self.projection_guide = pg.InfiniteLine(angle=90, pen=pen)
        self.ppm_guide = pg.InfiniteLine(angle=90, pen=pen)
        self.diffusion_guide = pg.InfiniteLine(angle=0, pen=pen)
        for plot, guide in ((self.projection, self.projection_guide),
                            (self.plot, self.ppm_guide), (self.plot, self.diffusion_guide)):
            plot.addItem(guide, ignoreBounds=True)
            guide.setZValue(10)
            guide.hide()
        self.trace = self.projection.plot(pen=pg.mkPen(self.trace_color))
        self.map_color_action.triggered.connect(lambda: self.choose_color(map_color=True))
        self.trace_color_action.triggered.connect(lambda: self.choose_color(map_color=False))
        self.grid_button.toggled.connect(self.set_grid)
        self.set_grid(self.grid_button.isChecked())
        self.update_map_color()
        self.rois = QtWidgets.QListWidget()
        self.rois.setMinimumWidth(195)
        self.rois.setToolTip('Odznacz ROI, aby wyłączyć je tylko z mapy.')
        self.rois.itemChanged.connect(self.roi_toggled)
        content.addWidget(self.rois)
        content.setSizes([1040, 250])
        content.setStretchFactor(0, 1)
        self.coordinate = QtWidgets.QLabel('Włącz Punkt / prowadnice, aby odczytać współrzędne.')
        layout.addWidget(self.coordinate)
        self.note = QtWidgets.QLabel()
        self.note.setWordWrap(True)
        layout.addWidget(self.note)
        self.set_tool(0)

    @staticmethod
    def spin_box(toolbar, name, low, high, initial, tip):
        toolbar.addWidget(QtWidgets.QLabel(name))
        control = QtWidgets.QDoubleSpinBox()
        control.setRange(low, high)
        control.setDecimals(1)
        control.setSuffix(' %')
        control.setValue(initial)
        control.setToolTip(tip)
        toolbar.addWidget(control)
        return control

    def set_tool(self, index):
        self.map_view.pick_point = index == 2
        self.map_view.setMouseMode(pg.ViewBox.PanMode if index == 1 else pg.ViewBox.RectMode)
        self.reference_view.mode = 'pan' if index == 1 else 'zoom'

    def eventFilter(self, watched, event):
        if (event.type() == QtCore.QEvent.Type.MouseButtonDblClick
                and event.button() == QtCore.Qt.MouseButton.LeftButton):
            self.reset_view()
            return True
        return super().eventFilter(watched, event)

    def sync_from_projection(self, axes):
        if axes[0]:
            self.map_view.setXRange(*self.reference_view.viewRange()[0], padding=0)

    def show_coordinates(self, ppm, diffusion):
        for guide, value in ((self.projection_guide, ppm), (self.ppm_guide, ppm),
                             (self.diffusion_guide, diffusion)):
            guide.setPos(value)
            guide.show()
        self.coordinate.setText(f'δ = {ppm:.2f} ppm    ·    D = {diffusion*self.D_SCALE:.3e} m²/s')

    def update_data(self, experiment_id, dataset, results):
        new_experiment = self.experiment_id != experiment_id
        self.experiment_id, self.dataset, self.results = experiment_id, dataset, results
        self.stale = False
        self.statusBar().clearMessage()
        if new_experiment:
            self.excluded.clear()
            for guide in (self.projection_guide, self.ppm_guide, self.diffusion_guide):
                guide.hide()
            self.coordinate.setText('Włącz Punkt / prowadnice, aby odczytać współrzędne.')
        self.reference = dataset.spectrum[dataset.options['reference_layer']].real
        self.trace.setData(dataset.ppm, self.reference)
        with QtCore.QSignalBlocker(self.rois):
            self.rois.clear()
            for result in results:
                sid = result['id']
                group = result.get('group')
                diffusion = fitted_diffusion(result)
                label = f'{sid} · {group}' if group else sid
                detail = f'{diffusion:.3e} m²/s' if diffusion is not None else 'brak dodatniego D'
                item = QtWidgets.QListWidgetItem(f'{label}\n{detail}')
                item.setData(QtCore.Qt.ItemDataRole.UserRole, sid)
                item.setToolTip('\n'.join((result.get('fit') or {}).get('warnings', [])))
                if diffusion is None:
                    item.setFlags(item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEnabled)
                else:
                    item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
                    item.setCheckState(QtCore.Qt.CheckState.Unchecked if sid in self.excluded else QtCore.Qt.CheckState.Checked)
                self.rois.addItem(item)
        first_map = not self.images
        self.render_map()
        if new_experiment or first_map:
            self.reset_view()

    def invalidate(self):
        self.stale = True
        self.save_button.setEnabled(False)
        self.statusBar().showMessage('Mapa pokazuje ostatnie ukończone przeliczenie.')

    def roi_toggled(self, item):
        sid = item.data(QtCore.Qt.ItemDataRole.UserRole)
        if item.checkState() == QtCore.Qt.CheckState.Checked:
            self.excluded.discard(sid)
        else:
            self.excluded.add(sid)
        self.render_map()

    def render_map(self):
        if self.dataset is None:
            return
        for image in self.images:
            self.plot.removeItem(image)
        self.images.clear()
        ppm = self.dataset.ppm
        layer = self.dataset.options['reference_layer'] + 1
        bands, overlaps = roi_bands(ppm, self.reference, self.results, self.excluded)
        bands = [band for band in bands if np.any(band.intensity > 0)]
        self.peak_intensity = max((np.max(band.intensity) for band in bands), default=0.)
        self.save_button.setEnabled(bool(bands) and not self.stale)
        if not bands:
            self.full_range = ((ppm[0], ppm[-1]), (0., 1.))
            self.note.setText('Brak dodatnich sygnałów do wyświetlenia. Sprawdź dopasowania D, fazę i wybór ROI.')
            return
        maximum = max(band.diffusion for band in bands) / self.D_SCALE
        fwhm = self.width.value() * maximum / 100
        offset = np.linspace(-2*fwhm, 2*fwhm, 81)
        kernel = np.exp(-4 * np.log(2) * (offset / fwhm)**2)
        d_step = offset[1] - offset[0]
        ppm_step = ppm[1] - ppm[0]
        for band in bands:
            center = band.diffusion / self.D_SCALE
            image = pg.ImageItem(axisOrder='row-major')
            image.setImage((kernel[:, None] * band.intensity).astype(np.float32), autoLevels=False)
            image.setLookupTable(self.lookup_table)
            image.setRect(QtCore.QRectF(band.ppm[0]-ppm_step/2, center+offset[0]-d_step/2,
                                       len(band.ppm)*ppm_step, len(offset)*d_step))
            self.plot.addItem(image)
            self.images.append(image)
        span = maximum - min(band.diffusion for band in bands) / self.D_SCALE
        d_margin = max(2.2*fwhm, .04*span)
        left = min(band.ppm[0] for band in bands)
        right = max(band.ppm[-1] for band in bands)
        x_margin = max(.02, .04*(right-left))
        self.full_range = (
            (max(ppm[0], left-x_margin), min(ppm[-1], right+x_margin)),
            (max(0., maximum - span - d_margin), maximum + d_margin),
        )
        self.note.setText(
            f'Widmo 1D: warstwa {layer} (referencja). Mapa: kształt widma w ROI przy jego D. '
            'Kolor pokazuje dodatnie intensywności. Szerokość w D jest ustawieniem graficznym.'
            + (' Wspólne punkty nakładających się ROI są pominięte.' if overlaps else ''))
        self.update_intensity()

    def update_intensity(self):
        gain = 10 ** (self.gain.value() / 20)
        self.gain_label.setText(f'{gain:.3g}×')
        threshold = self.peak_intensity * self.floor.value() / 100
        saturation = threshold + (self.peak_intensity - threshold) / gain
        for image in self.images:
            image.setLevels((threshold, saturation))

    def saved_color(self, key, default):
        color = QtGui.QColor(self.settings.value(key, default, type=str))
        return color if color.isValid() else QtGui.QColor(default)

    def choose_color(self, *, map_color):
        current = self.map_color if map_color else self.trace_color
        target = 'mapy 2D' if map_color else 'widma 1D'
        color = QtWidgets.QColorDialog.getColor(current, self, f'Kolor {target}')
        if not color.isValid():
            return
        if map_color:
            self.map_color = color
            self.settings.setValue('dosy_map/color', color.name())
            self.update_map_color()
        else:
            self.trace_color = color
            self.settings.setValue('dosy_map/trace_color', color.name())
            self.trace.setPen(pg.mkPen(color))

    @staticmethod
    def color_table(color):
        # ImageItem clips values below the lower level to the transparent first entry.
        colors = np.empty((256, 4), dtype=np.uint8)
        colors[:, :3] = color.getRgb()[:3]
        colors[:, 3] = np.round(255 * np.sqrt(np.linspace(0., 1., len(colors)))).astype(np.uint8)
        return colors

    def update_map_color(self):
        self.lookup_table = self.color_table(self.map_color)
        for image in self.images:
            image.setLookupTable(self.lookup_table)

    def set_grid(self, visible):
        for plot in (self.projection, self.plot):
            plot.showGrid(x=visible, y=visible, alpha=.17)
        self.settings.setValue('dosy_map/grid', visible)

    def reset_view(self):
        if self.full_range is None:
            return
        x, y = self.full_range
        self.plot.setRange(xRange=x, yRange=y, padding=0)
        self.projection.setXRange(*x, padding=0)
        self.reference_view.enableAutoRange(axis=pg.ViewBox.YAxis)

    def choose_save(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Zapisz aktywny widok DOSY 2D',
                                                       'DOSY_mapa.png', 'PNG (*.png)')
        if path:
            if not path.lower().endswith('.png'):
                path += '.png'
            if not self.save_view(path):
                QtWidgets.QMessageBox.warning(self, 'Zapis widoku', 'Nie udało się zapisać pliku PNG.')

    def save_view(self, path):
        # Render a separate copy so export never changes the live plots or settings.
        canvas = QtWidgets.QWidget()
        canvas.setAttribute(QtCore.Qt.WidgetAttribute.WA_DontShowOnScreen)
        canvas.setStyleSheet('background: white;')
        layout = QtWidgets.QVBoxLayout(canvas)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(3)
        try:
            spectrum = self.export_plot(self.projection, layout)
            map_plot = self.export_plot(self.plot, layout)
            navy = QtGui.QColor('#142b5c')
            spectrum.plot(*self.trace.getData(), pen=pg.mkPen(navy, width=1.3))
            lookup = self.color_table(navy)
            for image in self.images:
                copy = pg.ImageItem(axisOrder='row-major')
                copy.setImage(image.image, levels=image.getLevels(), lut=lookup)
                copy.setTransform(image.transform())
                map_plot.addItem(copy)
            pen = pg.mkPen('k', style=QtCore.Qt.PenStyle.DashLine)
            for plot, guide in ((spectrum, self.projection_guide), (map_plot, self.ppm_guide),
                                (map_plot, self.diffusion_guide)):
                if guide.isVisible():
                    line = pg.InfiniteLine(pos=guide.value(), angle=guide.angle, pen=pen)
                    line.setZValue(10)
                    plot.addItem(line, ignoreBounds=True)
            canvas.resize(self.graphics.width(), layout.sizeHint().height())
            canvas.show()
            for source, plot in ((self.projection, spectrum), (self.plot, map_plot)):
                x_range, y_range = source.getViewBox().viewRange()
                plot.setRange(xRange=x_range, yRange=y_range, padding=0)
            return canvas.grab().save(str(path), 'PNG')
        finally:
            canvas.close()
            canvas.deleteLater()

    def export_plot(self, source, layout):
        plot = pg.PlotWidget(background='w')
        layout.addWidget(plot)
        plot.setFixedHeight(source.height())
        plot.invertX(True)
        plot.hideButtons()
        plot.getViewBox().setBorder('k')
        for name in ('left', 'bottom'):
            original = source.getAxis(name)
            axis = plot.getAxis(name)
            axis.setPen('k')
            axis.setTickPen('k')
            axis.setTextPen('k')
            axis.setLabel(original.labelText, units=original.labelUnits, color='#000000')
            axis.setStyle(showValues=original.style['showValues'])
            axis.enableAutoSIPrefix(original.autoSIPrefix)
        plot.getAxis('left').setWidth(source.getAxis('left').width())
        grid = self.grid_button.isChecked()
        plot.showGrid(x=grid, y=grid, alpha=.17)
        return plot
