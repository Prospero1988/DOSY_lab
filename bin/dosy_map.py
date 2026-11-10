"""ROI-based DOSY projection from the reference spectrum and cached fits."""
from i18n import tr

from dataclasses import dataclass

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from d_statistics import finite, quality
from dosy_alignment_panel import AlignmentPanel, ILT_MODES, MODE_NAMES, ORIGINAL
from plot_export import NAVY, light_plot
from spectrum_view import SpectrumBox
from storage import source_save_name


EXPLORATORY, QC_FILTERED, MONO_REFERENCE = range(3)
POLICY_NAMES = ('Eksploracyjna', 'Filtrowana QC', 'Odniesienie monoeksponencjalne')


@dataclass(frozen=True)
class Band:
    ppm: np.ndarray
    intensity: np.ndarray
    diffusion: float
    sid: str = ''


def diffusion_projection(bands, offset, kernel, diffusion_scale):
    """Sum the positive map intensity along ppm on a common display D grid."""
    centers = np.array([band.diffusion/diffusion_scale for band in bands])
    step = offset[1]-offset[0]
    start, stop = max(0., centers.min()+offset[0]), centers.max()+offset[-1]
    grid = np.linspace(start, stop, int(np.ceil((stop-start)/step))+1)
    total = np.zeros_like(grid)
    for center, band in zip(centers, bands):
        row_sums = np.maximum(band.intensity, 0.).sum()*kernel
        total += np.interp(grid, center+offset, row_sums, left=0., right=0.)
    return total, grid


def fitted_diffusion(result):
    fit = result.get('fit') or {}
    diffusion = finite(fit.get('D_m2_s'))
    return diffusion if fit.get('ok') and diffusion is not None and diffusion > 0 else None


def roi_bands(ppm, reference, results, excluded=(), *, overrides=None, profiles=None, reconstructed_only=False):
    """Use each FFT bin at most once. Mask ambiguous bins in overlapping ROIs."""
    overrides, profiles = overrides or {}, profiles or {}
    intervals = []
    coverage = np.zeros(len(ppm) + 1, dtype=int)
    for result in results:
        diffusion = fitted_diffusion(result)
        if result['id'] in excluded or (not reconstructed_only and (
                diffusion is None or profiles.get(result['id']) == ())):
            continue
        start, stop = np.searchsorted(ppm, result['bounds'])
        if start < stop:
            intervals.append((start, stop, result['id'], diffusion))
            coverage[start] += 1
            coverage[stop] -= 1
    shared = np.cumsum(coverage[:-1]) > 1
    bands = []
    for start, stop, sid, diffusion in intervals:
        mask = shared[start:stop]
        if mask.all() or (reconstructed_only and sid not in profiles):
            continue
        for center, intensity in profiles.get(sid, ((overrides.get(sid, diffusion), reference[start:stop]),)):
            bands.append(Band(ppm[start:stop], np.where(mask, 0., intensity), center, sid))
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


class DiffusionProjectionBox(pg.ViewBox):
    def wheelEvent(self, event, axis=None):
        factor = float(np.exp(np.clip(-event.delta()/120*np.log(1.2), -4, 4)))
        self.scaleBy(x=factor, center=QtCore.QPointF(0, 0))
        event.accept()


class DosyMapWindow(QtWidgets.QMainWindow):
    D_SCALE = 1e-10

    def __init__(self, parent):
        super().__init__(parent, QtCore.Qt.WindowType.Window)
        self.setWindowTitle(tr('DOSY Lab · Mapa DOSY 2D'))
        self.resize(1320, 880)
        self.setStyleSheet(parent.styleSheet())
        self.dataset = None
        self.experiment_id = None
        self.results = []
        self.excluded = set()
        self.excluded_groups = set()
        self.images = []
        self.unresolved_curves = []
        self.native_guides = []
        self.overlay_items = []
        self.map_legend = ''
        self.represented_peaks = {}
        self.stale = True
        self.peak_intensity = 0.
        self.full_range = None
        self.rendered_mode = (ORIGINAL, EXPLORATORY)
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
            ('Zoom [Z]', tr('Przeciągnij prostokąt na mapie lub zakres ppm nad widmem 1D.')),
            (tr('Przesuwaj [P]'), tr('Przeciągnij, aby przesunąć widok.')),
            (tr('Punkt / prowadnice'), tr('Kliknij mapę, aby odczytać ppm i D. Delete usuwa prowadnice.')),
        )):
            button = QtWidgets.QPushButton(name)
            button.setToolTip(tip)
            button.setCheckable(True)
            self.tools.addButton(button, index)
            toolbar.addWidget(button)
        self.tools.button(0).setChecked(True)
        self.tools.idClicked.connect(self.set_tool)
        reset = QtWidgets.QPushButton(tr('Cały zakres [F]'))
        reset.clicked.connect(self.reset_view)
        toolbar.addWidget(reset)
        for key, callback in [('Z', self.tools.button(0).click), ('F', self.reset_view), ('P', self.tools.button(1).click)]:
            shortcut = QtGui.QShortcut(QtGui.QKeySequence(key), self)
            shortcut.setAutoRepeat(False)
            shortcut.activated.connect(callback)
        fullscreen = QtWidgets.QPushButton(tr('Pełny ekran'))
        fullscreen.clicked.connect(self.toggle_fullscreen)
        toolbar.addWidget(fullscreen)
        escape = QtGui.QShortcut(QtGui.QKeySequence('Esc'), self)
        escape.activated.connect(lambda: self.showNormal() if self.isFullScreen() else None)
        toolbar.addStretch()
        toolbar.addWidget(QtWidgets.QLabel(tr('Wzmocnienie 2D')))
        self.gain = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.gain.setRange(0, 40)
        self.gain.setFixedWidth(100)
        self.gain.setToolTip(tr('Wzmocnij jasność mapy od 1× do 100×. Próg podłogi pozostaje niezależny.'))
        self.gain.valueChanged.connect(self.update_intensity)
        toolbar.addWidget(self.gain)
        self.gain_label = QtWidgets.QLabel('1×')
        self.gain_label.setMinimumWidth(39)
        toolbar.addWidget(self.gain_label)
        self.floor = self.spin_box(toolbar, tr('Podłoga'), 0, 99, 0,
                                  tr('Próg w % najsilniejszej dodatniej intensywności przed wzmocnieniem. Rolka zmienia próg.'))
        self.floor.valueChanged.connect(self.update_intensity)
        self.width = self.spin_box(toolbar, tr('Szerokość pasma D'), .2, 15, 3,
                                  tr('Szerokość FWHM w % maksymalnego D widocznych ROI. Ustawienie graficzne.'))
        self.width.valueChanged.connect(self.render_map)
        view_button = QtWidgets.QToolButton()
        view_button.setText(tr('Widok'))
        view_button.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        appearance_menu = QtWidgets.QMenu(view_button)
        self.map_color_action = appearance_menu.addAction(tr('Kolor mapy 2D…'))
        self.trace_color_action = appearance_menu.addAction(tr('Kolor widma 1D…'))
        view_button.setMenu(appearance_menu)
        toolbar.addWidget(view_button)
        self.grid_button = QtWidgets.QPushButton(tr('Siatka'))
        self.grid_button.setCheckable(True)
        self.grid_button.setChecked(self.settings.value('dosy_map/grid', True, type=bool))
        self.grid_button.setToolTip(tr('Pokaż lub ukryj siatkę na mapie i widmie 1D.'))
        toolbar.addWidget(self.grid_button)
        self.diffusion_button = QtWidgets.QPushButton(tr('Projekcja D'))
        self.diffusion_button.setCheckable(True)
        self.diffusion_button.setToolTip(tr('Pokaż po lewej sumę dodatnich intensywności wzdłuż osi ppm dla każdego D.'))
        toolbar.addWidget(self.diffusion_button)
        delete = QtGui.QShortcut(QtGui.QKeySequence(QtCore.Qt.Key.Key_Delete), self)
        delete.activated.connect(self.clear_guides)
        self.save_button = QtWidgets.QPushButton(tr('Zapisz widok PNG'))
        self.save_button.setEnabled(False)
        self.save_button.clicked.connect(self.choose_save)
        toolbar.addWidget(self.save_button)

        alignment_toolbar = QtWidgets.QHBoxLayout()
        layout.addLayout(alignment_toolbar)
        adaptive_toolbar = QtWidgets.QHBoxLayout()
        layout.addLayout(adaptive_toolbar)
        content = QtWidgets.QSplitter()
        layout.addWidget(content, 1)
        self.graphics = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        content.addWidget(self.graphics)
        self.reference_view = SpectrumBox()
        self.projection = pg.PlotWidget(viewBox=self.reference_view)
        self.reference_view.bind_selection_feedback(self.projection)
        self.map_view = MapViewBox()
        self.plot = pg.PlotWidget(viewBox=self.map_view)
        self.diffusion_plot = pg.PlotWidget(viewBox=DiffusionProjectionBox(enableMenu=False))
        self.diffusion_plot.setFixedWidth(185)
        self.diffusion_spacer = QtWidgets.QWidget()
        self.diffusion_spacer.setFixedWidth(self.diffusion_plot.width())
        self.graphics.addWidget(self.plot_row(self.projection, self.diffusion_spacer))
        self.graphics.addWidget(self.plot_row(self.plot, self.diffusion_plot))
        self.graphics.setSizes([210, 540])
        self.graphics.setStretchFactor(1, 1)
        self.projection.setLabel('left', '1D · Re')
        self.projection.getAxis('bottom').setStyle(showValues=False)
        self.plot.setLabel('left', 'D / 10⁻¹⁰ [m²/s]')
        self.plot.setLabel('bottom', tr('¹H — przesunięcie chemiczne'), units='ppm')
        self.plot.getAxis('bottom').enableAutoSIPrefix(False)
        self.diffusion_plot.setLabel('bottom', 'Σ I')
        self.diffusion_plot.getAxis('bottom').enableAutoSIPrefix(False)
        self.diffusion_plot.hideAxis('left')
        for plot in (self.diffusion_plot, self.plot):
            plot.getAxis('bottom').setHeight(45)
        self.diffusion_plot.invertX(True)
        self.diffusion_plot.setMouseEnabled(x=False, y=False)
        self.diffusion_plot.hideButtons()
        self.diffusion_plot.setLimits(xMin=0)
        for plot in (self.projection, self.plot):
            plot.invertX(True)
            plot.getAxis('left').setWidth(110)
            plot.viewport().installEventFilter(self)
        self.map_view.setLimits(yMin=0)
        self.projection.setXLink(self.plot)
        self.diffusion_plot.setYLink(self.plot)
        self.reference_view.sigRangeChangedManually.connect(self.sync_from_projection)
        self.map_view.floor_stepped.connect(lambda steps: self.floor.setValue(self.floor.value() + steps))
        self.map_view.point_clicked.connect(self.show_coordinates)
        pen = pg.mkPen('#facc15', style=QtCore.Qt.PenStyle.DashLine)
        self.projection_guide = pg.InfiniteLine(angle=90, pen=pen)
        self.ppm_guide = pg.InfiniteLine(angle=90, pen=pen)
        self.diffusion_guide = pg.InfiniteLine(angle=0, pen=pen)
        self.diffusion_projection_guide = pg.InfiniteLine(angle=0, pen=pen)
        for plot, guide in ((self.projection, self.projection_guide),
                            (self.plot, self.ppm_guide), (self.plot, self.diffusion_guide),
                            (self.diffusion_plot, self.diffusion_projection_guide)):
            plot.addItem(guide, ignoreBounds=True)
            guide.setZValue(10)
            guide.hide()
        self.d_label = pg.TextItem(color='#facc15', fill=pg.mkBrush('#101827'))
        self.d_label.setZValue(11)
        self.plot.addItem(self.d_label, ignoreBounds=True)
        self.d_label.hide()
        self.map_view.sigRangeChanged.connect(self.position_d_label)
        self.trace = self.projection.plot(pen=pg.mkPen(self.trace_color))
        self.diffusion_trace = self.diffusion_plot.plot(pen=pg.mkPen(self.trace_color))
        self.diffusion_button.setChecked(self.settings.value('dosy_map/diffusion_projection', False, type=bool))
        self.set_diffusion_projection_visible(self.diffusion_button.isChecked())
        self.diffusion_button.toggled.connect(self.set_diffusion_projection_visible)
        self.map_color_action.triggered.connect(lambda: self.choose_color(map_color=True))
        self.trace_color_action.triggered.connect(lambda: self.choose_color(map_color=False))
        self.grid_button.toggled.connect(self.set_grid)
        self.set_grid(self.grid_button.isChecked())
        self.update_map_color()
        sidebar = QtWidgets.QWidget()
        sidebar.setMinimumWidth(195)
        sidebar_layout = QtWidgets.QVBoxLayout(sidebar)
        sidebar_layout.setContentsMargins(0, 0, 0, 0)
        self.groups_label = QtWidgets.QLabel(tr('Grupy / związki'))
        sidebar_layout.addWidget(self.groups_label)
        self.groups = QtWidgets.QTreeWidget()
        self.groups.setHeaderHidden(True)
        self.groups.setRootIsDecorated(False)
        self.groups.setMaximumHeight(140)
        self.groups.setToolTip(tr('Ukryj całą grupę na mapie. Ponowne włączenie zachowuje wybór pojedynczych sygnałów.'))
        self.groups.itemChanged.connect(self.group_toggled)
        sidebar_layout.addWidget(self.groups)
        sidebar_layout.addWidget(QtWidgets.QLabel(tr('Sygnały')))
        self.rois = QtWidgets.QListWidget()
        self.rois.setToolTip(tr('Odznacz ROI, aby wyłączyć je tylko z mapy.'))
        self.rois.itemChanged.connect(self.roi_toggled)
        self.rois.itemClicked.connect(lambda item: self.statusBar().showMessage(
            item.toolTip().replace('\n',' · ')))
        sidebar_layout.addWidget(self.rois, 1)
        self.policy = QtWidgets.QComboBox()
        self.policy.addItems([tr(label) for label in POLICY_NAMES])
        self.policy.setToolTip(tr('Polityka wizualizacji ILT'))
        self.policy.currentIndexChanged.connect(self.refresh_alignment)
        self.policy.hide()
        alignment_toolbar.addWidget(self.policy)
        self.overlay = QtWidgets.QComboBox()
        self.overlay.addItems([tr('Tylko profile na mapie'),
                               tr('ILT + nierozdzielone ROI na 1D'),
                               tr('ILT + nierozdzielone ROI i odniesienia mono D')])
        self.overlay.setCurrentIndex(1)
        self.overlay.currentIndexChanged.connect(self.render_map)
        self.overlay.hide()
        sidebar_layout.addWidget(self.overlay)
        self.alignment = AlignmentPanel(self, alignment_toolbar, adaptive_toolbar, sidebar_layout)
        content.addWidget(sidebar)
        content.setSizes([1040, 250])
        content.setStretchFactor(0, 1)
        self.coordinate = QtWidgets.QLabel()
        layout.addWidget(self.coordinate)
        self.clear_guides()
        self.note = QtWidgets.QLabel()
        self.note.setWordWrap(True)
        layout.addWidget(self.note)
        self.alignment.changed.connect(self.refresh_alignment)
        self.set_tool(0)

    @staticmethod
    def plot_row(plot, left):
        row = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(left)
        layout.addWidget(plot, 1)
        return row

    def set_diffusion_projection_visible(self, visible):
        self.diffusion_spacer.setVisible(visible)
        self.diffusion_plot.setVisible(visible)
        self.settings.setValue('dosy_map/diffusion_projection', visible)

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

    def toggle_fullscreen(self):
        if self.isFullScreen():
            self.showNormal()
        else:
            self.showFullScreen()

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
                             (self.diffusion_guide, diffusion), (self.diffusion_projection_guide, diffusion)):
            guide.setPos(value)
            guide.show()
        value = f'D = {diffusion*self.D_SCALE:.3e} m²/s'
        self.coordinate.setText(f'δ = {ppm:.2f} ppm    ·    {value}')
        self.d_label.setText(value)
        self.position_d_label()
        self.d_label.show()
        self.d_label.updateTransform(force=True)

    def clear_guides(self):
        for item in (self.projection_guide, self.ppm_guide, self.diffusion_guide,
                     self.diffusion_projection_guide, self.d_label):
            item.hide()
        self.coordinate.setText(tr('Włącz Punkt / prowadnice, aby odczytać współrzędne.'))

    def position_d_label(self):
        ppm, diffusion = self.ppm_guide.value(), self.diffusion_guide.value()
        (left, right), (bottom, top) = self.map_view.viewRange()
        # The ppm axis runs from high to low, so small ppm is near the right edge.
        x_anchor = 1.04 if ppm < left + .25*(right-left) else -.04
        y_anchor = -.12 if diffusion > top - .15*(top-bottom) else 1.12
        self.d_label.setAnchor((x_anchor, y_anchor))
        self.d_label.setPos(ppm, diffusion)

    def update_data(self, experiment_id, dataset, results):
        new_experiment = self.experiment_id != experiment_id
        self.experiment_id, self.dataset, self.results = experiment_id, dataset, results
        self.stale = False
        self.statusBar().clearMessage()
        if new_experiment:
            self.excluded.clear()
            self.excluded_groups.clear()
            self.clear_guides()
        self.reference = dataset.spectrum[dataset.options['reference_layer']].real
        self.trace.setData(dataset.ppm, self.reference)
        groups = sorted({r['group'] for r in results if r.get('group') and fitted_diffusion(r) is not None}, key=str.casefold)
        self.groups_label.setVisible(bool(groups))
        self.groups.setVisible(bool(groups))
        with QtCore.QSignalBlocker(self.groups):
            self.groups.clear()
            for group in groups:
                item = QtWidgets.QTreeWidgetItem([group])
                item.setCheckState(0, QtCore.Qt.CheckState.Unchecked if group in self.excluded_groups else QtCore.Qt.CheckState.Checked)
                self.groups.addTopLevelItem(item)
        with QtCore.QSignalBlocker(self.rois):
            self.rois.clear()
            for result in results:
                sid = result['id']
                group = result.get('group')
                diffusion = fitted_diffusion(result)
                label = f'{sid} · {group}' if group else sid
                detail = f'{diffusion:.3e} m²/s' if diffusion is not None else tr('brak dodatniego D')
                item = QtWidgets.QListWidgetItem(f'{label}\n{detail}')
                item.setData(QtCore.Qt.ItemDataRole.UserRole, sid)
                item.setToolTip('\n'.join((result.get('fit') or {}).get('warnings', [])))
                item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(QtCore.Qt.CheckState.Unchecked if sid in self.excluded else QtCore.Qt.CheckState.Checked)
                self.rois.addItem(item)
        first_map = not self.images
        with QtCore.QSignalBlocker(self.alignment):
            self.alignment.set_data(dataset, results, new_experiment)
        self.render_map()
        if new_experiment or first_map:
            self.reset_view()
            self.rendered_mode = (self.alignment.mode.currentIndex(), self.policy.currentIndex())

    def invalidate(self):
        self.stale = True
        self.save_button.setEnabled(False)
        self.statusBar().showMessage(tr('Mapa pokazuje ostatnie ukończone przeliczenie.'))

    def roi_toggled(self, item):
        sid = item.data(QtCore.Qt.ItemDataRole.UserRole)
        if item.checkState() == QtCore.Qt.CheckState.Checked:
            self.excluded.discard(sid)
        else:
            self.excluded.add(sid)
        self.render_map()

    def group_toggled(self, item, column):
        group = item.text(column)
        if item.checkState(column) == QtCore.Qt.CheckState.Checked:
            self.excluded_groups.discard(group)
        else:
            self.excluded_groups.add(group)
        self.render_map()

    def refresh_alignment(self):
        self.render_map()
        mode = self.alignment.mode.currentIndex(), self.policy.currentIndex()
        if self.alignment.ready and self.full_range is not None and mode != self.rendered_mode:
            self.plot.setYRange(*self.full_range[1], padding=0.)
            self.rendered_mode = mode

    def render_map(self):
        if self.dataset is None:
            return
        for image in self.images:
            self.plot.removeItem(image)
        self.images.clear()
        for plot, item in self.overlay_items:
            plot.removeItem(item)
        self.overlay_items.clear()
        self.unresolved_curves.clear()
        self.native_guides.clear()
        self.map_legend = ''
        self.represented_peaks.clear()
        ppm = self.dataset.ppm
        layer = self.dataset.options['reference_layer'] + 1
        visible = [result for result in self.results if result.get('group') not in self.excluded_groups]
        mode = self.alignment.mode.currentIndex()
        policy = self.policy.currentIndex()
        ilt_mode = mode in ILT_MODES
        mono_reference = ilt_mode and policy == MONO_REFERENCE
        reconstruct = ilt_mode and not mono_reference
        self.policy.setVisible(ilt_mode)
        self.overlay.setVisible(reconstruct)
        overrides, hidden, profiles = self.alignment.map_data(exploratory=reconstruct and policy==EXPLORATORY)
        if mono_reference:
            overrides, profiles = {}, {}
        selected = [r for r in visible if r['id'] not in self.excluded | hidden
                    and (not reconstruct or profiles.get(r['id']) != ())]
        bands, overlaps = roi_bands(ppm, self.reference, selected, overrides=overrides,
                                    profiles=profiles, reconstructed_only=reconstruct)
        bands = [band for band in bands if np.any(band.intensity > 0)]
        for band in bands:
            self.represented_peaks[band.sid] = max(self.represented_peaks.get(band.sid,0.),
                                                  float(np.max(band.intensity)))
        if ilt_mode and (self.alignment.ready or mono_reference):
            self.render_ilt_overlays(selected, set(self.represented_peaks), profiles, overlaps, policy)
        else:
            warnings = {row['id']: '\n'.join((row.get('fit') or {}).get('warnings', []))
                        for row in self.results}
            with QtCore.QSignalBlocker(self.rois):
                for index in range(self.rois.count()):
                    item = self.rois.item(index)
                    item.setToolTip(warnings.get(item.data(QtCore.Qt.ItemDataRole.UserRole), ''))
        self.peak_intensity = max((np.max(band.intensity) for band in bands), default=0.)
        self.save_button.setEnabled(bool(bands or self.unresolved_curves) and not self.stale
                                    and (self.alignment.ready or mono_reference))
        if not bands:
            self.diffusion_trace.setData([], [])
            upper = max((record[2] for record in self.native_guides),default=.95)
            self.full_range = ((ppm[0], ppm[-1]), (0., 1.05*upper))
            note = self.map_legend or self.alignment.error
            if not note and not self.alignment.ready:
                note = tr('Obliczanie rekonstrukcji ILT…')
            self.note.setText(note or tr('Brak dodatnich sygnałów do wyświetlenia. Sprawdź dopasowania D, fazę i wybór ROI.'))
            return
        maximum = max(band.diffusion for band in bands) / self.D_SCALE
        fwhm = self.width.value() * maximum / 100
        offset = np.linspace(-2*fwhm, 2*fwhm, 81)
        kernel = np.exp(-4 * np.log(2) * (offset / fwhm)**2)
        self.diffusion_trace.setData(*diffusion_projection(bands, offset, kernel, self.D_SCALE))
        d_step = offset[1] - offset[0]
        ppm_step = ppm[1] - ppm[0]
        for band in bands:
            center = band.diffusion / self.D_SCALE
            image = pg.ImageItem(axisOrder='row-major')
            image.setImage((kernel[:, None] * band.intensity).astype(np.float32), autoLevels=False)
            image.setLookupTable(self.lookup_table)
            if (self.alignment.highlighted and band.sid not in self.alignment.highlighted
                    and (band.sid, band.diffusion) not in self.alignment.highlighted_components):
                image.setOpacity(.25)
            image.setRect(QtCore.QRectF(band.ppm[0]-ppm_step/2, center+offset[0]-d_step/2,
                                       len(band.ppm)*ppm_step, len(offset)*d_step))
            self.plot.addItem(image)
            self.images.append(image)
        span = maximum - min(band.diffusion for band in bands) / self.D_SCALE
        d_margin = max(2.2*fwhm, .04*span)
        left = min(band.ppm[0] for band in bands)
        right = max(band.ppm[-1] for band in bands)
        x_margin = max(.02, .04*(right-left))
        x_left = min([left]+[curve[0][0] for curve in self.unresolved_curves])
        x_right = max([right]+[curve[0][-1] for curve in self.unresolved_curves])
        native = [record[2] for record in self.native_guides]
        self.full_range = (
            (max(ppm[0], x_left-x_margin), min(ppm[-1], x_right+x_margin)),
            (max(0., min([maximum-span]+native)-d_margin), max([maximum]+native)+d_margin),
        )
        self.note.setText(
            tr(f'Widmo 1D: warstwa {layer} (referencja). Mapa: kształt widma w ROI przy jego D. '
            'Kolor pokazuje dodatnie intensywności. Szerokość w D jest ustawieniem graficznym.')
            + (tr(' Wspólne punkty nakładających się ROI są pominięte.') if overlaps else ''))
        if mode != ORIGINAL and not mono_reference:
            detail = self.map_legend if ilt_mode else tr('Poza zaakceptowanymi grupami lub profilami zachowano oryginalne D.')
            if self.alignment.error:
                detail = tr('Błąd wyrównania') + ' · ' + tr(self.alignment.error)
            elif not self.alignment.ready:
                detail = tr('Obliczanie wyrównania. Tymczasowo wyświetlana jest mapa oryginalna.')
            self.note.setText(tr(MODE_NAMES[self.alignment.mode.currentIndex()]) + ' · ' + detail)
        self.update_intensity()

    def render_ilt_overlays(self, selected, represented, profiles, overlap, policy):
        ids = {row['id'] for row in selected}
        unassigned = ids-represented
        explored = self.alignment.exploratory
        filtered = self.alignment.profiles
        records = self.alignment.exploratory_diagnostics
        analyzed = sum(bool(records.get(sid,{}).get('solved')) for sid in ids)
        attempted = sum(bool(records.get(sid,{}).get('optimized')) for sid in ids)
        quality_exclusions = (set(explored)-set(filtered)) & ids
        unresolved = unassigned-quality_exclusions if policy==QC_FILTERED else unassigned
        by_id = {row['id']:row for row in selected}
        with QtCore.QSignalBlocker(self.rois):
            for index in range(self.rois.count()):
                item = self.rois.item(index)
                sid = item.data(QtCore.Qt.ItemDataRole.UserRole)
                if sid not in ids:
                    detail = tr('ROI wyłączone z widoku.')
                else:
                    detail = self.alignment.profile_detail(sid)
                    fit = by_id[sid]['fit']
                    _, level, explanation = quality(fit)
                    if explanation:
                        detail += '\n' + tr(level) + ' · ' + tr(explanation)
                    if sid in unassigned and sid in profiles:
                        detail += '\n' + tr('Profil nie daje dodatniej mapy po maskowaniu wspólnych binów.')
                item.setToolTip(detail)
        for row in selected:
            sid = row['id']
            if policy == MONO_REFERENCE or sid not in unassigned or not self.overlay.currentIndex():
                continue
            lo, hi = sorted(row['bounds'])
            mask = (self.dataset.ppm>=lo)&(self.dataset.ppm<=hi)
            x = self.dataset.ppm[mask]
            if len(x):
                y = self.reference[mask]
                self.unresolved_curves.append((x,y))
                item = self.projection.plot(x,y,pen=pg.mkPen('#f59e0b',width=2))
                self.overlay_items.append((self.projection,item))
            native = fitted_diffusion(row)
            if self.overlay.currentIndex()==2 and native is not None:
                d = native/self.D_SCALE
                self.native_guides.append((lo,hi,d))
                item = self.plot.plot([lo,hi],[d,d],pen=pg.mkPen('#c084fc',width=2,
                    style=QtCore.Qt.PenStyle.DashLine))
                self.overlay_items.append((self.plot,item))
        if policy == MONO_REFERENCE:
            self.map_legend = tr('Odniesienie monoeksponencjalne') + f" · {tr('ROI na mapie')}: {len(represented)}/{len(self.results)}"
        else:
            self.map_legend = ' · '.join(f'{tr(label)}: {count}' for label,count in (
                ('Wszystkie ROI',len(self.results)), ('Próba ILT',attempted), ('Zbieżne numerycznie',analyzed),
                ('ROI w rekonstrukcji ILT',len(represented)), ('ROI pomijane przez filtr QC',len(quality_exclusions)),
                ('nierozdzielone / bez przypisania',len(unresolved)), ('wyłączone z widoku',len(self.results)-len(ids))))
            self.map_legend = tr(POLICY_NAMES[policy]) + ' · ' + self.map_legend
            if self.overlay.currentIndex():
                self.map_legend += ' · ' + tr('Bursztyn na 1D = zmierzony profil bez przypisania')
            if self.overlay.currentIndex()==2:
                self.map_legend += ' · ' + tr('Fioletowa kreska = pozycja mono D, tylko odniesienie')
        if overlap:
            self.map_legend += ' · ' + tr('Wspólne biny pominięto w macierzy.')
        self.map_legend += ' · ' + tr('Wybierz ROI z listy, aby odczytać przyczynę.')

    def update_intensity(self):
        gain = 10 ** (self.gain.value() / 20)
        self.gain_label.setText(f'{gain:.3g}×')
        threshold = self.peak_intensity * self.floor.value() / 100
        saturation = threshold + (self.peak_intensity - threshold) / gain
        for image in self.images:
            image.setLevels((threshold, saturation))
        if self.map_legend:
            self.note.setText(self.reconstruction_title() + self.ilt_legend())

    def reconstruction_title(self):
        mode = self.alignment.mode.currentIndex()
        return '' if mode in ILT_MODES and self.policy.currentIndex()==MONO_REFERENCE else tr(MODE_NAMES[mode])+' · '

    def ilt_legend(self):
        threshold = self.peak_intensity*self.floor.value()/100
        hidden = sum(value<=threshold for value in self.represented_peaks.values())
        return self.map_legend + (f" · {tr('ROI poniżej progu wyświetlania')}: {hidden}" if hidden else '')

    def saved_color(self, key, default):
        color = QtGui.QColor(self.settings.value(key, default, type=str))
        return color if color.isValid() else QtGui.QColor(default)

    def choose_color(self, *, map_color):
        current = self.map_color if map_color else self.trace_color
        title = tr('Kolor mapy 2D') if map_color else tr('Kolor widma 1D')
        color = QtWidgets.QColorDialog.getColor(current, self, title)
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
            self.diffusion_trace.setPen(pg.mkPen(color))

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
        for plot in (self.projection, self.plot, self.diffusion_plot):
            plot.showGrid(x=visible, y=visible, alpha=.17)
        self.settings.setValue('dosy_map/grid', visible)

    def reset_view(self):
        if self.full_range is None:
            return
        x, y = self.full_range
        self.plot.setRange(xRange=x, yRange=y, padding=0)
        self.projection.setXRange(*x, padding=0)
        self.reference_view.enableAutoRange(axis=pg.ViewBox.YAxis)
        self.diffusion_plot.enableAutoRange(axis=pg.ViewBox.XAxis)

    def choose_save(self):
        filename = source_save_name(self.dataset.name, '.png')
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, tr('Zapisz aktywny widok DOSY 2D'),
                                                       filename, 'PNG (*.png)')
        if path:
            if not path.lower().endswith('.png'):
                path += '.png'
            if not self.save_view(path):
                QtWidgets.QMessageBox.warning(self, tr('Zapis widoku'), tr('Nie udało się zapisać pliku PNG.'))

    def save_view(self, path):
        # Render a separate copy so export never changes the live plots or settings.
        mono_reference = self.alignment.mode.currentIndex() in ILT_MODES and self.policy.currentIndex()==MONO_REFERENCE
        if self.stale or (not self.alignment.ready and not mono_reference):
            return False
        canvas = QtWidgets.QWidget()
        canvas.setAttribute(QtCore.Qt.WidgetAttribute.WA_DontShowOnScreen)
        canvas.setStyleSheet('background: white;')
        layout = QtWidgets.QVBoxLayout(canvas)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(3)
        try:
            mode = self.alignment.mode.currentIndex()
            title = tr(MODE_NAMES[mode])
            if mode in ILT_MODES:
                title = self.reconstruction_title() + self.ilt_legend()
            elif mode != ORIGINAL:
                overrides, _, profiles = self.alignment.map_data()
                title += f' · {tr("Przetworzone ROI")} {len(overrides) + len(profiles)}'
            heading = QtWidgets.QLabel(title)
            heading.setStyleSheet('color: black; background: white; font-weight: bold;')
            heading.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
            heading.setWordWrap(True)
            layout.addWidget(heading)
            spectrum = light_plot(self.projection)
            map_plot = light_plot(self.plot)
            plots = [(self.projection, spectrum), (self.plot, map_plot)]
            guides = [(spectrum, self.projection_guide), (map_plot, self.ppm_guide), (map_plot, self.diffusion_guide)]
            navy = QtGui.QColor(NAVY)
            for source, plot in plots:
                plot.setFixedHeight(source.height())
            if self.diffusion_button.isChecked():
                diffusion = light_plot(self.diffusion_plot)
                diffusion.hideAxis('left')
                for plot in (diffusion, map_plot):
                    plot.getAxis('bottom').setHeight(self.plot.getAxis('bottom').height())
                diffusion.setFixedSize(self.diffusion_plot.size())
                diffusion.plot(*self.diffusion_trace.getData(), pen=pg.mkPen(navy, width=1.3))
                spacer = QtWidgets.QWidget()
                spacer.setFixedWidth(diffusion.width())
                layout.addWidget(self.plot_row(spectrum, spacer))
                layout.addWidget(self.plot_row(map_plot, diffusion))
                plots.append((self.diffusion_plot, diffusion))
                guides.append((diffusion, self.diffusion_projection_guide))
            else:
                layout.addWidget(spectrum)
                layout.addWidget(map_plot)
            spectrum.plot(*self.trace.getData(), pen=pg.mkPen(navy, width=1.3))
            for x,y in self.unresolved_curves:
                spectrum.plot(x,y,pen=pg.mkPen('#b45309',width=2))
            for lo,hi,d in self.native_guides:
                map_plot.plot([lo,hi],[d,d],pen=pg.mkPen('#7e22ce',width=2,
                    style=QtCore.Qt.PenStyle.DashLine))
            lookup = self.color_table(navy)
            for image in self.images:
                copy = pg.ImageItem(axisOrder='row-major')
                copy.setImage(image.image, levels=image.getLevels(), lut=lookup)
                copy.setTransform(image.transform())
                copy.setOpacity(image.opacity())
                map_plot.addItem(copy)
            pen = pg.mkPen('k', style=QtCore.Qt.PenStyle.DashLine)
            for plot, guide in guides:
                if guide.isVisible():
                    line = pg.InfiniteLine(pos=guide.value(), angle=guide.angle, pen=pen)
                    line.setZValue(10)
                    plot.addItem(line, ignoreBounds=True)
            if self.d_label.isVisible():
                label = pg.TextItem(self.d_label.toPlainText(), color='k', anchor=self.d_label.anchor,
                                    fill=pg.mkBrush('w'))
                label.setPos(self.d_label.pos())
                label.setZValue(11)
                map_plot.addItem(label, ignoreBounds=True)
            canvas.resize(self.graphics.width(), layout.sizeHint().height())
            canvas.show()
            for source, plot in plots:
                x_range, y_range = source.getViewBox().viewRange()
                plot.setRange(xRange=x_range, yRange=y_range, padding=0)
            image = canvas.grab().toImage()
            image.setText('DOSY reconstruction',title)
            image.setDotsPerMeterX(11811)
            image.setDotsPerMeterY(11811)
            return image.save(str(path), 'PNG')
        finally:
            canvas.close()
            canvas.deleteLater()
