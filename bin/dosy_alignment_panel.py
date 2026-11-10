"""Experimental map alignment controls and background calculation."""
from PySide6 import QtCore, QtWidgets

from dosy_alignment import alignment_plan, eligible_diffusion, ilt_profiles
from i18n import tr
from ilt_tools import SolveJob


ORIGINAL, ALIGNED, SHARED, ILT = range(4)
MODE_NAMES = ('Oryginalna mapa DOSY', 'Wyrównanie D · tylko widok', 'Wspólny fit D', 'Global ILT · profile w ROI')
CHECKED = QtCore.Qt.CheckState.Checked
UNCHECKED = QtCore.Qt.CheckState.Unchecked
USER_ROLE = QtCore.Qt.ItemDataRole.UserRole


class AlignmentPanel(QtCore.QObject):
    changed = QtCore.Signal()

    def __init__(self, owner, toolbar, sidebar):
        super().__init__(owner)
        self.dataset, self.rows, self.snapshot = None, (), None
        self.groups, self.profiles, self.rejected = (), {}, ()
        self.manual, self.hidden = set(), set()
        self.job, self.revision = None, 0
        self.ready = True
        self.error = ''
        self.timer = QtCore.QTimer(self, singleShot=True, interval=180)
        self.timer.timeout.connect(self.start)

        toolbar.addWidget(QtWidgets.QLabel(tr('Tryb mapy')))
        self.mode = QtWidgets.QComboBox()
        self.mode.addItems([tr(name) for name in MODE_NAMES])
        self.mode.model().item(ILT).setEnabled(False)
        self.mode.currentIndexChanged.connect(self.queue)
        toolbar.addWidget(self.mode)
        toolbar.addWidget(QtWidgets.QLabel(tr('Tolerancja D')))
        self.slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.slider.setRange(10, 200)
        self.slider.setValue(50)
        self.slider.setMaximumWidth(150)
        toolbar.addWidget(self.slider, 1)
        self.tolerance = QtWidgets.QDoubleSpinBox()
        self.tolerance.setRange(1., 20.)
        self.tolerance.setDecimals(2)
        self.tolerance.setSingleStep(.1)
        self.tolerance.setSuffix(' %')
        self.tolerance.setValue(5.)
        self.tolerance.setKeyboardTracking(False)
        self.slider.valueChanged.connect(self.slider_changed)
        self.tolerance.valueChanged.connect(self.tolerance_changed)
        toolbar.addWidget(self.tolerance)
        apply = QtWidgets.QPushButton(tr('Przelicz wyrównanie'))
        apply.clicked.connect(self.recalculate)
        toolbar.addWidget(apply)
        reset = QtWidgets.QPushButton(tr('Wróć do oryginału'))
        reset.clicked.connect(lambda: (self.mode.setCurrentIndex(ORIGINAL), owner.reset_view()))
        toolbar.addWidget(reset)
        self.tree_label = QtWidgets.QLabel(tr('Grupy wyrównania D'))
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setMaximumHeight(190)
        self.tree.setToolTip(tr('Grupa: widoczność na mapie. ROI: udział we wspólnym D.'))
        self.tree.itemChanged.connect(self.item_changed)
        sidebar.addWidget(self.tree_label)
        sidebar.addWidget(self.tree)
        self.update_tree()

    def slider_changed(self, value):
        self.tolerance.setValue(value/10.)

    def tolerance_changed(self, value):
        with QtCore.QSignalBlocker(self.slider):
            self.slider.setValue(round(10*value))
        if self.mode.currentIndex() in (ALIGNED, SHARED):
            self.queue()

    def set_data(self, dataset, rows, new_experiment):
        if self.dataset is dataset and self.rows is rows:
            return
        self.dataset, self.rows = dataset, rows
        if new_experiment:
            self.manual.clear()
            self.hidden.clear()
            self.mode.setCurrentIndex(ORIGINAL)
        self.queue()

    def set_snapshot(self, snapshot):
        if self.snapshot is snapshot:
            return
        self.snapshot = snapshot
        self.mode.model().item(ILT).setEnabled(snapshot is not None)
        if self.mode.currentIndex() == ILT:
            if snapshot is None:
                self.mode.setCurrentIndex(ORIGINAL)
            else:
                self.queue()

    def queue(self, *_, refresh_tree=True):
        self.revision += 1
        self.error = ''
        self.groups, self.profiles, self.rejected = (), {}, ()
        self.ready = self.mode.currentIndex() == ORIGINAL
        self.timer.stop()
        if not self.ready and self.dataset is not None:
            self.timer.start()
        if refresh_tree:
            self.update_tree()
        self.changed.emit()

    def recalculate(self):
        self.queue()
        self.timer.stop()
        self.start()

    def start(self):
        mode = self.mode.currentIndex()
        if self.job is not None or self.dataset is None or mode == ORIGINAL:
            return
        dataset, rows = self.dataset, self.rows
        excluded, tolerance = set(self.manual), self.tolerance.value()
        snapshot = self.snapshot

        def compute(_):
            if mode == ILT and snapshot is not None:
                ids, result = snapshot
                profiles, rejected = ilt_profiles(dataset.b, dataset.ppm, dataset.spectrum.real,
                                                  dataset.options['reference_layer'], rows, ids, result)
                return (), profiles, rejected
            return alignment_plan(dataset.b, rows, tolerance, excluded, shared=mode == SHARED), {}, ()

        self.job = SolveJob(self.revision, mode, compute)
        self.job.signals.finished.connect(self.solved)
        QtCore.QThreadPool.globalInstance().start(self.job)

    @QtCore.Slot(int, object, object, str)
    def solved(self, revision, _, payload, error):
        self.job = None
        if revision != self.revision:
            if not self.timer.isActive():
                self.start()
            return
        self.error = error
        if not error:
            self.groups, self.profiles, self.rejected = payload
            self.ready = True
        self.update_tree()
        self.changed.emit()

    def active_group_d(self, group):
        mode = self.mode.currentIndex()
        if mode == ALIGNED:
            return group.representative
        return group.shared if not group.incompatible else None

    def map_data(self):
        overrides = {}
        hidden = set()
        for group in self.groups:
            d = self.active_group_d(group)
            if d is not None:
                overrides.update(dict.fromkeys(group.ids, d))
            if group.ids in self.hidden:
                hidden.update(group.ids)
        return overrides, hidden, self.profiles

    def item_changed(self, item, column):
        info = item.data(column, USER_ROLE)
        if info is None:
            return
        kind, key = info
        selected = self.hidden if kind == 'group' else self.manual
        if item.checkState(column) == CHECKED:
            selected.discard(key)
        else:
            selected.add(key)
        if kind == 'group':
            self.changed.emit()
        else:
            self.queue(refresh_tree=False)

    def update_tree(self):
        show = self.mode.currentIndex() in (ALIGNED, SHARED)
        self.tree_label.setVisible(show)
        self.tree.setVisible(show)
        with QtCore.QSignalBlocker(self.tree):
            self.tree.clear()
            if not show:
                return
            by_id = {row['id']: row for row in self.rows}
            assigned = set()
            for group in self.groups:
                assigned.update(group.ids)
                d = self.active_group_d(group)
                label = f'D = {d:.4e} m²/s' if d is not None else tr('⚠ Wspólny D odrzucony')
                root = QtWidgets.QTreeWidgetItem([label])
                root.setData(0, USER_ROLE, ('group', group.ids))
                root.setCheckState(0, UNCHECKED if group.ids in self.hidden else CHECKED)
                if group.incompatible:
                    root.setToolTip(0, tr('Niezgodne zaniki') + ' · ' + ', '.join(group.incompatible))
                self.tree.addTopLevelItem(root)
                root.setExpanded(True)
                diagnostics = dict(zip(group.ids, group.diagnostics))
                for sid in group.ids:
                    self.add_roi(root, by_id[sid], group.shared or d, diagnostics.get(sid))
            outside = [row for row in self.rows if row['id'] not in assigned and eligible_diffusion(row) is not None]
            if outside:
                root = QtWidgets.QTreeWidgetItem([tr('ROI poza grupami')])
                self.tree.addTopLevelItem(root)
                root.setExpanded(True)
                for row in outside:
                    self.add_roi(root, row)

    def add_roi(self, root, row, shared=None, rms=None):
        sid = row['id']
        original = row['fit']['D_m2_s']
        text = f'{sid} · {original:.3e}'
        if shared is not None:
            text += f' → {shared:.3e}'
        item = QtWidgets.QTreeWidgetItem([text])
        item.setData(0, USER_ROLE, ('roi', sid))
        item.setCheckState(0, UNCHECKED if sid in self.manual else CHECKED)
        tip = tr('D oryginalne → D grupy [m²/s]')
        if rms is not None:
            tip += f'\nRMS: {rms[0]:.4g} → {rms[1]:.4g}'
        item.setToolTip(0, tip)
        root.addChild(item)
