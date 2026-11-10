"""Experimental map alignment controls and background calculation."""
from threading import Event
import math

from PySide6 import QtCore, QtWidgets

from adaptive_ilt import adaptive_local_ilt
from dosy_alignment import (alignment_plan, eligible_diffusion, exploratory_profiles, ilt_profiles,
                            individual_map_solution)
from i18n import tr
from ilt_tools import SolveJob
from roi_quality import estimate_spectral_noise


ORIGINAL, ALIGNED, SHARED, ILT, ADAPTIVE, INDIVIDUAL = range(6)
ILT_MODES = (ILT, ADAPTIVE, INDIVIDUAL)
MODE_NAMES = ('Oryginalna mapa DOSY', 'Wyrównanie D · tylko widok', 'Wspólny fit D', 'Global ILT · profile w ROI',
              'Adaptive Local ILT', 'Indywidualne ILT · profile w ROI')
CHECKED = QtCore.Qt.CheckState.Checked
UNCHECKED = QtCore.Qt.CheckState.Unchecked
USER_ROLE = QtCore.Qt.ItemDataRole.UserRole


class AlignmentPanel(QtCore.QObject):
    changed = QtCore.Signal()

    def __init__(self, owner, toolbar, adaptive_toolbar, sidebar):
        super().__init__(owner)
        self.dataset, self.rows, self.snapshot = None, (), None
        self.individual = {}
        self.groups, self.profiles, self.rejected = (), {}, ()
        self.withheld = {}
        self.profile_diagnostics = {}
        self.exploratory, self.exploratory_diagnostics = {}, {}
        self.local_groups = ()
        self.highlighted = set()
        self.highlighted_components = set()
        self.component_membership = {}
        self.manual, self.hidden = set(), set()
        self.job, self.revision = None, 0
        self.cancel = None
        self.ready = True
        self.error = ''
        self.timer = QtCore.QTimer(self, singleShot=True, interval=180)
        self.timer.timeout.connect(self.start)

        toolbar.addWidget(QtWidgets.QLabel(tr('Tryb mapy')))
        self.mode = QtWidgets.QComboBox()
        self.mode.addItems([tr(name) for name in MODE_NAMES])
        self.mode.model().item(ILT).setEnabled(False)
        self.mode.model().item(INDIVIDUAL).setEnabled(False)
        self.mode.currentIndexChanged.connect(self.queue)
        toolbar.addWidget(self.mode)
        toolbar.addWidget(QtWidgets.QLabel(tr('Tolerancja D')))
        self.slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.slider.setRange(10, 300)
        self.slider.setValue(50)
        self.slider.setMaximumWidth(150)
        toolbar.addWidget(self.slider, 1)
        self.tolerance = QtWidgets.QDoubleSpinBox()
        self.tolerance.setRange(1., 30.)
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
        self.adaptive_controls = QtWidgets.QWidget()
        adaptive_toolbar.addWidget(self.adaptive_controls)
        row = QtWidgets.QHBoxLayout(self.adaptive_controls)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(QtWidgets.QLabel(tr('Balans amplitud')))
        self.balance_slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.balance_slider.setRange(0, 300)
        self.balance_slider.setValue(100)
        self.balance_slider.setMaximumWidth(160)
        self.balance_slider.valueChanged.connect(lambda value: self.balance.setValue(10.**(value/100.)))
        row.addWidget(self.balance_slider)
        self.balance = QtWidgets.QDoubleSpinBox()
        self.balance.setRange(1., 1000.)
        self.balance.setDecimals(2)
        self.balance.setValue(10.)
        self.balance.setSuffix(' ×')
        self.balance.setKeyboardTracking(False)
        self.balance.setToolTip(tr('Stosunek największej do słabszej amplitudy. Słabsze ROI otrzymują niższą wagę tylko przy gorszym S/N.'))
        self.balance.valueChanged.connect(self.balance_changed)
        row.addWidget(self.balance)
        row.addWidget(QtWidgets.QLabel(tr('Lokalne λ')))
        self.regularization = QtWidgets.QLineEdit('1e-6')
        self.regularization.setMaximumWidth(90)
        self.regularization.editingFinished.connect(self.queue)
        row.addWidget(self.regularization)
        self.group_count = QtWidgets.QLabel()
        row.addWidget(self.group_count)
        diagnostics = QtWidgets.QPushButton(tr('Diagnostyka grup'))
        diagnostics.clicked.connect(self.show_diagnostics)
        row.addWidget(diagnostics)
        defaults = QtWidgets.QPushButton(tr('Przywróć ustawienia'))
        defaults.clicked.connect(self.reset_defaults)
        row.addWidget(defaults)
        row.addStretch()
        self.tree_label = QtWidgets.QLabel(tr('Grupy wyrównania D'))
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setMaximumHeight(190)
        self.tree.setToolTip(tr('Grupa: widoczność na mapie. ROI: udział we wspólnym D.'))
        self.tree.itemChanged.connect(self.item_changed)
        self.tree.currentItemChanged.connect(self.highlight_group)
        sidebar.addWidget(self.tree_label)
        sidebar.addWidget(self.tree)
        self.update_tree()

    def slider_changed(self, value):
        self.tolerance.setValue(value/10.)

    def tolerance_changed(self, value):
        with QtCore.QSignalBlocker(self.slider):
            self.slider.setValue(round(10*value))
        if self.mode.currentIndex() in (ALIGNED, SHARED, ADAPTIVE):
            self.queue()

    def balance_changed(self, value):
        with QtCore.QSignalBlocker(self.balance_slider):
            self.balance_slider.setValue(round(100*math.log10(value)))
        self.queue()

    def reset_defaults(self):
        self.tolerance.setValue(5.)
        self.balance.setValue(10.)
        self.regularization.setText('1e-6')
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

    def set_individual(self, snapshot):
        if self.individual is snapshot:
            return
        self.individual = snapshot
        self.mode.model().item(INDIVIDUAL).setEnabled(bool(snapshot))
        if self.mode.currentIndex() == INDIVIDUAL:
            self.queue()

    def queue(self, *_, refresh_tree=True):
        self.revision += 1
        if self.cancel is not None:
            self.cancel.set()
        self.error = ''
        self.groups, self.profiles, self.rejected = (), {}, ()
        self.withheld = {}
        self.profile_diagnostics = {}
        self.exploratory, self.exploratory_diagnostics = {}, {}
        self.local_groups = ()
        self.highlighted.clear()
        self.highlighted_components.clear()
        self.component_membership = {}
        self.adaptive_controls.setVisible(self.mode.currentIndex() == ADAPTIVE)
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
        snapshot, individual = self.snapshot, self.individual
        try:
            regularization = float(self.regularization.text().strip().replace(',', '.')) if mode == ADAPTIVE else 0.
        except ValueError:
            self.error = tr('Wpisz prawidłową wartość liczbową.')
            self.changed.emit()
            return
        balance = self.balance.value()
        self.cancel = cancel = Event()

        def compute(_):
            if mode in (ILT, INDIVIDUAL):
                if mode == ILT:
                    solutions = [snapshot] if snapshot is not None else []
                else:
                    solutions = [((sid,),individual_map_solution(result)) for sid,result in individual.items()]
                withheld, diagnostics, explore_diagnostics = {}, {}, {}
                profiles, exploratory, rejected = {}, {}, []
                noise = estimate_spectral_noise(dataset.phased)
                for ids, result in solutions:
                    if cancel.is_set():
                        break
                    args = (dataset.b,dataset.ppm,dataset.spectrum.real,dataset.options['reference_layer'],
                            rows,ids,result)
                    accepted, missing = ilt_profiles(*args, withheld=withheld, diagnostics=diagnostics,
                                                     noise_sigma=noise)
                    recovered, records = exploratory_profiles(*args)
                    profiles.update(accepted)
                    rejected.extend(missing)
                    exploratory.update(recovered)
                    explore_diagnostics.update(records)
                return (), profiles, tuple(rejected), (), {}, withheld, diagnostics, exploratory, explore_diagnostics
            if mode == ADAPTIVE:
                local = adaptive_local_ilt(dataset, rows, tolerance, balance, regularization, excluded, cancel.is_set)
                return ((), local.profiles, local.unresolved, local.groups, local.membership, local.withheld,
                        local.diagnostics, local.exploratory, local.exploratory_diagnostics)
            return alignment_plan(dataset.b, rows, tolerance, excluded, shared=mode == SHARED), {}, (), (), {}, {}, {}, {}, {}

        self.job = SolveJob(self.revision, mode, compute)
        self.job.signals.finished.connect(self.solved)
        QtCore.QThreadPool.globalInstance().start(self.job)

    @QtCore.Slot(int, object, object, str)
    def solved(self, revision, _, payload, error):
        self.job = None
        self.cancel = None
        if revision != self.revision:
            if not self.timer.isActive():
                self.start()
            return
        self.error = error
        if not error:
            (self.groups, self.profiles, self.rejected, self.local_groups,
             self.component_membership, self.withheld, self.profile_diagnostics,
             self.exploratory, self.exploratory_diagnostics) = payload
            self.ready = True
        self.update_tree()
        self.changed.emit()

    def active_group_d(self, group):
        mode = self.mode.currentIndex()
        if mode == ALIGNED:
            return group.representative
        return group.shared if not group.incompatible else None

    def map_data(self, exploratory=False):
        overrides = {}
        hidden = set()
        for group in self.groups:
            d = self.active_group_d(group)
            if d is not None:
                overrides.update(dict.fromkeys(group.ids, d))
            if group.ids in self.hidden:
                hidden.update(group.ids)
        for group in self.local_groups:
            if group.ids in self.hidden:
                hidden.update(group.ids)
        hidden_components = {item for ids in self.hidden for item in self.component_membership.get(ids, ())}
        hidden.update(self.manual if self.mode.currentIndex() == ADAPTIVE else ())
        profiles = self.exploratory if exploratory else self.profiles
        if hidden_components:
            profiles = {sid: tuple((d, y) for d, y in components if (sid, d) not in hidden_components)
                        for sid, components in profiles.items()}
        return overrides, hidden, profiles

    def profile_detail(self, sid):
        record = self.profile_diagnostics.get(sid)
        exploratory = self.exploratory_diagnostics.get(sid)
        if record is None and exploratory is None:
            return tr('ROI nie było wybrane do tego obliczenia ILT.')
        reasons = {
            'assigned':'Przypisanie ILT zaakceptowane przez kryteria rekonstrukcji.',
            'discarded_component':'Pominięta składowa nie może zostać przypisana do pozostałego maksimum.',
            'spectral_mismatch':'Profil widmowy nie zgadza się z zanikiem.',
            'profile_dependent_diffusion':'Alternatywne odczyty kształtu lub fazy dają istotnie różne D dla tego ROI.',
            'peak_mismatch':'Zanik ROI nie zgadza się z D maksimum, mimo możliwego dopasowania całego przedziału ILT.',
            'integral_mismatch':'Przewidywane całki ILT nie zgadzają się z zanikiem ROI.',
            'boundary_peak':'Maksimum na granicy siatki D.',
            'no_component':'Brak dodatniej składowej powyżej progu udziału.',
            'no_positive_reconstructed_profile':'Odtworzone składowe nie mają dodatniego pola na warstwie referencyjnej.',
            'insufficient_profile':'Za mało punktów widmowych do rekonstrukcji.',
            'insufficient_spectral_snr':'Za mały S/N profilu referencyjnego.',
            'ill_conditioned_basis':'Za mało niezależnych danych do odtworzenia składowych.',
            'spectral_solver_failed':'Nie udało się odtworzyć profilu widmowego.',
            'solver_not_converged':'Optymalizacja ILT nie osiągnęła zbieżności.',
            'solver_failed':'Nie udało się przeprowadzić optymalizacji ILT.',
            'insufficient_attenuation_information':'Zanik nie spełnia kryteriów wejściowych AILT.',
            'no_group_or_independent_split_evidence':'Brak zgodnej grupy i podstaw do niezależnego rozdzielenia.',
            'no_stable_independent_split':'Brak stabilnego niezależnego rozdzielenia ROI.',
            'unexplained_reference':'Część profilu referencyjnego nie ma wyjaśnienia w składowych.',
        }
        record = record or {}
        lines = []
        if record:
            lines.append(tr('Filtr QC') + ' · ' + tr(reasons.get(record['reason'],record['reason'])))
        if exploratory is not None:
            reason = exploratory['reason']
            label = 'Profil pokazany w widoku eksploracyjnym.' if reason == 'assigned' else reasons.get(reason,reason)
            lines.append(tr('Eksploracyjna') + ' · ' + tr(label))
            if exploratory.get('independent'):
                lines.append(tr('Niezależne ILT dla ROI poza zgodną grupą.'))
            if exploratory.get('candidate_D'):
                peaks = ', '.join(f'{d:.4e}' for d in exploratory['candidate_D'])
                lines.append(f"{tr('Maksima przedziałów ILT')}: {peaks} m²/s")
            if exploratory.get('boundary_peak'):
                lines.append(tr(reasons['boundary_peak']))
            if exploratory.get('unexplained_positive_fraction') is not None:
                label = tr('Niewyjaśnione dodatnie pole referencji')
                lines.append(f"{tr('Eksploracyjna')} · {label}: {exploratory['unexplained_positive_fraction']:.2%}")
            if exploratory.get('relative_residual') is not None:
                lines.append(f"{tr('Eksploracyjna względna reszta widmowa')}: {exploratory['relative_residual']:.2%}")
        if record.get('omitted_fraction'):
            lines.append(f"{tr('Pominięta masa modelu na warstwie referencyjnej')}: {record['omitted_fraction']:.2%}")
        for diffusion, fraction in record.get('omitted_components',()):
            lines.append(f"{tr('Pominięta składowa')}: {diffusion:.4e} m²/s · {fraction:.2%}")
        if record.get('unexplained_positive_fraction') is not None:
            lines.append(f"{tr('Filtr QC')} · {tr('Niewyjaśnione dodatnie pole referencji')}: "
                         f"{record['unexplained_positive_fraction']:.2%}")
        for key, limit, label in (
            ('peak_rms','peak_limit','RMS dla D maksimum'),
            ('integral_rms','integral_limit','RMS całek modelu ILT'),
            ('excess_residual','spectral_limit','Reszta widmowa ponad szum'),
        ):
            if record.get(key) is not None:
                lines.append(f'{tr("Filtr QC")} · {tr(label)}: {record[key]:.2%} / {record[limit]:.2%}')
        return '\n'.join(lines)

    def highlight_group(self, item, _):
        entry = item.data(0, USER_ROLE) if item is not None else None
        self.highlighted = set(entry[1]) if entry is not None and entry[0] == 'group' else set()
        self.highlighted_components = (set(self.component_membership.get(entry[1], ())) if self.highlighted else set())
        self.changed.emit()

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
        adaptive = self.mode.currentIndex() == ADAPTIVE
        self.adaptive_controls.setVisible(adaptive)
        self.group_count.setText(tr('Grupy lokalne') + f' · {sum(not g.independent for g in self.local_groups)}')
        show = self.mode.currentIndex() in (ALIGNED, SHARED, ADAPTIVE)
        self.tree_label.setVisible(show)
        self.tree.setVisible(show)
        with QtCore.QSignalBlocker(self.tree):
            self.tree.clear()
            if not show:
                return
            by_id = {row['id']: row for row in self.rows}
            assigned = set()
            for index, group in enumerate(self.local_groups, 1):
                assigned.update(group.ids)
                root = self.add_group(f'L{index} · {self.local_label(group)}', group.ids)
                for sid in group.ids:
                    self.add_roi(root, by_id[sid], group.dominant)
                for sid, d in self.component_membership.get(group.ids, ()):
                    child = QtWidgets.QTreeWidgetItem([f'{sid} · {d:.3e} · {tr("składowa")}'])
                    child.setToolTip(0, tr('Składowa rozdzielona niezależnie. Dopasowanie do grupy dotyczy widoku.'))
                    root.addChild(child)
            for group in self.groups:
                assigned.update(group.ids)
                d = self.active_group_d(group)
                label = f'D = {d:.4e} m²/s' if d is not None else tr('⚠ Wspólny D odrzucony')
                root = self.add_group(label, group.ids)
                if group.incompatible:
                    root.setToolTip(0, tr('Niezgodne zaniki') + ' · ' + ', '.join(group.incompatible))
                diagnostics = dict(zip(group.ids, group.diagnostics))
                for sid in group.ids:
                    self.add_roi(root, by_id[sid], group.shared or d, diagnostics.get(sid))
            outside = [row for row in self.rows if row['id'] not in assigned
                       and (row['id'] in self.manual or eligible_diffusion(row) is not None)]
            if outside:
                root = QtWidgets.QTreeWidgetItem([tr('ROI poza grupami')])
                self.tree.addTopLevelItem(root)
                root.setExpanded(True)
                for row in outside:
                    self.add_roi(root, row)

    def add_group(self, label, ids):
        root = QtWidgets.QTreeWidgetItem([label])
        root.setData(0, USER_ROLE, ('group', ids))
        root.setCheckState(0, UNCHECKED if ids in self.hidden else CHECKED)
        self.tree.addTopLevelItem(root)
        root.setExpanded(True)
        return root

    @staticmethod
    def local_label(group):
        if group.independent:
            return tr('Niezależne rozdzielenie ROI')
        if not group.converged:
            return tr('Brak zbieżności')
        label = f'D = {group.dominant:.3e}' if group.dominant is not None else 'D = —'
        return label + (f' · ⚠ {len(group.rejected)}' if group.rejected else '')

    def show_diagnostics(self):
        dialog = QtWidgets.QDialog(self.parent())
        dialog.setWindowTitle(tr('Diagnostyka grup'))
        dialog.resize(1000, 380)
        layout = QtWidgets.QVBoxLayout(dialog)
        table = QtWidgets.QTableWidget(len(self.local_groups), 8)
        table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setHorizontalHeaderLabels([tr('Grupa'), 'ROI', tr('Zakres D [m²/s]'), tr('D dominujące [m²/s]'),
                                         tr('Stosunek amplitud'), tr('Status'), 'RMS / Iₛ', tr('Nierozdzielone')])
        table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        for index, group in enumerate(self.local_groups):
            status = ('Niezależne rozdzielenie ROI' if group.independent
                      else 'Zbieżny' if group.converged else 'Brak zbieżności')
            values = (f'L{index+1}', str(len(group.ids)), ' – '.join(f'{d:.3e}' for d in group.d_range),
                      f'{group.dominant:.3e}' if group.dominant else '—', f'{group.amplitude_ratio:.3g}×',
                      tr(status), f'{group.rms:.4g}' if group.rms is not None else '—', str(len(group.rejected)))
            for column, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(value)
                members = ', '.join(group.rejected if column == 7 else group.ids)
                if column == 1 and group.ids in self.component_membership:
                    members += '\n' + ', '.join(f'{sid} · {d:.3e}' for sid, d in self.component_membership[group.ids])
                item.setToolTip(members)
                table.setItem(index, column, item)
        layout.addWidget(table)
        if self.rejected:
            note = QtWidgets.QLabel(tr('ROI bez przypisania D w tej rekonstrukcji') + ' · ' + ', '.join(self.rejected))
            note.setWordWrap(True)
            layout.addWidget(note)
        close = QtWidgets.QPushButton(tr('Zamknij'))
        close.clicked.connect(dialog.close)
        layout.addWidget(close)
        dialog.exec()

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
