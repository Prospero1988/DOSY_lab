"""Compact shared 1D processing controls; polynomial baseline is previewed first."""
from i18n import tr

from copy import deepcopy
from pathlib import Path
import csv
import numpy as np
from PySide6 import QtCore,QtWidgets
from core import baseline_mask, load_profile, save_profile, solvent_center_ppm, solvent_wing_ranges
from jeol import validate_b
from notices import NoticeLabel
from storage import PROCESSING, PROFILE_FILTER, choose_save_file, dialog_directory, source_save_name


class ProcessingPanel(QtWidgets.QDockWidget):
    changed=QtCore.Signal()
    mode_requested=QtCore.Signal(str)
    solvent_wings_requested=QtCore.Signal()
    def __init__(self,parent):
        super().__init__(tr('Processing 1D — wspólny dla warstw'),parent)
        self.controls={};self.regions=[];self._setting=True
        self.source_name = None
        self._dataset = None
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setWidget(scroll)
        widget=QtWidgets.QWidget();scroll.setWidget(widget);layout=QtWidgets.QVBoxLayout(widget)
        form=QtWidgets.QFormLayout();layout.addLayout(form)
        self.reference_target = QtWidgets.QDoubleSpinBox()
        self.reference_target.setRange(-1000., 1000.)
        self.reference_target.setDecimals(6)
        self.reference_target.setSingleStep(.01)
        self.reference_target.setSuffix(' ppm')
        self.reference_button = QtWidgets.QPushButton(tr('Ustaw odniesienie z widma'))
        self.reference_button.setCheckable(True)
        self.reference_button.setEnabled(False)
        self.reference_button.setToolTip(tr('Kliknij sygnał na widmie, aby przypisać mu wpisaną wartość ppm.'))
        self.reference_button.clicked.connect(lambda checked: self.mode_requested.emit('reference' if checked else 'zoom'))
        for key,title,lo,hi,step,dec in [
            ('lb_Hz','LB [Hz]',0,100,.1,3),('gaussian_FWHM_Hz',tr('Gauss: FWHM [Hz]'),0,100,.1,3),
            ('phase0_deg','PH0 [°]',-3600,3600,.5,6),('phase1_deg',tr('PH1 [° / całe SW]'),-100000,100000,10,6),
            ('pivot_ppm','Pivot [ppm]',-1000,1000,.001,6),('ppm_shift',tr('Referencja: przesunięcie ppm'),-100,100,.001,6),
            ('filter_delay_samples',tr('Opóźnienie filtra [próbki]'),0,250,.01,6)]:
            if key == 'ppm_shift':
                form.addRow(tr('Docelowa wartość referencji'), self.reference_target)
                form.addRow('', self.reference_button)
            box=QtWidgets.QDoubleSpinBox();box.setRange(lo,hi);box.setSingleStep(step);box.setDecimals(dec)
            box.valueChanged.connect(self.emit_change);self.controls[key]=box;form.addRow(title,box)
        self.controls['gaussian_FWHM_Hz'].setToolTip(tr('Szerokość FWHM jądra Gaussa w Hz: exp(-(π GB t)²/(4 ln 2)). GB=0 wyłącza funkcję. Można łączyć z LB.'))
        self.controls['phase1_deg'].setToolTip(tr('Stopnie pomiędzy końcami całej akwizycji SW, w kierunku rosnących ppm. Pivot ustala punkt PH0. Bez niezależnych faz warstw.'))
        self.controls['filter_delay_samples'].setToolTip(tr('Wspólna rampa fazowa filtra, oddzielna od ręcznej PH1. 0 wyłącza ją. Dane FID pozostają niezmienione.'))
        zf=QtWidgets.QComboBox()
        for n in (1,2,4,8):zf.addItem(f'{n}×'+(tr(' — bez dopełnienia') if n==1 else ''),n)
        zf.currentIndexChanged.connect(self.emit_change);form.addRow('Zero filling',zf);self.controls['zero_fill_factor']=zf
        ref=QtWidgets.QSpinBox();ref.setRange(1,1);ref.valueChanged.connect(self.emit_change)
        self.controls['reference_layer']=ref;form.addRow(tr('Warstwa referencyjna'),ref)
        crop=QtWidgets.QCheckBox(tr('Przytnij pasmo filtra według JDF'));crop.toggled.connect(self.emit_change)
        self.controls['crop_filter']=crop;layout.addWidget(crop)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        for text,fn in [(tr('Reset faz'),self.reset_phase),(tr('Pivot: klik'),lambda:self.mode_requested.emit('pivot'))]:
            b=QtWidgets.QPushButton(text);b.clicked.connect(fn);row.addWidget(b)
        layout.addWidget(QtWidgets.QLabel(tr('Baza wielomianowa: tylko wskazane obszary')))
        form=QtWidgets.QFormLayout();layout.addLayout(form)
        degree=QtWidgets.QSpinBox();degree.setRange(1,5);degree.valueChanged.connect(self.emit_change)
        self.controls['baseline_degree']=degree;form.addRow(tr('Rząd 1–5'),degree)
        for key,title in [('baseline_per_layer',tr('Osobny wielomian dla każdej warstwy')),
                          ('baseline_preview',tr('Pokaż dopasowaną bazę')),('baseline_enabled',tr('Zatwierdź odejmowanie bazy'))]:
            b=QtWidgets.QCheckBox(title)
            if key=='baseline_enabled':b.toggled.connect(self.update_baseline_mode)
            b.toggled.connect(self.emit_change);self.controls[key]=b;layout.addWidget(b)
        per_layer=self.controls['baseline_per_layer']
        per_layer.setChecked(True)
        per_layer.setToolTip(tr('Przełącza sposób podglądu. Zatwierdzona korekcja zawsze odejmuje osobną bazę od każdej warstwy.'))
        self.show_baseline_regions = QtWidgets.QCheckBox(tr('Pokaż obszary bazy'))
        self.show_baseline_regions.setChecked(True)
        self.show_baseline_regions.setToolTip(tr('Pokaż lub ukryj niebieskie obszary bazy na widmie. Korekcja nadal korzysta z zapisanych obszarów.'))
        layout.addWidget(self.show_baseline_regions)
        self.base_list=QtWidgets.QListWidget();self.base_list.setMaximumHeight(95);layout.addWidget(self.base_list)
        actions = QtWidgets.QGridLayout()
        layout.addLayout(actions)
        self.baseline_button = QtWidgets.QPushButton(tr('Zaznacz bazę'))
        self.baseline_button.setCheckable(True)
        self.baseline_button.setToolTip(tr('Zaznacz fragmenty bez pików, z których program wyznaczy linię bazową. Żółta linia to podgląd.'))
        self.baseline_button.clicked.connect(lambda checked: self.mode_requested.emit('baseline' if checked else 'zoom'))
        actions.addWidget(self.baseline_button, 0, 0)
        remove_region = QtWidgets.QPushButton(tr('Usuń obszar'))
        remove_region.clicked.connect(self.remove_region)
        actions.addWidget(remove_region, 0, 1)
        self.clear_baseline_button = QtWidgets.QPushButton(tr('Wyczyść bazę'))
        self.clear_baseline_button.setToolTip(tr('Usuń wszystkie niebieskie obszary bazy i wyłącz jej odejmowanie.'))
        self.clear_baseline_button.setEnabled(False)
        actions.addWidget(self.clear_baseline_button, 1, 0)
        undo_baseline = QtWidgets.QPushButton(tr('Cofnij korekcję'))
        undo_baseline.clicked.connect(self.undo_baseline)
        actions.addWidget(undo_baseline, 1, 1)
        solvent_panel = QtWidgets.QFrame()
        solvent_panel.setObjectName('solventTailPanel')
        solvent_panel.setStyleSheet('''
            QFrame#solventTailPanel { background:#302619; border:1px solid #805826; border-radius:5px; }
            QFrame#solventTailPanel QLabel, QFrame#solventTailPanel QCheckBox { background:transparent; color:#f9deb3; }
            QFrame#solventTailPanel QCheckBox:disabled { color:#ac9a83; }
        ''')
        solvent_layout = QtWidgets.QVBoxLayout(solvent_panel)
        solvent_layout.setContentsMargins(8, 6, 8, 6)
        solvent_tail = QtWidgets.QCheckBox(tr('Ogon rozpuszczalnika (presaturacja)'))
        solvent_tail.toggled.connect(self.emit_change)
        self.controls['baseline_solvent_tail'] = solvent_tail
        solvent_layout.addWidget(solvent_tail)
        self.solvent_hint = NoticeLabel()
        solvent_layout.addWidget(self.solvent_hint)
        self.solvent_wings_button = QtWidgets.QPushButton(tr('Pokaż skrzydła presaturacji'))
        self.solvent_wings_button.clicked.connect(self.solvent_wings_requested)
        self.solvent_wings_button.setToolTip(tr('Pokaż lub ukryj pomarańczowe pasma wokół IRR. W pasmach zaznacz fragmenty bez pików po obu stronach rozpuszczalnika.'))
        solvent_layout.addWidget(self.solvent_wings_button)
        layout.insertWidget(layout.indexOf(per_layer), solvent_panel)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        for text,fn in [(tr('Zapisz profil'),self.save),(tr('Wczytaj profil'),self.load)]:
            b=QtWidgets.QPushButton(text);b.clicked.connect(fn);row.addWidget(b)
        note=QtWidgets.QLabel(tr('FID → apodyzacja → zero filling → FFT → wspólna faza → zatwierdzona baza → całki.\nLB=GB=0: brak apodyzacji. Niebieskie obszary bazy nie mogą przecinać ROI przy zatwierdzaniu.'))
        note.setWordWrap(True);layout.addWidget(note);layout.addStretch()
        scroll.setMinimumWidth(widget.minimumSizeHint().width() + scroll.verticalScrollBar().sizeHint().width())
        self._setting = False
        self.update_solvent_setup()
        self.hide()

    def emit_change(self,*args):
        if not self._setting:
            self.update_solvent_setup()
            self.changed.emit()

    def set_solvent_context(self, dataset):
        self._dataset = dataset
        self.update_solvent_setup()

    def solvent_ranges(self):
        if self._dataset is None:
            raise ValueError('Wczytaj JDF z presaturacją, aby ustawić korekcję ogona.')
        return solvent_wing_ranges(solvent_center_ppm(self._dataset.metadata, self.controls['ppm_shift'].value()))

    def update_solvent_setup(self):
        control = self.controls['baseline_solvent_tail']
        self.solvent_wings_button.setEnabled(False)
        try:
            wings = self.solvent_ranges()
            self.solvent_wings_button.setEnabled(True)
            shift = self.controls['ppm_shift'].value() - self._dataset.options['ppm_shift']
            ppm = self._dataset.ppm + shift
            center = (wings[0][1] + wings[1][0]) / 2
            baseline_mask(ppm, self.regions, tail_center=center)
        except ValueError as exc:
            control.setEnabled(control.isChecked())
            self.solvent_hint.set_notice(tr(str(exc)), 'warning' if self.solvent_wings_button.isEnabled() else None)
            return
        control.setEnabled(True)
        text = ('Maska gotowa. Włącz ogon i zatwierdź odejmowanie bazy.'
                if not control.isChecked() else
                'Ogon wybrany. Zatwierdź odejmowanie bazy, aby przeliczyć całki.'
                if not self.controls['baseline_enabled'].isChecked() else
                'Ogon wybrany. Sprawdź żółty podgląd na kilku warstwach.')
        self.solvent_hint.set_notice(tr(text))

    def set_options(self,options,nlayers):
        self._setting=True
        for key,control in self.controls.items():
            if key=='reference_layer':control.setRange(1,nlayers);control.setValue(options[key]+1)
            elif isinstance(control,QtWidgets.QCheckBox):control.setChecked(options[key])
            elif isinstance(control,QtWidgets.QComboBox):control.setCurrentIndex(control.findData(options[key]))
            else:control.setValue(options[key])
        self.regions=deepcopy(options['baseline_regions']);self.refresh_regions()
        self.update_baseline_mode(self.controls['baseline_enabled'].isChecked());self._setting=False

    def options(self):
        o={}
        for k,c in self.controls.items():
            if isinstance(c,QtWidgets.QCheckBox):o[k]=c.isChecked()
            elif isinstance(c,QtWidgets.QComboBox):o[k]=c.currentData()
            else:o[k]=c.value()-(1 if k=='reference_layer' else 0)
        o['baseline_regions']=deepcopy(self.regions)
        return o

    def update_baseline_mode(self, enabled):
        per_layer=self.controls['baseline_per_layer']
        if enabled:
            with QtCore.QSignalBlocker(per_layer):
                per_layer.setChecked(True)
        per_layer.setEnabled(not enabled)

    def refresh_regions(self):
        self.base_list.clear()
        for a,b in self.regions:self.base_list.addItem(f'{min(a,b):.5f} … {max(a,b):.5f} ppm')
        self.clear_baseline_button.setEnabled(bool(self.regions))
        self.update_solvent_setup()

    def add_region(self,bounds):
        self.regions.append(sorted(bounds));self.refresh_regions();self.changed.emit()

    def remove_region(self):
        i=self.base_list.currentRow()
        if i<0:return
        if len(self.regions)==1:
            self.clear_regions()
            return
        self.regions.pop(i);self.refresh_regions();self.changed.emit()

    def clear_regions(self):
        self.regions.clear()
        with QtCore.QSignalBlocker(self.controls['baseline_enabled']):
            self.controls['baseline_enabled'].setChecked(False)
        self.update_baseline_mode(False)
        self.refresh_regions()
        self.changed.emit()

    def undo_baseline(self):self.controls['baseline_enabled'].setChecked(False)
    def reset_phase(self):
        self._setting=True
        for key in ('phase0_deg','phase1_deg'):self.controls[key].setValue(0)
        self._setting=False;self.changed.emit()

    def save(self):
        try:
            filename = source_save_name(self.source_name, '.list', 'processing')
            path=choose_save_file(self,tr('Zapisz profil'),PROCESSING,filename,tr(PROFILE_FILTER))
            if path:save_profile(path,self.options())
        except (OSError,ValueError) as exc:
            QtWidgets.QMessageBox.warning(self,tr('Profil processingu'),tr(str(exc)))

    def load(self):
        try:
            path,_=QtWidgets.QFileDialog.getOpenFileName(self,tr('Wczytaj profil'),dialog_directory(PROCESSING),tr(PROFILE_FILTER))
            if path:self.set_options(load_profile(path),self.controls['reference_layer'].maximum());self.changed.emit()
        except (OSError,ValueError) as exc:
            QtWidgets.QMessageBox.warning(self,tr('Profil processingu'),tr(str(exc)))


class BDialog(QtWidgets.QDialog):
    def __init__(self,raw,parent=None):
        super().__init__(parent);self.setWindowTitle(tr('Lista b — wartości w s/m² i jawne źródło'));self.resize(640,600)
        self.raw=raw;self.values=None;layout=QtWidgets.QVBoxLayout(self)
        source=QtWidgets.QLabel(tr(raw.metadata['b_source']['description']));source.setWordWrap(True);layout.addWidget(source)
        self.table=QtWidgets.QTableWidget(len(raw.fid),3);self.table.setHorizontalHeaderLabels([tr('Warstwa akwizycji'),'g [T/m]','b [s/m²]'])
        layout.addWidget(self.table)
        for i,g in enumerate(raw.gradients):
            for j,value in enumerate((str(i+1),f'{g:.8g}' if np.isfinite(g) else '—',f'{raw.b[i]:.12g}' if raw.b is not None else '')):
                item=QtWidgets.QTableWidgetItem(value)
                if j<2:item.setFlags(item.flags()&~QtCore.Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(i,j,item)
        self.reason=QtWidgets.QLineEdit();self.reason.setPlaceholderText(tr('Źródło / opis zmiany, np. niezależna lista z parametrów eksperymentu'))
        layout.addWidget(self.reason)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        load=QtWidgets.QPushButton('Import CSV');load.clicked.connect(self.import_csv);row.addWidget(load)
        apply=QtWidgets.QPushButton(tr('Zastosuj jako ręczną listę'));apply.clicked.connect(self.apply);row.addWidget(apply)
        cancel=QtWidgets.QPushButton(tr('Anuluj'));cancel.clicked.connect(self.reject);row.addWidget(cancel)
        self.error=NoticeLabel();layout.addWidget(self.error)

    def apply(self):
        try:
            if not self.reason.text().strip():raise ValueError(tr('Podaj źródło listy. Nie kalibruj jej oczekiwaną wartością D.'))
            self.values=validate_b([float(self.table.item(i,2).text()) for i in range(self.table.rowCount())],len(self.raw.fid))
            self.accept()
        except Exception as exc:self.error.set_notice(tr(str(exc)), 'error')

    def import_csv(self):
        path,_=QtWidgets.QFileDialog.getOpenFileName(self,tr('Lista b'),'','CSV (*.csv *.txt)')
        if not path:return
        try:
            text=Path(path).read_text(encoding='utf-8-sig');dialect=csv.Sniffer().sniff(text,delimiters=',;\t') if any(c in text for c in ',;\t') else csv.excel
            rows=list(csv.reader(text.splitlines(),dialect));index=0
            if rows[0][0].strip().lower() in ('b','b_s_m2','b_s_per_m2'):rows=rows[1:]
            elif len(rows[0])>1:
                keys=[x.strip().lower() for x in rows[0]]
                index=next(i for i,k in enumerate(keys) if k in ('b','b_s_m2','b_s_per_m2'));rows=rows[1:]
            values=validate_b([float(row[index]) for row in rows if row],len(self.raw.fid))
            for i,v in enumerate(values):self.table.item(i,2).setText(f'{v:.12g}')
            self.reason.setText('Import CSV: '+str(Path(path).resolve()))
            self.error.set_notice('')
        except Exception as exc:self.error.set_notice(tr('Nieprawidłowy CSV: ')+str(exc), 'error')
