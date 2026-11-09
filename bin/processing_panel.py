"""Compact shared 1D processing controls; polynomial baseline is previewed first."""
from i18n import tr

from copy import deepcopy
from pathlib import Path
import csv
import numpy as np
from PySide6 import QtCore,QtWidgets
from core import save_profile,load_profile
from jeol import validate_b
from storage import PROCESSING, PROFILE_FILTER, choose_save_file, dialog_directory, source_save_name


class ProcessingPanel(QtWidgets.QDockWidget):
    changed=QtCore.Signal()
    mode_requested=QtCore.Signal(str)
    def __init__(self,parent):
        super().__init__(tr('Processing 1D — wspólny dla warstw'),parent)
        self.controls={};self.regions=[];self._setting=False
        self.source_name = None
        widget=QtWidgets.QWidget();self.setWidget(widget);layout=QtWidgets.QVBoxLayout(widget)
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
        self.base_list=QtWidgets.QListWidget();self.base_list.setMaximumHeight(95);layout.addWidget(self.base_list)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        for text,fn in [(tr('Zaznacz bazę'),lambda:self.mode_requested.emit('baseline')),(tr('Usuń obszar'),self.remove_region),(tr('Cofnij korekcję'),self.undo_baseline)]:
            b=QtWidgets.QPushButton(text);b.clicked.connect(fn);row.addWidget(b)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        for text,fn in [(tr('Zapisz profil'),self.save),(tr('Wczytaj profil'),self.load)]:
            b=QtWidgets.QPushButton(text);b.clicked.connect(fn);row.addWidget(b)
        note=QtWidgets.QLabel(tr('FID → apodyzacja → zero filling → FFT → wspólna faza → zatwierdzona baza → całki.\nLB=GB=0: brak apodyzacji. Niebieskie obszary bazy nie mogą przecinać ROI przy zatwierdzaniu.'))
        note.setWordWrap(True);layout.addWidget(note);layout.addStretch()
        self.hide()

    def emit_change(self,*args):
        if not self._setting:self.changed.emit()

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
        self.error=QtWidgets.QLabel('');self.error.setWordWrap(True);layout.addWidget(self.error)

    def apply(self):
        try:
            if not self.reason.text().strip():raise ValueError(tr('Podaj źródło listy. Nie kalibruj jej oczekiwaną wartością D.'))
            self.values=validate_b([float(self.table.item(i,2).text()) for i in range(self.table.rowCount())],len(self.raw.fid))
            self.accept()
        except Exception as exc:self.error.setText(tr(str(exc)))

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
        except Exception as exc:self.error.setText(tr('Nieprawidłowy CSV: ')+str(exc))
