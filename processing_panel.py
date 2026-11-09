"""Compact shared 1D processing controls; polynomial baseline is previewed first."""
from copy import deepcopy
from pathlib import Path
import csv
import numpy as np
from PySide6 import QtCore,QtWidgets
from core import save_profile,load_profile
from jeol import validate_b


class ProcessingPanel(QtWidgets.QDockWidget):
    changed=QtCore.Signal()
    mode_requested=QtCore.Signal(str)
    def __init__(self,parent):
        super().__init__('Processing 1D — wspólny dla warstw',parent)
        self.controls={};self.regions=[];self._setting=False
        widget=QtWidgets.QWidget();self.setWidget(widget);layout=QtWidgets.QVBoxLayout(widget)
        form=QtWidgets.QFormLayout();layout.addLayout(form)
        for key,title,lo,hi,step,dec in [
            ('lb_Hz','LB [Hz]',0,100,.1,3),('gaussian_FWHM_Hz','Gauss: FWHM [Hz]',0,100,.1,3),
            ('phase0_deg','PH0 [°]',-3600,3600,.5,6),('phase1_deg','PH1 [° / całe SW]',-100000,100000,10,6),
            ('pivot_ppm','Pivot [ppm]',-1000,1000,.001,6),('ppm_shift','Referencja: przesunięcie ppm',-100,100,.001,6),
            ('filter_delay_samples','Opóźnienie filtra [próbki]',0,250,.01,6)]:
            box=QtWidgets.QDoubleSpinBox();box.setRange(lo,hi);box.setSingleStep(step);box.setDecimals(dec)
            box.valueChanged.connect(self.emit_change);self.controls[key]=box;form.addRow(title,box)
        self.controls['gaussian_FWHM_Hz'].setToolTip('Szerokość FWHM jądra Gaussa w Hz: exp(-(π GB t)²/(4 ln 2)). GB=0 wyłącza funkcję. Można łączyć z LB.')
        self.controls['phase1_deg'].setToolTip('Stopnie pomiędzy końcami całej akwizycji SW, w kierunku rosnących ppm. Pivot ustala punkt PH0. Bez niezależnych faz warstw.')
        self.controls['filter_delay_samples'].setToolTip('Wspólna rampa fazowa filtra, oddzielna od ręcznej PH1. 0 wyłącza ją. Dane FID pozostają niezmienione.')
        zf=QtWidgets.QComboBox()
        for n in (1,2,4,8):zf.addItem(f'{n}×'+(' — bez dopełnienia' if n==1 else ''),n)
        zf.currentIndexChanged.connect(self.emit_change);form.addRow('Zero filling',zf);self.controls['zero_fill_factor']=zf
        ref=QtWidgets.QSpinBox();ref.setRange(1,1);ref.valueChanged.connect(self.emit_change)
        self.controls['reference_layer']=ref;form.addRow('Warstwa referencyjna',ref)
        crop=QtWidgets.QCheckBox('Przytnij pasmo filtra według JDF');crop.toggled.connect(self.emit_change)
        self.controls['crop_filter']=crop;layout.addWidget(crop)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        for text,fn in [('Reset faz',self.reset_phase),('Pivot: klik',lambda:self.mode_requested.emit('pivot'))]:
            b=QtWidgets.QPushButton(text);b.clicked.connect(fn);row.addWidget(b)
        layout.addWidget(QtWidgets.QLabel('Baza wielomianowa: tylko wskazane obszary'))
        form=QtWidgets.QFormLayout();layout.addLayout(form)
        degree=QtWidgets.QSpinBox();degree.setRange(1,5);degree.valueChanged.connect(self.emit_change)
        self.controls['baseline_degree']=degree;form.addRow('Rząd 1–5',degree)
        for key,title in [('baseline_per_layer','Osobny wielomian dla każdej warstwy'),
                          ('baseline_preview','Pokaż dopasowaną bazę'),('baseline_enabled','Zatwierdź odejmowanie bazy')]:
            b=QtWidgets.QCheckBox(title)
            if key=='baseline_enabled':b.toggled.connect(self.update_baseline_mode)
            b.toggled.connect(self.emit_change);self.controls[key]=b;layout.addWidget(b)
        per_layer=self.controls['baseline_per_layer']
        per_layer.setChecked(True)
        per_layer.setToolTip('Przełącza sposób podglądu. Zatwierdzona korekcja zawsze odejmuje osobną bazę od każdej warstwy.')
        self.base_list=QtWidgets.QListWidget();self.base_list.setMaximumHeight(95);layout.addWidget(self.base_list)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        for text,fn in [('Zaznacz bazę',lambda:self.mode_requested.emit('baseline')),('Usuń obszar',self.remove_region),('Cofnij korekcję',self.undo_baseline)]:
            b=QtWidgets.QPushButton(text);b.clicked.connect(fn);row.addWidget(b)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        for text,fn in [('Zapisz profil',self.save),('Wczytaj profil',self.load)]:
            b=QtWidgets.QPushButton(text);b.clicked.connect(fn);row.addWidget(b)
        note=QtWidgets.QLabel('FID → apodyzacja → zero filling → FFT → wspólna faza → zatwierdzona baza → całki.\nLB=GB=0: brak apodyzacji. Niebieskie obszary bazy nie mogą przecinać ROI przy zatwierdzaniu.')
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
        path,_=QtWidgets.QFileDialog.getSaveFileName(self,'Zapisz processing','processing.json','JSON (*.json)')
        if path:save_profile(path,self.options())

    def load(self):
        path,_=QtWidgets.QFileDialog.getOpenFileName(self,'Wczytaj processing','','JSON (*.json)')
        if path:
            try:self.set_options(load_profile(path),self.controls['reference_layer'].maximum());self.changed.emit()
            except Exception as e:QtWidgets.QMessageBox.warning(self,'Profil processingu',str(e))


class BDialog(QtWidgets.QDialog):
    def __init__(self,raw,parent=None):
        super().__init__(parent);self.setWindowTitle('Lista b — wartości w s/m² i jawne źródło');self.resize(640,600)
        self.raw=raw;self.values=None;layout=QtWidgets.QVBoxLayout(self)
        source=QtWidgets.QLabel(raw.metadata['b_source']['description']);source.setWordWrap(True);layout.addWidget(source)
        self.table=QtWidgets.QTableWidget(len(raw.fid),3);self.table.setHorizontalHeaderLabels(['Warstwa akwizycji','g [T/m]','b [s/m²]'])
        layout.addWidget(self.table)
        for i,g in enumerate(raw.gradients):
            for j,value in enumerate((str(i+1),f'{g:.8g}' if np.isfinite(g) else '—',f'{raw.b[i]:.12g}' if raw.b is not None else '')):
                item=QtWidgets.QTableWidgetItem(value)
                if j<2:item.setFlags(item.flags()&~QtCore.Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(i,j,item)
        self.reason=QtWidgets.QLineEdit();self.reason.setPlaceholderText('Źródło / opis zmiany, np. niezależna lista z parametrów eksperymentu')
        layout.addWidget(self.reason)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        load=QtWidgets.QPushButton('Import CSV');load.clicked.connect(self.import_csv);row.addWidget(load)
        apply=QtWidgets.QPushButton('Zastosuj jako ręczną listę');apply.clicked.connect(self.apply);row.addWidget(apply)
        cancel=QtWidgets.QPushButton('Anuluj');cancel.clicked.connect(self.reject);row.addWidget(cancel)
        self.error=QtWidgets.QLabel('');self.error.setWordWrap(True);layout.addWidget(self.error)

    def apply(self):
        try:
            if not self.reason.text().strip():raise ValueError('Podaj źródło listy. Nie kalibruj jej oczekiwaną wartością D.')
            self.values=validate_b([float(self.table.item(i,2).text()) for i in range(self.table.rowCount())],len(self.raw.fid))
            self.accept()
        except Exception as exc:self.error.setText(str(exc))

    def import_csv(self):
        path,_=QtWidgets.QFileDialog.getOpenFileName(self,'Lista b','','CSV (*.csv *.txt)')
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
        except Exception as exc:self.error.setText('Nieprawidłowy CSV: '+str(exc))
