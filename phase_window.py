"""Manual shared phasing: cached FFT preview, mouse gestures, explicit commit."""
from copy import deepcopy
import numpy as np
from PySide6 import QtCore, QtWidgets
import pyqtgraph as pg
from core import phase_radians
from spectrum_view import PhaseBox
from spectrum_settings import SpectrumAppearance


class PhaseWindow(QtWidgets.QDialog):
    applied = QtCore.Signal(dict)

    def __init__(self, dataset, parent=None, appearance=SpectrumAppearance()):
        super().__init__(parent)
        self.setWindowTitle('DOSY Lab · Interaktywne fazowanie wszystkich warstw')
        self.setWindowModality(QtCore.Qt.WindowModality.WindowModal)
        self.resize(1350,850)
        if parent is not None:self.setStyleSheet(parent.styleSheet())
        self.dataset=dataset
        self.appearance = appearance
        self.initial=deepcopy(dataset.options)
        self.previous_pivot=self.initial['pivot_ppm']
        # Undo only the old phase, preserving the exact common apodization/FFT.
        self.unphased=dataset.phased*np.exp(-1j*phase_radians(dataset.ppm,dataset.metadata,self.initial))[None,:]
        self.preview=None
        self.last_drag_dy = 0.
        layout=QtWidgets.QVBoxLayout(self)
        help_label=QtWidgets.QLabel(
            '<b>Lewy + góra/dół</b> — PH0 &nbsp;·&nbsp; <b>Prawy + góra/dół</b> — PH1 &nbsp;·&nbsp; '
            '<b>Shift + przeciąganie</b> — precyzyjne dostrajanie<br>'
            '<b>Lewy klik</b> — pivot &nbsp;·&nbsp; <b>Środkowy + ruch</b> — przesuwanie &nbsp;·&nbsp; '
            '<b>Rolka</b> — wysokość widma')
        help_label.setWordWrap(True);layout.addWidget(help_label)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        self.controls={}
        for key,label,value,low,high,step,decimals in [
            ('phase0_deg','PH0 [°]',self.initial['phase0_deg'],-3600,3600,.5,6),
            ('phase1_deg','PH1 [° / SW]',self.initial['phase1_deg'],-100000,100000,10,6),
            ('pivot_ppm','Pivot [ppm]',self.initial['pivot_ppm'],-1000,1000,.001,6)]:
            row.addWidget(QtWidgets.QLabel(label))
            control=QtWidgets.QDoubleSpinBox();control.setDecimals(decimals);control.setRange(low,high)
            control.setSingleStep(step);control.setValue(value);self.controls[key]=control
            control.valueChanged.connect(self.pivot_changed if key=='pivot_ppm' else self.refresh);row.addWidget(control)
        row.addWidget(QtWidgets.QLabel('Podgląd:'))
        self.layer=QtWidgets.QComboBox()
        for i in range(len(dataset.raw.fid)):
            b='brak b' if dataset.b is None else f'b={dataset.b[i]:.3e}'
            self.layer.addItem(f'{i+1}: {b}',i)
        self.layer.setCurrentIndex(self.initial['reference_layer']);self.layer.currentIndexChanged.connect(self.draw)
        row.addWidget(self.layer,1)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        row.addWidget(QtWidgets.QLabel('Czułość PH0 [°/px]:'))
        self.ph0_sensitivity=self.sensitivity(row,[.05,.1,.5,1,2],.5)
        row.addWidget(QtWidgets.QLabel('PH1 [°/px]:'))
        self.ph1_sensitivity=self.sensitivity(row,[.5,2,10,50,100],10)
        self.imaginary=QtWidgets.QCheckBox('Pokaż Im');self.imaginary.toggled.connect(self.draw);row.addWidget(self.imaginary)
        self.before=QtWidgets.QCheckBox('Pokaż fazę sprzed edycji');self.before.toggled.connect(self.draw);row.addWidget(self.before)
        self.zoom_button=QtWidgets.QPushButton('Zoom zakresu');self.zoom_button.setCheckable(True)
        self.zoom_button.setToolTip('Włącz, aby przeciągnięciem powiększyć zakres. Wyłącz, aby wrócić do gestów PH0/PH1.')
        row.addWidget(self.zoom_button)
        self.pan_button = QtWidgets.QPushButton('Przesuń widmo')
        self.pan_button.setCheckable(True)
        self.pan_button.setToolTip('Przeciągnij widmo lewym przyciskiem w pionie lub poziomie. Wyłącz, aby wrócić do fazowania.')
        row.addWidget(self.pan_button)
        self.zoom_button.toggled.connect(lambda checked:self.set_navigation(self.pan_button,checked))
        self.pan_button.toggled.connect(lambda checked:self.set_navigation(self.zoom_button,checked))
        row.addStretch()
        for text,callback in [('Cały zakres',lambda:self.plot.autoRange()),('Pełny ekran',self.toggle_fullscreen),('Reset PH0/PH1',self.reset_phase)]:
            button=QtWidgets.QPushButton(text);button.clicked.connect(callback);row.addWidget(button)
        self.view=PhaseBox();self.plot=pg.PlotWidget(viewBox=self.view);layout.addWidget(self.plot,1)
        self.view.bind_selection_feedback(self.plot)
        self.plot.invertX(True)
        self.appearance.apply_grid(self.plot)
        self.plot.setLabel('bottom','¹H [ppm]');self.plot.setLabel('left','Widmo przed korekcją bazy',units='j.u.')
        self.view.phase_dragged.connect(self.drag_phase);self.view.clicked.connect(self.set_pivot)
        self.phase_label=QtWidgets.QLabel();self.phase_label.setWordWrap(True);layout.addWidget(self.phase_label)
        # The raw FID is diagnostic only. Never guess a vendor delay from its maximum.
        corrected_delay=self.initial['filter_delay_samples']
        self.filter_note=QtWidgets.QLabel(
            f'Korekcja filtra: {corrected_delay:.4f} próbki. PH0/PH1 dostrajają już skorygowane widmo. '
            'Opóźnienie można zmienić w Processing. Poniżej: początek surowego FID referencji.'
            if corrected_delay else
            'Początek surowego FID referencji. Jeśli potrzebna jest ręczna korekcja opóźnienia, każda próbka odpowiada rampie 360° / SW.')
        self.filter_note.setWordWrap(True)
        self.filter_note.setStyleSheet('color:#fbbf24;');layout.addWidget(self.filter_note)
        self.fid_plot=pg.PlotWidget();self.fid_plot.setMaximumHeight(120);self.fid_plot.setMinimumHeight(90)
        self.fid_plot.setLabel('bottom','Indeks próbki FID (od 0)');self.fid_plot.setLabel('left','|FID|')
        n=min(100,dataset.raw.fid.shape[1]);ref=self.initial['reference_layer']
        self.fid_plot.plot(np.arange(n),abs(dataset.raw.fid[ref,:n]),pen='#fbbf24',symbol='o',symbolSize=3)
        layout.addWidget(self.fid_plot)
        row=QtWidgets.QHBoxLayout();layout.addLayout(row)
        row.addWidget(QtWidgets.QLabel('Ręcznie: opóźnienie [próbki]'))
        self.delay=QtWidgets.QDoubleSpinBox();self.delay.setRange(-250,250);self.delay.setDecimals(4);self.delay.setValue(0)
        row.addWidget(self.delay)
        button=QtWidgets.QPushButton('Dodaj rampę 360° × opóźnienie');button.clicked.connect(self.add_delay);row.addWidget(button)
        self.delay.setEnabled(not corrected_delay);button.setEnabled(not corrected_delay)
        row.addStretch()
        cancel=QtWidgets.QPushButton('Anuluj');cancel.clicked.connect(self.reject);row.addWidget(cancel)
        self.apply_button=QtWidgets.QPushButton('Zastosuj wspólną fazę do wszystkich warstw')
        self.apply_button.clicked.connect(self.commit);row.addWidget(self.apply_button)
        self.refresh();self.plot.autoRange()

    @staticmethod
    def sensitivity(row,values,default):
        combo=QtWidgets.QComboBox()
        for value in values:combo.addItem(str(value),value)
        combo.setCurrentIndex(values.index(default));row.addWidget(combo);return combo

    def options(self):
        o=deepcopy(self.initial)
        for key,control in self.controls.items():o[key]=control.value()
        return o

    def refresh(self,*args):
        if not hasattr(self,'plot'):return
        o=self.options()
        self.preview=self.unphased*np.exp(1j*phase_radians(self.dataset.ppm,self.dataset.metadata,o))[None,:]
        self.draw()
        self.phase_label.setText(f"PH0={o['phase0_deg']:.3f}° · PH1={o['phase1_deg']:.3f}° / całe SW · pivot={o['pivot_ppm']:.6f} ppm. "
            'Podgląd bez odejmowania bazy. Zmiany dotyczą wszystkich warstw. Ustawienia Delty nie są kopiowane.')

    def draw(self,*args):
        if self.preview is None:return
        self.plot.clear();i=self.layer.currentIndex();x=self.dataset.ppm
        self.plot.plot(x,self.preview[i].real,pen=self.appearance.pen())
        if self.imaginary.isChecked():self.plot.plot(x,self.preview[i].imag,pen=self.appearance.pen('#fb923c',QtCore.Qt.PenStyle.DotLine))
        if self.before.isChecked():self.plot.plot(x,self.dataset.phased[i].real,pen=self.appearance.pen('#64748b',QtCore.Qt.PenStyle.DashLine))
        self.pivot=pg.InfiniteLine(self.controls['pivot_ppm'].value(),movable=not self.pan_button.isChecked(),pen=pg.mkPen('#c084fc',width=2))
        self.pivot.sigPositionChangeFinished.connect(lambda line:self.set_pivot(line.value()))
        self.plot.addItem(self.pivot)

    def set_navigation(self, other_button, checked):
        if checked:
            with QtCore.QSignalBlocker(other_button):
                other_button.setChecked(False)
        self.view.hide_selection_preview()
        panning = self.pan_button.isChecked()
        self.view.mode = 'pan' if panning else 'zoom'
        self.view.phasing = not (panning or self.zoom_button.isChecked())
        self.pivot.setMovable(not panning)

    def set_pivot(self,value):self.controls['pivot_ppm'].setValue(value)

    def pivot_changed(self,value):
        # Moving the pivot is a change of coordinates, not an extra phase edit.
        width=self.dataset.metadata['x_sweep_Hz']/(self.dataset.metadata['x_frequency_Hz']/1e6)
        control=self.controls['phase0_deg']
        ph0=control.value()+self.controls['phase1_deg'].value()*(value-self.previous_pivot)/width
        if abs(ph0)>3600:ph0=(ph0+180)%360-180
        control.blockSignals(True);control.setValue(ph0);control.blockSignals(False)
        self.previous_pivot=value;self.refresh()

    def drag_phase(self, which, dy, start, finish, fine):
        key = 'phase0_deg' if which == 0 else 'phase1_deg'
        sensitivity = self.ph0_sensitivity if which == 0 else self.ph1_sensitivity
        degrees_per_pixel = float(sensitivity.currentData())
        if fine:
            degrees_per_pixel = min(degrees_per_pixel / 20, .25)
        # Scale only the movement since the last event, so Shift can change mid-drag.
        movement = dy if start else dy - self.last_drag_dy
        control = self.controls[key]
        control.setValue(control.value() + movement * degrees_per_pixel)
        self.last_drag_dy = 0. if finish else dy

    def add_delay(self):
        control=self.controls['phase1_deg'];control.setValue(control.value()+360*self.delay.value())

    def reset_phase(self):
        for key in ('phase0_deg','phase1_deg'):self.controls[key].blockSignals(True);self.controls[key].setValue(0);self.controls[key].blockSignals(False)
        self.refresh()

    def toggle_fullscreen(self):
        if self.isFullScreen():self.showNormal()
        else:self.showFullScreen()

    def commit(self):
        self.applied.emit({key:control.value() for key,control in self.controls.items()})
        self.accept()
