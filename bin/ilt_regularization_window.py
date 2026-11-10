"""Shared regularization diagnostics and explicit application for the ILT views."""
from threading import Event

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from i18n import tr
from ilt import D_UNIT
from ilt_help import ILT_WINDOW_FLAGS
from ilt_regularization import SweepOptions, scan
from ilt_tools import ILTPlotTools, SolveJob


class RegularizationControls:
    def __init__(self, parent, layout, snapshot, token, apply, restore, on_manual=None):
        self.parent = parent
        self.snapshot, self.token = snapshot, token
        self.apply, self.restore, self.on_manual = apply, restore, on_manual
        self.previous = None
        self.dialog = None
        self.mode = QtWidgets.QComboBox()
        self.mode.addItems([tr('Ręczna λ'),tr('Automatyczna λ')])
        self.mode.currentIndexChanged.connect(self.change_mode)
        layout.addWidget(self.mode)
        button = QtWidgets.QPushButton(tr('Auto λ / stabilność'))
        button.clicked.connect(self.open)
        layout.addWidget(button)
        self.label = QtWidgets.QLabel()
        layout.addWidget(self.label)

    def set_mode(self,index):
        with QtCore.QSignalBlocker(self.mode):
            self.mode.setCurrentIndex(index)

    def change_mode(self,index):
        if index:
            self.set_mode(int(bool(self.label.text())))
            self.open()
        else:
            if self.on_manual is not None:
                self.on_manual()
            self.manual_override()

    def manual_override(self):
        self.set_mode(0)
        self.label.clear()

    def open(self):
        if self.dialog is not None:
            self.dialog.raise_()
            self.dialog.activateWindow()
            return
        try:
            token, current, factory, context = self.snapshot()
        except ValueError as exc:
            QtWidgets.QMessageBox.information(self.parent,tr('Auto λ / stabilność'),tr(str(exc)))
            return
        self.dialog = RegularizationWindow(self.parent,self,token,current,factory,context)
        self.dialog.finished.connect(self.closed)
        self.dialog.show()

    def dismiss(self):
        if self.dialog is not None:
            self.dialog.close()

    def closed(self):
        self.dialog = None

    def accept(self,token,recommendations):
        if self.token() != token:
            return False
        label = self.label.text()
        self.previous = self.apply(recommendations),label
        self.set_mode(1)
        methods = {tr(result.method) for result in recommendations.values()}
        summary = next(iter(methods)) if len(methods)==1 else tr('Różne metody')
        values = [result.recommendation for result in recommendations.values()]
        self.label.setText((f'λ = {values[0]:.3g}' if len(values)==1 else f'λ · {len(values)} {tr("grupy")}')
                           + ' · ' + summary)
        return True

    def restore_previous(self):
        if self.previous is not None:
            value,label = self.previous
            self.previous = None
            self.restore(value)
            self.label.setText(label)
            self.set_mode(int(bool(label)))


class RegularizationWindow(QtWidgets.QDialog, ILTPlotTools):
    def __init__(self,parent,controls,token,current,factory,context):
        super().__init__(parent,ILT_WINDOW_FLAGS)
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)
        self.setWindowTitle(tr('Stabilność regularyzacji ILT'))
        self.resize(1200,760)
        self.controls,self.token,self.current,self.factory = controls,token,current,factory
        self.context = context
        self.job,self.cancel = None,None
        self.results,self.sources,self.cache = {},{},{}
        self.preview_index = None
        self.export_buttons = []
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(QtWidgets.QLabel(context))
        bar = QtWidgets.QHBoxLayout()
        layout.addLayout(bar)
        self.minimum = QtWidgets.QLineEdit(f'{max(current/1000,1e-12):.6g}')
        self.maximum = QtWidgets.QLineEdit(f'{min(current*1000,100):.6g}')
        self.count = QtWidgets.QSpinBox()
        self.count.setRange(7,41)
        self.count.setValue(13)
        self.method = QtWidgets.QComboBox()
        for label,value in (('Automatyczny wybór kryterium','auto'),('Krzywa L / warunkowa','corner'),
                            ('Zasada rozbieżności','noise')):
            self.method.addItem(tr(label),value)
        for label,widget in (('λ min',self.minimum),('λ max',self.maximum),('Liczba prób',self.count),
                             ('Kryterium',self.method)):
            bar.addWidget(QtWidgets.QLabel(tr(label)))
            bar.addWidget(widget)
        self.noise = QtWidgets.QCheckBox(tr('Potwierdzam wiarygodne σ(I) z pomiaru'))
        self.noise.setToolTip(tr('Użyj wyłącznie niezależnie oszacowanych błędów standardowych całek.'))
        self.grid = QtWidgets.QCheckBox(tr('Sprawdź siatkę D'))
        bar.addWidget(self.noise)
        bar.addWidget(self.grid)
        self.run = QtWidgets.QPushButton(tr('Przelicz skan λ'))
        self.run.clicked.connect(self.start)
        bar.addWidget(self.run)
        self.stop = QtWidgets.QPushButton(tr('Anuluj skan'))
        self.stop.clicked.connect(self.cancel_scan)
        self.stop.setEnabled(False)
        bar.addWidget(self.stop)
        self.inputs = (self.minimum,self.maximum,self.count,self.method,self.noise,self.grid)
        self.target = QtWidgets.QComboBox()
        self.target.currentIndexChanged.connect(self.render)
        layout.addWidget(self.target)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        tabs = QtWidgets.QTabWidget()
        layout.addWidget(tabs,1)
        plots = []
        for tab,title,x,y in (
            ('Krzywa','A · Krzywa L / warunkowa','Norma reszty funkcji celu','Norma kary λ'),
            ('Reszty','B · Reszty a λ','λ','Znormalizowana reszta'),
            ('Kary','C · Kary a λ','λ','Norma kary λ'),
            ('Maksima','D · Maksima D a λ','λ','D [m²/s]'),
            ('Podgląd','E · Podgląd rozkładu','D [m²/s]','Względna masa komórki'),
        ):
            page = QtWidgets.QWidget()
            tabs.addTab(page,tr(tab))
            plots.append(self.create_plot(QtWidgets.QVBoxLayout(page),title,tr(x),tr(y)))
        self.plots = tuple(plots)
        for plot in self.plots[1:4]:
            plot.setLogMode(x=True,y=False)
        self.plots[3].setLogMode(x=True,y=True)
        self.plots[4].setLogMode(x=True,y=False)
        self.details = QtWidgets.QPlainTextEdit(readOnly=True)
        self.details.setMaximumHeight(155)
        layout.addWidget(self.details)
        actions = QtWidgets.QHBoxLayout()
        layout.addLayout(actions)
        self.apply_button = QtWidgets.QPushButton(tr('Zastosuj rekomendację'))
        self.apply_button.clicked.connect(self.apply)
        self.apply_button.setEnabled(False)
        self.restore_button = QtWidgets.QPushButton(tr('Przywróć poprzednie λ'))
        self.restore_button.setEnabled(controls.previous is not None)
        self.restore_button.clicked.connect(self.restore)
        export = QtWidgets.QPushButton(tr('Zapisz wszystkie wykresy PNG'))
        export.clicked.connect(lambda: self.export_png(self.plots,'ILT_lambda'))
        actions.addWidget(self.apply_button)
        actions.addWidget(self.restore_button)
        actions.addWidget(export)
        self.export_buttons.append(export)
        for signal in (self.minimum.textChanged,self.maximum.textChanged,self.count.valueChanged,
                       self.method.currentIndexChanged,self.noise.toggled,self.grid.toggled):
            signal.connect(self.settings_changed)
        self.start()

    def settings_changed(self,*_):
        self.apply_button.setEnabled(False)
        if self.results:
            self.status.setText(tr('Ustawienia skanu zmienione. Przelicz skan, aby zastosować wynik.'))

    def cancel_scan(self):
        if self.cancel is not None:
            self.cancel.set()
            self.stop.setEnabled(False)

    def start(self):
        if self.job is not None:
            return
        try:
            options = SweepOptions(*(float(w.text().strip().replace(',','.')) for w in (self.minimum,self.maximum)),
                                   self.count.value(),self.method.currentData(),self.noise.isChecked(),self.grid.isChecked())
        except ValueError:
            self.status.setText(tr('Wpisz prawidłową wartość liczbową.'))
            return
        if options in self.cache:
            self.show_results(*self.cache[options])
            return
        self.cancel = cancel = Event()
        self.run.setEnabled(False)
        self.stop.setEnabled(True)
        for widget in self.inputs:
            widget.setEnabled(False)
        self.apply_button.setEnabled(False)
        self.status.setText(tr('Obliczanie skanu λ…'))
        for button in self.export_buttons:
            button.setEnabled(False)

        def compute(settings):
            sources = self.factory()
            return sources,{key:scan(source,settings,cancel.is_set) for key,source in sources.items()}

        self.job = SolveJob(0,options,compute)
        self.job.signals.finished.connect(self.solved)
        QtCore.QThreadPool.globalInstance().start(self.job)

    @QtCore.Slot(int,object,object,str)
    def solved(self,_,options,payload,error):
        self.job,self.cancel = None,None
        self.run.setEnabled(True)
        self.stop.setEnabled(False)
        for widget in self.inputs:
            widget.setEnabled(True)
        if error:
            self.status.setText(tr(error))
            return
        if len(self.cache)>=4:
            self.cache.pop(next(iter(self.cache)))
        self.cache[options] = payload
        self.show_results(*payload)

    def show_results(self,sources,results):
        self.sources,self.results = sources,results
        with QtCore.QSignalBlocker(self.target):
            self.target.clear()
            for key in results:
                self.target.addItem(', '.join(key),key)
        self.apply_button.setEnabled(any(r.recommendation is not None for r in results.values()))
        for button in self.export_buttons:
            button.setEnabled(bool(results))
        if not results:
            self.status.setText(tr('Brak zgodnych grup ROI do skanu λ.'))
        self.render()

    def scatter(self,plot,x,y,indices=None):
        item = pg.PlotDataItem(x,y,pen=pg.mkPen('#67e8f9',width=1.2),
            symbol='o',symbolSize=7,symbolBrush='#0891b2',symbolPen=None,data=indices)
        item.sigPointsClicked.connect(lambda _,points,*event: self.preview(int(points[0].data())))
        plot.addItem(item)

    def render(self,*_):
        result = self.results.get(self.target.currentData())
        if result is None:
            return
        c = result.candidates
        valid = [i for i,item in enumerate(c) if item.valid]
        lams = np.array([item.regularization for item in c])
        self.plots[0].setLogMode(x=True,y=True)
        for plot in self.plots:
            plot.clear()
        if valid:
            self.scatter(self.plots[0],[c[i].objective_residual for i in valid],[c[i].penalty for i in valid],valid)
            self.scatter(self.plots[1],lams[valid],[c[i].residual for i in valid],valid)
            self.plots[1].plot(lams,[i.rms for i in c],pen=pg.mkPen('#fb923c',width=1.2))
            self.scatter(self.plots[2],lams[valid],[c[i].penalty for i in valid],valid)
        for _,track in result.tracks:
            self.scatter(self.plots[3],[lams[p[0]] for p in track],[p[1] for p in track],[p[0] for p in track])
        if result.expected:
            for level in (result.expected,*result.noise_band):
                self.plots[1].addItem(pg.InfiniteLine(level,angle=0,pen=pg.mkPen('#a7f3d0',style=QtCore.Qt.PenStyle.DotLine)))
        current = self.sources[self.target.currentData()].current
        for lam,color in ((current,'#dbe4ef'),(result.recommendation,'#4ade80')):
            if valid and lam is not None and lams[0]<=lam<=lams[-1]:
                i = min(valid,key=lambda i:abs(np.log(lams[i]/lam)))
                if c[i].objective_residual>0 and c[i].penalty>0:
                    self.plots[0].plot([c[i].objective_residual],[c[i].penalty],pen=None,
                        symbol='o',symbolSize=13,symbolPen=pg.mkPen(color,width=2),symbolBrush=None)
        for plot in self.plots[1:4]:
            for lam,color in ((current,'#dbe4ef'),(result.recommendation,'#4ade80')):
                if lam is not None and lams[0]<=lam<=lams[-1]:
                    plot.addItem(pg.InfiniteLine(np.log10(lam),angle=90,pen=pg.mkPen(color,style=QtCore.Qt.PenStyle.DashLine)))
        self.selection = []
        for plot in self.plots[1:4]:
            line = pg.InfiniteLine(angle=90,pen=pg.mkPen('#facc15',style=QtCore.Qt.PenStyle.DotLine))
            plot.addItem(line,ignoreBounds=True)
            self.selection.append(line)
        rec = f'{result.recommendation:.5g}' if result.recommendation else '—'
        band = f'{result.interval[0]:.3g}–{result.interval[1]:.3g}' if result.interval else '—'
        joint = f' · μ = {result.fixed_mu:.4g}' if result.fixed_mu is not None else ''
        self.status.setText(f'{tr(result.method)} · λ = {rec}{joint} · {tr("Zakres kandydatów")} {band}'
            f' · {tr("Próby")} {lams[0]:.3g}–{lams[-1]:.3g} ({len(valid)}/{len(c)}) · {tr("Bieżące λ")} {current:.4g}\n'
            + ' '.join(map(tr,result.warnings)))
        preferred = result.recommendation or current
        self.preview(min(valid,key=lambda i:abs(np.log(lams[i]/preferred))) if valid else 0)

    def preview(self,index):
        result = self.results.get(self.target.currentData())
        if result is None:
            return
        c = result.candidates[index]
        self.preview_index = index
        for line in self.selection:
            line.setValue(np.log10(c.regularization))
        plot = self.plots[4]
        plot.clear()
        plot.setLabel('left',tr('Względna masa komórki' if result.fixed_mu is None
                                else 'Znormalizowana norma L₂ wiersza'))
        if c.valid:
            solution = c.solution
            weights = np.asarray(solution.weights)
            weights = weights if weights.ndim==1 else np.linalg.norm(weights,axis=1)
            if weights.sum()>0:
                plot.plot(solution.diffusion,weights/weights.sum(),pen=pg.mkPen('#67e8f9',width=2))
        details = [f'λ = {c.regularization:.5g} · {tr(c.status)} · RMS = {c.rms:.4g} · '
                   f'{tr("Ważona norma reszty")} {c.residual:.4g} · ‖Lf‖ = {c.penalty:.4g}']
        if result.fixed_mu is not None:
            details.append(f'μ √M ‖F‖₂,₁ = {c.sparsity_cost:.4g} · {tr("Krzywa używa reszty z kosztem μ.")}')
        if result.corner is not None:
            details.append(f"{tr('Załamanie krzywej')}: λ = {result.corner:.4g}")
        if result.expected is not None:
            details.append(f'{tr("Oczekiwana norma szumu")} {result.expected:.4g} '
                           f'[{result.noise_band[0]:.4g}–{result.noise_band[1]:.4g}]')
        if c.solution is not None and hasattr(c.solution,'stationarity'):
            details.append(f"{tr('Iteracje')}: {c.solution.iterations} · {tr('Reszta stacjonarności')}: "
                           f'{c.solution.stationarity:.4g}')
        details.append(f"{tr('Liczba maksimów')}: {len(c.maxima)}")
        for status,track in result.tracks:
            for position,d,fraction,_ in track:
                if position==index:
                    details.append(f'{d:.4e} m²/s · {fraction:.2%} · {tr(status)}')
                    masses = c.masses[c.maxima.index(d)]
                    ids = self.target.currentData()
                    details.append('I(b=0) [j.u.·ppm] · ' +
                                   ', '.join(f'{sid}={mass:.4g}' for sid,mass in zip(ids,masses)))
        failed = [f'{point.regularization:.3g} ({tr(point.status)})' for point in result.candidates if not point.valid]
        if failed:
            details.append(f"{tr('Nieudane próby λ')}: " + ', '.join(failed))
        for label,record in result.grid_checks:
            d = ', '.join(f'{d:.4e}' for d in record.maxima) or '—'
            details.append(f'{tr(label)} · λ = {record.regularization:.4g} · {d} m²/s · {tr(record.status)}')
        details.append(tr('Biała kreska: bieżące λ. Zielona: sugestia. Żółta: podgląd. Pomarańczowa: RMS.'))
        self.details.setPlainText('\n'.join(details))

    def apply(self):
        recommendations = {key:value for key,value in self.results.items() if value.recommendation is not None}
        if not recommendations or not self.apply_button.isEnabled():
            return
        if self.controls.accept(self.token,recommendations):
            self.accept()
        else:
            self.status.setText(tr('Dane lub ustawienia się zmieniły. Otwórz nowy skan λ.'))
            self.apply_button.setEnabled(False)

    def restore(self):
        self.controls.restore_previous()
        self.reject()

    def done(self,result):
        if self.cancel is not None:
            self.cancel.set()
        super().done(result)
