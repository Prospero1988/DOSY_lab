"""Shared Qt worker and plot tools for experimental ILT windows."""
import html
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from fullscreen_plot import enable_fullscreen
from i18n import tr
from plot_export import save_light_plot


class SolveSignals(QtCore.QObject):
    finished = QtCore.Signal(int, object, object, str)


class SolveJob(QtCore.QRunnable):
    def __init__(self, revision, parameters, compute):
        super().__init__()
        self.revision, self.parameters = revision, parameters
        self.compute = compute
        self.signals = SolveSignals()

    def run(self):
        result, error = None, ''
        try:
            result = self.compute(self.parameters)
        except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as exc:
            error = str(exc)
        self.signals.finished.emit(self.revision, self.parameters, result, error)


class ILTPlotTools:
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
                formatter = getattr(plot, 'format_coordinate', None)
                text = formatter(x, point.y()) if formatter else f'x = {x:.5g}   y = {point.y():.5g}'
                coordinate.setText(text)
                if plot.fullscreen_controller.window is not None:
                    plot.fullscreen_controller.window.statusBar().showMessage(text)

        plot.scene().sigMouseMoved.connect(read_coordinate)
        return plot

    def export_png(self, plots, name):
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, tr('Zapisz PNG'), f'{name}.png', 'PNG (*.png)')
        if not filename:
            return
        path = Path(filename if filename.lower().endswith('.png') else filename+'.png')
        try:
            if len(plots) == 1:
                save_light_plot(plots[0], path, legend=True, dpi=300, **getattr(plots[0], 'export_options', {}))
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
                save_light_plot(plot, panel, legend=True, **getattr(plot, 'export_options', {}))
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
