"""Spectrum mouse tools. Display scaling never changes spectral data or ROI areas."""
import numpy as np
from PySide6 import QtCore
from shiboken6 import isValid
import pyqtgraph as pg


class SpectrumBox(pg.ViewBox):
    mouse_help = 'Przeciągnięcie: zoom zakresu / aktywne narzędzie · Rolka: wysokość widma · Środkowy przycisk: przesuwanie · F11: tryb okna'
    selected = QtCore.Signal(float, float)
    clicked = QtCore.Signal(float)
    split_requested = QtCore.Signal(str, float)

    def __init__(self):
        super().__init__(enableMenu=False)
        self.mode = 'zoom'
        self.selection_preview = None
        self.selection_labels = []
        self.selection_press = None
        self.selection_plot = None
        self.selection_viewport = None
        self.split_regions = {}
        self.split_drag = None
        self.setMouseMode(self.PanMode)

    def bind_selection_feedback(self, plot):
        """Observe the initial press too, before Qt's drag-distance threshold."""
        self.selection_plot = plot
        self.selection_viewport = plot.viewport()
        self.selection_viewport.installEventFilter(self)

    def selection_kind(self, button):
        if getattr(self, 'phasing', False):
            return None
        if button == QtCore.Qt.MouseButton.RightButton:
            return 'zoom'
        if button == QtCore.Qt.MouseButton.LeftButton and self.mode in ('zoom', 'roi', 'baseline'):
            return self.mode
        return None

    def show_selection_preview(self, start, end, kind):
        if self.selection_preview is None or not isValid(self.selection_preview):
            self.selection_labels = []
            self.selection_preview = pg.LinearRegionItem((start,end), movable=False)
            self.selection_preview.setZValue(1e8)
            self.selection_preview.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
            for i,line in enumerate(self.selection_preview.lines):
                line.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
                label = pg.InfLineLabel(line, text='{value:.3f} ppm', position=.96-i*.08,
                                       fill=(16,24,39,220), color='#e2e8f0')
                label.setAcceptedMouseButtons(QtCore.Qt.MouseButton.NoButton)
                self.selection_labels.append(label)
        item = self.selection_preview
        if item.parentItem() is not self.childGroup:
            self.addItem(item, ignoreBounds=True)
        color = {'roi':(103,232,249), 'baseline':(96,165,250), 'zoom':(250,204,21), 'split':(244,114,182)}[kind]
        item.setBrush(pg.mkBrush(*color,32))
        for line in item.lines:line.setPen(pg.mkPen(color,width=1.5))
        item.setRegion((start,end))
        item.show()
        item.lines[1].setVisible(abs(start-end)>1e-12)
        # InfLineLabel suppresses updates while hidden; refresh after revealing it.
        for label in self.selection_labels:label.valueChanged()

    def hide_selection_preview(self):
        self.selection_press = None
        self.split_drag = None
        if self.selection_preview is not None and isValid(self.selection_preview):self.selection_preview.hide()

    def eventFilter(self, watched, event):
        plot = self.selection_plot
        if watched is self.selection_viewport and plot is not None and isValid(plot):
            kind = event.type()
            if kind == QtCore.QEvent.Type.MouseButtonPress:
                tool = self.selection_kind(event.button())
                pos = plot.mapToScene(event.position().toPoint())
                if tool and self.sceneBoundingRect().contains(pos):
                    # Existing movable ROI/pivot handles keep their normal live editing.
                    editing = event.button()==QtCore.Qt.MouseButton.LeftButton and any(
                        isinstance(item,(pg.LinearRegionItem,pg.InfiniteLine)) and item.movable
                        for item in plot.scene().items(pos))
                    if not editing:
                        x = self.mapSceneToView(pos).x()
                        self.selection_press = (event.button(),x,tool)
                        self.show_selection_preview(x,x,tool)
            elif kind == QtCore.QEvent.Type.MouseMove and self.selection_press is not None:
                button,start,tool = self.selection_press
                if event.buttons() & button:
                    pos = plot.mapToScene(event.position().toPoint())
                    self.show_selection_preview(start,self.mapSceneToView(pos).x(),tool)
            elif kind in (QtCore.QEvent.Type.MouseButtonDblClick,QtCore.QEvent.Type.Hide,QtCore.QEvent.Type.WindowDeactivate):
                self.hide_selection_preview()
            elif kind==QtCore.QEvent.Type.MouseButtonRelease and self.split_drag is None:
                self.hide_selection_preview()
        return super().eventFilter(watched,event)

    def wheelEvent(self, event, axis=None):
        # Change vertical display gain around intensity zero; preserve ppm range.
        factor = float(np.exp(np.clip(-event.delta()/120 * np.log(1.2), -4, 4)))
        self.scaleBy(y=factor, center=QtCore.QPointF(0, 0))
        self.sigRangeChangedManually.emit([False, True])
        event.accept()

    def mouseDragEvent(self, event, axis=None):
        left = event.button() == QtCore.Qt.MouseButton.LeftButton
        right = event.button() == QtCore.Qt.MouseButton.RightButton
        if axis is None and left and self.mode=='split':
            event.accept()
            self.drag_split(event)
        elif axis is None and ((left and self.mode in ('zoom', 'roi', 'baseline')) or right):
            event.accept()
            start, end = event.buttonDownPos(), event.pos()
            a, b = self.mapToView(start).x(), self.mapToView(end).x()
            if not event.isFinish():
                self.show_selection_preview(a,b,'zoom' if right else self.mode)
                return
            self.hide_selection_preview()
            if abs(a-b) < 1e-8:
                return
            if left and self.mode in ('roi', 'baseline'):
                self.selected.emit(a, b)
            else:
                self.setXRange(min(a,b), max(a,b), padding=.02)
                self.sigRangeChangedManually.emit([True, False])
        elif left and self.mode in ('pivot','modify','reference'):
            event.accept()
        else:
            super().mouseDragEvent(event, axis)

    def drag_split(self,event):
        if event.isStart():
            self.split_drag = None
            start = self.mapToView(event.buttonDownPos()).x()
            tolerance = 8 * abs(self.viewPixelSize()[0])
            nearest = tolerance
            for sid,item in self.split_regions.items():
                lo,hi = item.getRegion()
                for edge in (lo,hi):
                    distance = abs(edge-start)
                    if distance<nearest:
                        nearest = distance
                        self.split_drag = (sid,lo,hi,edge)
        if self.split_drag is None:
            return
        sid,lo,hi,edge = self.split_drag
        cut = self.mapToView(event.pos()).x()
        if event.isFinish():
            self.hide_selection_preview()
            if lo<cut<hi:
                self.split_requested.emit(sid,cut)
        else:
            self.show_selection_preview(edge,min(hi,max(lo,cut)),'split')

    def mouseClickEvent(self, event):
        self.hide_selection_preview()
        if self.mode in ('pivot', 'reference') and event.button() == QtCore.Qt.MouseButton.LeftButton:
            event.accept()
            self.clicked.emit(self.mapToView(event.pos()).x())
        else:
            super().mouseClickEvent(event)


class PhaseBox(SpectrumBox):
    phase_dragged = QtCore.Signal(int, float, bool, bool, bool)
    phasing = True

    def mouseDragEvent(self, event, axis=None):
        if self.phasing and axis is None and event.button() in (QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.MouseButton.RightButton):
            event.accept()
            if event.isStart():self.disableAutoRange()
            which = 0 if event.button() == QtCore.Qt.MouseButton.LeftButton else 1
            dy = event.buttonDownScreenPos().y() - event.screenPos().y()
            fine = bool(event.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier)
            self.phase_dragged.emit(which, dy, event.isStart(), event.isFinish(), fine)
        else:
            super().mouseDragEvent(event, axis)

    def mouseClickEvent(self, event):
        if self.mode != 'pan' and event.button() == QtCore.Qt.MouseButton.LeftButton:
            event.accept()
            self.clicked.emit(self.mapToView(event.pos()).x())
        else:
            super().mouseClickEvent(event)
