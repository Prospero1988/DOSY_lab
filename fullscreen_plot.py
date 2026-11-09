"""Temporarily move a live plot into a fullscreen window; never copy/refit data."""
from copy import deepcopy
import weakref

from PySide6 import QtCore, QtGui, QtWidgets
from shiboken6 import isValid


_controllers = weakref.WeakSet()


class FullscreenPlotWindow(QtWidgets.QMainWindow):
    def __init__(self, controller):
        super().__init__(controller.owner, QtCore.Qt.WindowType.Window)
        self.controller = controller
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)
        self.setWindowTitle(controller.title+' · Pełny ekran')
        self.setStyleSheet(controller.owner.styleSheet())
        central = QtWidgets.QWidget(); self.setCentralWidget(central)
        self.box = QtWidgets.QVBoxLayout(central)
        row = QtWidgets.QHBoxLayout(); self.box.addLayout(row)
        back = QtWidgets.QPushButton('Wróć do małego widoku (Esc)'); back.clicked.connect(self.close); row.addWidget(back)
        all_data = QtWidgets.QPushButton('Cały zakres'); all_data.clicked.connect(controller.reset_view or controller.plot.autoRange); row.addWidget(all_data)
        hint=getattr(controller.plot.getViewBox(),'mouse_help','Zoom: kółko myszy · Przesuwanie: przeciągnięcie · F11: tryb okna')
        row.addStretch(); row.addWidget(QtWidgets.QLabel(hint))
        self.context = QtWidgets.QLabel(); self.context.setWordWrap(True)
        self.context.setTextFormat(QtCore.Qt.TextFormat.RichText); self.box.addWidget(self.context)
        self.box.addWidget(controller.plot, 1)
        self.shortcuts = []
        for key, callback in [('Esc', self.close), ('F11', self.toggle_fullscreen)]:
            shortcut = QtGui.QShortcut(QtGui.QKeySequence(key), self)
            shortcut.setContext(QtCore.Qt.ShortcutContext.WindowShortcut); shortcut.activated.connect(callback)
            self.shortcuts.append(shortcut)
        self.context_timer = QtCore.QTimer(self); self.context_timer.timeout.connect(self.update_context)
        self.context_timer.start(200); self.update_context()

    def update_context(self):
        text = self.controller.context_provider() if self.controller.context_provider else ''
        if text != self.context.text(): self.context.setText(text)
        self.context.setVisible(bool(text))

    def toggle_fullscreen(self):
        if self.isFullScreen(): self.showMaximized()
        else: self.showFullScreen()

    def closeEvent(self, event):
        self.context_timer.stop(); self.controller.restore()
        super().closeEvent(event)


class FullscreenPlotController(QtCore.QObject):
    def __init__(self, plot, title, context_provider=None, *, reset_view=None):
        super().__init__(plot)
        self.plot = plot; self.title = title; self.context_provider = context_provider
        self.reset_view = reset_view
        self.owner = plot.window(); self.window = None; self.placeholder = None
        self.viewport = plot.viewport()
        self.original_parent = None; self.original_layout = None; self.splitter_sizes = None
        self.saved_views = []
        plot.viewport().installEventFilter(self); self.owner.installEventFilter(self)
        if reset_view is None:
            plot.setToolTip('Dwuklik lewym przyciskiem otwiera pełny ekran. Esc lub ponowny dwuklik przywraca okno.')
        else:
            plot.setToolTip('Przeciągnij, aby powiększyć obszar. Dwuklik przywraca cały zakres.')
        _controllers.add(self)

    def eventFilter(self, watched, event):
        if not isValid(self.plot): return False
        if watched is self.viewport and event.type() == QtCore.QEvent.Type.MouseButtonDblClick:
            if event.button() == QtCore.Qt.MouseButton.LeftButton:
                if self.reset_view is not None: self.reset_view()
                elif self.window is None: self.open()
                else: self.window.close()
                event.accept(); return True
        if watched is self.owner and event.type() in (QtCore.QEvent.Type.Close, QtCore.QEvent.Type.Hide, QtCore.QEvent.Type.DeferredDelete):
            if self.window is not None: self.window.close()
        return super().eventFilter(watched, event)

    def suspend_links(self):
        # ViewBox screen-geometry alignment must not connect windows of
        # different sizes. Keep the declarations, suspend callbacks and restore
        # the original small-view zoom/auto-range modes on return.
        views = {self.plot.getViewBox()}
        available = [c.plot.getViewBox() for c in list(_controllers) if isValid(c) and isValid(c.plot)]
        changed = True
        while changed:
            changed = False
            for view in available:
                targets = [view.linkedView(axis) for axis in (0, 1)]
                if view in views or any(target in views for target in targets):
                    for linked in [view]+[t for t in targets if t is not None]:
                        if linked not in views: views.add(linked); changed = True
        self.saved_views = [(v, deepcopy(v.viewRange()), list(v.autoRangeEnabled()), v.linksBlocked) for v in views]
        for view in views: view.blockLink(True)

    def open(self):
        if self.window is not None:
            self.window.showFullScreen(); self.window.raise_(); self.window.activateWindow(); return
        current = getattr(self.owner, '_fullscreen_plot_window', None)
        if current is not None: current.close()
        self.original_parent = self.plot.parentWidget()
        placeholder = QtWidgets.QFrame(self.original_parent); box = QtWidgets.QVBoxLayout(placeholder)
        label = QtWidgets.QLabel(self.title+' — otwarty w osobnym oknie.'); label.setWordWrap(True); box.addWidget(label)
        show = QtWidgets.QPushButton('Pokaż pełny ekran'); show.clicked.connect(self.open); box.addWidget(show)
        back = QtWidgets.QPushButton('Przywróć wykres tutaj'); back.clicked.connect(self.close); box.addWidget(back)
        if isinstance(self.original_parent, QtWidgets.QSplitter):
            self.splitter_sizes = self.original_parent.sizes()
            self.original_parent.replaceWidget(self.original_parent.indexOf(self.plot), placeholder)
        else:
            self.original_layout = self.original_parent.layout()
            if self.original_layout is None or self.original_layout.replaceWidget(self.plot, placeholder) is None:
                placeholder.deleteLater(); raise RuntimeError('Nie można odnaleźć miejsca wykresu w interfejsie.')
        self.placeholder = placeholder
        self.suspend_links()
        self.window = FullscreenPlotWindow(self)
        self.owner._fullscreen_plot_window = self.window
        screen = self.owner.screen()
        if screen is not None: self.window.setGeometry(screen.availableGeometry())
        self.window.showFullScreen(); self.window.raise_(); self.window.activateWindow()
        self.plot.setFocus()

    def close(self):
        if self.window is not None: self.window.close()

    def restore(self):
        window = self.window
        if window is None: return
        self.window = None
        if getattr(self.owner, '_fullscreen_plot_window', None) is window: self.owner._fullscreen_plot_window = None
        window.box.removeWidget(self.plot)
        if isValid(self.original_parent) and isValid(self.placeholder):
            if isinstance(self.original_parent, QtWidgets.QSplitter):
                self.original_parent.replaceWidget(self.original_parent.indexOf(self.placeholder), self.plot)
                self.original_parent.setSizes(self.splitter_sizes)
            else:
                self.original_layout.replaceWidget(self.placeholder, self.plot)
                self.plot.setParent(self.original_parent)
                self.original_layout.activate()
            self.plot.show(); self.placeholder.hide(); self.placeholder.deleteLater()
        for view, ranges, automatic, blocked in self.saved_views:
            if not isValid(view): continue
            view.setRange(xRange=ranges[0], yRange=ranges[1], padding=0, disableAutoRange=False)
            for axis in (0, 1): view.enableAutoRange(axis, automatic[axis])
            view.blockLink(blocked)
        self.saved_views = []; self.placeholder = None


def enable_fullscreen(plot, title, context_provider=None, *, reset_view=None):
    if getattr(plot, 'fullscreen_controller', None) is None:
        plot.fullscreen_controller = FullscreenPlotController(plot, title, context_provider, reset_view=reset_view)
    return plot.fullscreen_controller
