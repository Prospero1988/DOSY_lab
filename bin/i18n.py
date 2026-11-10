"""Presentation translations. Analysis data and project formats are language independent."""
import json
from pathlib import Path
import re

from PySide6 import QtCore, QtGui, QtWidgets
import pyqtgraph as pg


_english = json.loads(Path(__file__).with_name('english.json').read_text(encoding='utf-8'))
_catalogs = {'en': _english, 'pl': {value: key for key, value in _english.items()}}
_language = 'pl'
_qt_translator = None


def _phrase_pattern(catalog):
    # Longer phrases also appear inside formatted summaries and error messages.
    phrases = sorted((key for key in catalog if len(key) >= 12 or key != key.strip()), key=len, reverse=True)
    return re.compile('|'.join(
        (r'(?<!\w)' if key[0].isalnum() else '') + re.escape(key)
        + (r'(?!\w)' if key[-1].isalnum() else '') for key in phrases
    ))


_patterns = {code: _phrase_pattern(catalog) for code, catalog in _catalogs.items()}


def _convert(text, language):
    if not text:
        return text
    catalog = _catalogs[language]
    if text in catalog:
        return catalog[text]
    return _patterns[language].sub(lambda match: catalog[match[0]], text)


def tr(text):
    """Render a Polish source string in the active interface language."""
    return _convert(text, 'en') if _language == 'en' else text


def current_language():
    return _language


def set_language(code):
    global _language, _qt_translator
    _language = 'en' if code == 'en' else 'pl'
    app = QtCore.QCoreApplication.instance()
    if app is None:
        return
    if _qt_translator is not None:
        app.removeTranslator(_qt_translator)
        _qt_translator.deleteLater()
        _qt_translator = None
    if _language == 'pl':
        translator = QtCore.QTranslator(app)
        directory = QtCore.QLibraryInfo.path(QtCore.QLibraryInfo.LibraryPath.TranslationsPath)
        if translator.load('qtbase_pl', directory):
            app.installTranslator(translator)
            _qt_translator = translator


def retranslate(root, code):
    """Translate widgets in place, preserving analysis state and plot ranges."""
    def update(getter, setter):
        old = getter()
        new = _convert(old, code)
        if new != old:
            setter(new)

    for widget in [root, *root.findChildren(QtWidgets.QWidget)]:
        update(widget.toolTip, widget.setToolTip)
        if widget.isWindow() or isinstance(widget, QtWidgets.QDockWidget):
            update(widget.windowTitle, widget.setWindowTitle)
        if isinstance(widget, (QtWidgets.QAbstractButton, QtWidgets.QLabel)):
            update(widget.text, widget.setText)
        elif isinstance(widget, QtWidgets.QTabWidget):
            for i in range(widget.count()):
                update(lambda: widget.tabText(i), lambda value: widget.setTabText(i, value))
        elif isinstance(widget, QtWidgets.QComboBox):
            with QtCore.QSignalBlocker(widget):
                for i in range(widget.count()):
                    if widget.itemText(i) != widget.itemData(i):
                        update(lambda: widget.itemText(i), lambda value: widget.setItemText(i, value))
        elif isinstance(widget, QtWidgets.QLineEdit):
            update(widget.placeholderText, widget.setPlaceholderText)
        elif isinstance(widget, QtWidgets.QStatusBar):
            update(widget.currentMessage, widget.showMessage)
        elif isinstance(widget, QtWidgets.QTableWidget):
            with QtCore.QSignalBlocker(widget):
                for i in range(widget.columnCount()):
                    header = widget.horizontalHeaderItem(i)
                    if header is not None:
                        update(header.text, header.setText)
                        update(header.toolTip, header.setToolTip)
                    for row in range(widget.rowCount()):
                        item = widget.item(row, i)
                        if item is not None:
                            if header is not None and header.text() == 'QC':
                                update(item.text, item.setText)
                            update(item.toolTip, item.setToolTip)
        elif isinstance(widget, QtWidgets.QListWidget):
            with QtCore.QSignalBlocker(widget):
                for i in range(widget.count()):
                    item = widget.item(i)
                    update(item.text, item.setText)
                    update(item.toolTip, item.setToolTip)
        elif isinstance(widget, pg.PlotWidget):
            for name in ('bottom', 'left', 'top', 'right'):
                axis = widget.getAxis(name)
                label = _convert(axis.labelText, code)
                units = _convert(axis.labelUnits, code)
                if (label, units) != (axis.labelText, axis.labelUnits):
                    axis.setLabel(label, units=units, **axis.labelStyle)
            title = widget.plotItem.titleLabel
            if title.isVisible():
                update(lambda: title.text, widget.setTitle)
    for action in root.findChildren(QtGui.QAction):
        update(action.text, action.setText)
        update(action.toolTip, action.setToolTip)
