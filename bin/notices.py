"""Consistent visual treatment of actionable messages and QC states."""
from PySide6 import QtCore, QtGui, QtWidgets


COLORS = {
    'error': ('#fecaca', '#381d28', '#a34552'),
    'warning': ('#fde7b0', '#332719', '#94703e'),
}
QC_LEVEL = {'NIEWAŻNY': 'error', 'UWAGA': 'warning'}


def notice_icon(level):
    icon = (QtWidgets.QStyle.StandardPixmap.SP_MessageBoxCritical if level == 'error'
            else QtWidgets.QStyle.StandardPixmap.SP_MessageBoxWarning)
    return QtWidgets.QApplication.style().standardIcon(icon)


class NoticeLabel(QtWidgets.QWidget):
    def __init__(self, text=''):
        super().__init__()
        self.setObjectName('dosyNotice')
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground)
        self.icon = QtWidgets.QLabel()
        self.message = QtWidgets.QLabel()
        self.message.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        self.message.setWordWrap(True)
        row = QtWidgets.QHBoxLayout(self)
        row.setSpacing(6)
        row.addWidget(self.icon, alignment=QtCore.Qt.AlignmentFlag.AlignTop)
        row.addWidget(self.message, 1)
        self.timer = QtCore.QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(lambda: self.set_notice(''))
        self.set_notice(text)

    def set_notice(self, text, level=None, timeout=0):
        self.timer.stop()
        level = level if text else None
        self.setProperty('noticeLevel', level)
        if level:
            color, background, border = COLORS[level]
            self.setStyleSheet(f'QWidget#dosyNotice {{ background:{background}; border:1px solid {border}; border-radius:4px; }}')
            self.message.setStyleSheet(f'color:{color}; background:transparent; font-weight:bold;')
            self.icon.setStyleSheet('background:transparent;')
            self.icon.setPixmap(notice_icon(level).pixmap(16, 16))
            self.layout().setContentsMargins(8, 5, 8, 5)
        else:
            self.setStyleSheet('')
            self.message.setStyleSheet('')
            self.layout().setContentsMargins(0, 0, 0, 0)
        self.icon.setVisible(bool(level))
        self.message.setText(text)
        self.setVisible(bool(text))
        if timeout and text:
            self.timer.start(timeout)

    def text(self):
        return self.message.text()


def style_qc_item(item, qc):
    level = QC_LEVEL.get(qc)
    if level is None:
        return
    item.setIcon(notice_icon(level))
    item.setForeground(QtGui.QColor(COLORS[level][0]))
    font = item.font()
    font.setBold(True)
    item.setFont(font)
