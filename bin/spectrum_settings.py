"""Per-user spectrum appearance and its editor."""
from i18n import tr

from dataclasses import dataclass

from PySide6 import QtCore, QtGui, QtWidgets
import pyqtgraph as pg


@dataclass(frozen=True)
class SpectrumAppearance:
    color: str = '#67e8f9'
    width: float = 1.3
    grid: bool = True

    @classmethod
    def load(cls, settings):
        default = cls()
        color = QtGui.QColor(settings.value('spectrum/color', default.color, type=str))
        width = settings.value('spectrum/width', default.width, type=float)
        grid = settings.value('spectrum/grid', default.grid, type=bool)
        return cls(color.name() if color.isValid() else default.color,
                   width if .5 <= width <= 6 else default.width, grid)

    def save(self, settings):
        settings.setValue('spectrum/color', self.color)
        settings.setValue('spectrum/width', self.width)
        settings.setValue('spectrum/grid', self.grid)
        settings.sync()

    def pen(self, color=None, style=QtCore.Qt.PenStyle.SolidLine):
        return pg.mkPen(color or self.color, width=self.width, style=style)

    def apply_grid(self, plot):
        plot.showGrid(x=self.grid, y=self.grid, alpha=.15)


class SpectrumSettingsDialog(QtWidgets.QDialog):
    def __init__(self, appearance, parent=None):
        super().__init__(parent)
        self.setWindowTitle(tr('Ustawienia widma'))
        self.color = appearance.color
        form = QtWidgets.QFormLayout(self)
        self.color_button = QtWidgets.QPushButton()
        self.color_button.clicked.connect(self.choose_color)
        self.show_color()
        form.addRow(tr('Kolor głównej linii'), self.color_button)
        self.width = QtWidgets.QDoubleSpinBox()
        self.width.setRange(.5, 6)
        self.width.setDecimals(1)
        self.width.setSingleStep(.1)
        self.width.setSuffix(' px')
        self.width.setValue(appearance.width)
        form.addRow(tr('Grubość linii'), self.width)
        self.grid = QtWidgets.QCheckBox(tr('Pokaż siatkę'))
        self.grid.setChecked(appearance.grid)
        form.addRow(self.grid)
        buttons = QtWidgets.QDialogButtonBox()
        buttons.addButton(tr('Zapisz'), QtWidgets.QDialogButtonBox.ButtonRole.AcceptRole)
        buttons.addButton(tr('Anuluj'), QtWidgets.QDialogButtonBox.ButtonRole.RejectRole)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

    def show_color(self):
        swatch = QtGui.QPixmap(22, 14)
        swatch.fill(QtGui.QColor(self.color))
        self.color_button.setIcon(QtGui.QIcon(swatch))
        self.color_button.setText(self.color)

    def choose_color(self):
        selected = QtWidgets.QColorDialog.getColor(QtGui.QColor(self.color), self, tr('Kolor widma'))
        if selected.isValid():
            self.color = selected.name()
            self.show_color()

    def appearance(self):
        return SpectrumAppearance(self.color, self.width.value(), self.grid.isChecked())
