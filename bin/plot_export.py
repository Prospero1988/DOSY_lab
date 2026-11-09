"""Light figure exports rendered separately from the live plots."""
from pathlib import Path

from PySide6 import QtCore, QtGui
import pyqtgraph as pg
from pyqtgraph.exporters import ImageExporter

from statistics_svg import export_svg


NAVY = '#142b5c'


def light_plot(source):
    plot = pg.PlotWidget(background='w')
    plot.invertX(source.getViewBox().xInverted())
    plot.hideButtons()
    plot.getViewBox().setBorder('k')
    for name in ('left', 'bottom'):
        original = source.getAxis(name)
        axis = plot.getAxis(name)
        axis.setPen('k')
        axis.setTickPen('k')
        axis.setTextPen('k')
        axis.setLabel(original.labelText, units=original.labelUnits, color='#000000')
        axis.setStyle(showValues=original.style['showValues'])
        axis.setLogMode(original.logMode)
        axis.enableAutoSIPrefix(original.autoSIPrefix)
        axis.setGrid(original.grid)
    plot.getAxis('left').setWidth(source.getAxis('left').width())
    title = source.plotItem.titleLabel
    if title.isVisible():
        plot.setTitle(title.text, color='#000000')
        label = plot.plotItem.titleLabel
        height = max(30, label.item.boundingRect().height())
        label.setMaximumHeight(height)
        plot.plotItem.layout.setRowFixedHeight(0, height)
    return plot


def _colored_pen(pen, color):
    result = QtGui.QPen(pg.mkPen(pen))
    result.setColor(QtGui.QColor(color))
    return result


def _light_item(item):
    if isinstance(item, pg.PlotDataItem):
        x, y = item.getData()
        if x is None or len(x) == 0:
            return None
        return pg.PlotDataItem(x, y, pen=_colored_pen(item.opts['pen'], NAVY),
                               symbol=item.opts['symbol'], symbolSize=item.opts['symbolSize'],
                               symbolPen=None, symbolBrush=NAVY, name=item.name())
    if isinstance(item, pg.ScatterPlotItem):
        spots = [
            dict(pos=point.pos(), size=point.size(), symbol=point.symbol(), pen=pg.mkPen(NAVY),
                 brush=pg.mkBrush(None if point.brush().style() == QtCore.Qt.BrushStyle.NoBrush else NAVY))
            for point in item.points()
        ]
        return pg.ScatterPlotItem(spots=spots, useCache=False)
    if isinstance(item, pg.ErrorBarItem):
        return pg.ErrorBarItem(**{**item.opts, 'pen': pg.mkPen('k')})
    if isinstance(item, pg.LinearRegionItem):
        alpha = item.brush.color().alpha()
        return pg.LinearRegionItem(item.getRegion(), orientation=item.orientation, movable=False,
                                   pen=_colored_pen(item.lines[0].pen, '#000000'),
                                   brush=pg.mkBrush(0, 0, 0, alpha))
    if isinstance(item, pg.InfiniteLine):
        return pg.InfiniteLine(item.value(), angle=item.angle, pen=_colored_pen(item.pen, '#000000'))
    if isinstance(item, pg.TextItem):
        label = pg.TextItem(item.toPlainText(), color='k', anchor=item.anchor)
        label.setPos(item.pos())
        return label
    return None


def save_light_plot(source, path, *, bottom_ticks=None, legend=False, dpi=None):
    path = Path(path)
    plot = light_plot(source)
    plot.setAttribute(QtCore.Qt.WidgetAttribute.WA_DontShowOnScreen)
    try:
        if legend:
            plot.addLegend(labelTextColor='k', brush='w')
        if bottom_ticks is not None:
            plot.getAxis('bottom').setTicks([bottom_ticks])
        for item in source.plotItem.items:
            if not item.isVisible():
                continue
            copy = _light_item(item)
            if copy is not None:
                copy.setZValue(item.zValue())
                if legend and isinstance(item, pg.InfiniteLine) and getattr(item, 'label', None) is not None:
                    pg.InfLineLabel(copy, text=item.label.format, position=item.label.orthoPos, color='k')
                plot.addItem(copy)
        plot.resize(max(source.width(), 600), max(source.height(), 300))
        plot.show()
        x_range, y_range = source.getViewBox().viewRange()
        plot.setRange(xRange=x_range, yRange=y_range, padding=0)
        if path.suffix.lower() == '.svg':
            export_svg(plot.plotItem, path)
        elif path.suffix.lower() == '.png':
            exporter = ImageExporter(plot.plotItem)
            exporter.parameters()['width'] = 1600
            if dpi is None:
                exporter.export(str(path))
            else:
                image = exporter.export(toBytes=True)
                image.setDotsPerMeterX(round(dpi / .0254))
                image.setDotsPerMeterY(round(dpi / .0254))
                if not image.save(str(path), 'PNG'):
                    raise OSError(f'Could not save PNG: {path}')
        else:
            raise ValueError('Wybierz rozszerzenie PNG lub SVG.')
    finally:
        plot.close()
        plot.deleteLater()
