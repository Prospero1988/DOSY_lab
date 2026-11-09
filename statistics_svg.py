"""Vector export compatibility for Qt 6.12 close-path tokens in pyqtgraph 0.14."""
import importlib


def export_svg(plot_item, path):
    """Use pyqtgraph's item transforms, retaining SVG closed subpaths.

    Qt 6.12 emits an explicit Z token. pyqtgraph 0.14 assumes every path
    token contains coordinates. Temporarily provide a coordinate-bearing
    closing segment to its transformer, then restore Z at the same position.
    Export is synchronous on the GUI thread; the library is restored on error.
    """
    module = importlib.import_module('pyqtgraph.exporters.SVGExporter')
    original = module.correctCoordinates

    def compatible_coordinates(node, *args):
        closures = []
        for element in node.getElementsByTagName('path'):
            tokens = element.getAttribute('d').strip().split()
            start = None
            positions = []
            for i, token in enumerate(tokens):
                if token.startswith('M'):
                    start = token[1:]
                elif token == 'Z' and start is not None:
                    tokens[i] = 'L' + start
                    positions.append(i)
            if positions:
                element.setAttribute('d', ' '.join(tokens))
                closures.append((element, positions))
        original(node, *args)
        for element, positions in closures:
            tokens = element.getAttribute('d').strip().split()
            for i in positions:
                tokens[i] = 'Z'
            element.setAttribute('d', ' '.join(tokens))

    module.correctCoordinates = compatible_coordinates
    try:
        exporter = module.SVGExporter(plot_item)
        exporter.parameters()['width'] = 1600
        exporter.export(str(path))
    finally:
        module.correctCoordinates = original
