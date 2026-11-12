"""Read-only details from the JDF and current processing."""
from i18n import tr

from html import escape

from PySide6 import QtWidgets

from version import VERSION


def parameter_text(parameter):
    value = parameter.get('value')
    if value is None:
        return '—'
    text = f'{value:.9g}' if isinstance(value, float) else str(value)
    if parameter.get('decimal_scaler'):
        text += f" × 10^{parameter['decimal_scaler']}"
    names = {4:'°C', 13:'Hz', 19:'m', 26:'ppm', 28:'s', 31:'T'}
    prefixes = {-3:'G', -2:'M', -1:'k', 0:'', 1:'m', 2:'µ', 3:'n'}
    for unit in parameter.get('units', []):
        code, power, prefix = unit['base_code'], unit['power'], unit['prefix_code']
        if not code:
            continue
        symbol = names.get(code, tr(f'[jedn. {code}]'))
        text += ' ' + prefixes.get(prefix, f'10^{(-3*prefix)}·') + symbol
        if power!=1:
            text += f'^{power}'
    return text


def section(title, rows):
    cells = ''.join(f'<tr><td><b>{escape(label)}</b></td><td>{escape(str(value))}</td></tr>'
                    for label,value in rows)
    return f'<h3>{escape(title)}</h3><table cellspacing="6">{cells}</table>'


class ExperimentInfo(QtWidgets.QDialog):
    def __init__(self,parent=None):
        super().__init__(parent)
        self.setWindowTitle(tr('Info o eksperymencie DOSY'))
        self.resize(850,620)
        layout = QtWidgets.QVBoxLayout(self)
        tabs = QtWidgets.QTabWidget()
        layout.addWidget(tabs)
        self.experiment = QtWidgets.QTextBrowser()
        self.processing = QtWidgets.QTextBrowser()
        self.parameters = QtWidgets.QTextBrowser()
        tabs.addTab(self.experiment,tr('Eksperyment'))
        tabs.addTab(self.processing,'Processing')
        tabs.addTab(self.parameters,tr('Parametry JDF'))
        about = QtWidgets.QWidget()
        about_layout = QtWidgets.QVBoxLayout(about)
        about_layout.addWidget(QtWidgets.QLabel(f'DOSY Lab {VERSION}'))
        summary = QtWidgets.QLabel(tr('Analiza klasycznych eksperymentów dyfuzyjnych JEOL.'))
        summary.setWordWrap(True)
        about_layout.addWidget(summary)
        about_layout.addStretch()
        tabs.addTab(about,tr('O programie'))
        close = QtWidgets.QPushButton(tr('Zamknij'))
        close.clicked.connect(self.close)
        close.setToolTip(tr('Wróć do analizy widma.'))
        layout.addWidget(close)

    def set_dataset(self,dataset,integral_scale):
        if dataset is None:
            for tab in (self.experiment,self.processing,self.parameters):
                tab.setPlainText(tr('Otwórz plik DOSY, aby zobaczyć informacje o pomiarze.'))
            return
        raw = dataset.raw
        meta = dataset.metadata
        params = meta['parameters']
        source = meta['b_source']
        acquisition = [(tr('Warstwy gradientowe'),len(raw.fid)),(tr('Punkty FID na warstwę'),raw.fid.shape[1]),
                       (tr('Częstotliwość'),f"{meta['x_frequency_Hz']/1e6:.6f} MHz"),
                       (tr('Szerokość akwizycji'),f"{meta['x_sweep_Hz']:.6g} Hz"),
                       ('Offset',f"{meta['offset_ppm']:.6g} ppm")]
        for key,label in (
            ('actual_start_timestamp',tr('Początek pomiaru (JDF)')),('solvent',tr('Rozpuszczalnik')),
            ('temp_get',tr('Temperatura odczytana')),('probe_attributes',tr('Sonda')),('scans',tr('Skany na warstwę')),
            ('x_prescans',tr('Skany wstępne')),('x_acq_time',tr('Czas akwizycji FID')),
            ('relaxation_delay',tr('Opóźnienie relaksacyjne')),('repetition_time',tr('Czas repetycji')),
            ('x_pulse',tr('Impuls X')),('recvr_gain',tr('Wzmocnienie odbiornika')),
        ):
            if key in params:
                acquisition.append((label,parameter_text(params[key])))
        diffusion = [(tr('Źródło b'),tr(source['description']))]
        for key,label in (('grad_shape',tr('Kształt gradientu')),('delta','δ'),('delta_large','Δ')):
            if key in params:
                diffusion.append((label,parameter_text(params[key])))
        for label,values,unit in ((tr('Zakres g'),dataset.gradients,'T/m'),(tr('Zakres b'),dataset.b,'s/m²')):
            if values is not None:
                finite = [v for v in values if float('-inf')<v<float('inf')]
                if finite:
                    diffusion.append((label,f'{min(finite):.6g} – {max(finite):.6g} {unit}'))
        self.experiment.setHtml(
            section(tr('Plik i sekwencja'),[(tr('Nazwa'),dataset.name),(tr('Ścieżka'),meta['source_path']),
                                       (tr('Sekwencja'),meta['sequence']),(tr('Program JXP'),meta.get('pulse_program') or '—')])
            + section(tr('Akwizycja'),acquisition) + section(tr('Kodowanie dyfuzji'),diffusion))
        o = dataset.options
        delay = o['filter_delay_samples']
        delay_model = meta.get('filter_delay_model')
        delay_text = tr(f'{delay:.6g} próbki')
        if delay and delay_model and abs(delay-delay_model['samples'])<1e-6:
            delay_text += ' (model JDF)'
        filter_band = tr('Pełny FID mimo flagi clipped w JDF') if meta.get('full_fid_despite_clip_flag') else tr('Zgodnie z metadanymi JDF')
        self.processing.setHtml(section(tr('Widmo i integracja'),[
            (tr('Rozmiar FFT'),dataset.diagnostics['nfft']),(tr('Warstwa referencyjna'),o['reference_layer']+1),
            ('LB',f"{o['lb_Hz']:g} Hz"),('Gauss',f"{o['gaussian_FWHM_Hz']:g} Hz"),
            ('PH0',f"{o['phase0_deg']:.6g}°"),('PH1',f"{o['phase1_deg']:.6g}° / SW"),
            ('Pivot',f"{o['pivot_ppm']:.6g} ppm"),
            (tr('Korekcja opóźnienia filtra'),delay_text),
            (tr('Pasmo filtra'),filter_band),
            (tr('Korekcja bazy'),tr('Włączona') if o['baseline_enabled'] else tr('Wyłączona')),
            (tr('Mnożnik wyświetlanych całek'),f'{integral_scale:.9g}'),
        ]))
        self.parameters.setHtml(section(tr('Pola zapisane w JDF'),
                                        [(key,parameter_text(value)) for key,value in sorted(params.items())]))
