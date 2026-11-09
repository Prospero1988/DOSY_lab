"""Read-only details from the JDF and current processing."""
from html import escape

from PySide6 import QtWidgets


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
        symbol = names.get(code, f'[jedn. {code}]')
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
        self.setWindowTitle('Info o eksperymencie DOSY')
        self.resize(850,620)
        layout = QtWidgets.QVBoxLayout(self)
        tabs = QtWidgets.QTabWidget()
        layout.addWidget(tabs)
        self.experiment = QtWidgets.QTextBrowser()
        self.processing = QtWidgets.QTextBrowser()
        self.parameters = QtWidgets.QTextBrowser()
        tabs.addTab(self.experiment,'Eksperyment')
        tabs.addTab(self.processing,'Processing')
        tabs.addTab(self.parameters,'Parametry JDF')
        close = QtWidgets.QPushButton('Zamknij')
        close.clicked.connect(self.close)
        close.setToolTip('Wróć do analizy widma.')
        layout.addWidget(close)

    def set_dataset(self,dataset,integral_scale):
        if dataset is None:
            for tab in (self.experiment,self.processing,self.parameters):
                tab.setPlainText('Otwórz plik DOSY, aby zobaczyć informacje o pomiarze.')
            return
        raw = dataset.raw
        meta = dataset.metadata
        params = meta['parameters']
        source = meta['b_source']
        acquisition = [('Warstwy gradientowe',len(raw.fid)),('Punkty FID na warstwę',raw.fid.shape[1]),
                       ('Częstotliwość',f"{meta['x_frequency_Hz']/1e6:.6f} MHz"),
                       ('Szerokość akwizycji',f"{meta['x_sweep_Hz']:.6g} Hz"),
                       ('Offset',f"{meta['offset_ppm']:.6g} ppm")]
        for key,label in (
            ('actual_start_timestamp','Początek pomiaru (JDF)'),('solvent','Rozpuszczalnik'),
            ('temp_get','Temperatura odczytana'),('probe_attributes','Sonda'),('scans','Skany na warstwę'),
            ('x_prescans','Skany wstępne'),('x_acq_time','Czas akwizycji FID'),
            ('relaxation_delay','Opóźnienie relaksacyjne'),('repetition_time','Czas repetycji'),
            ('x_pulse','Impuls X'),('recvr_gain','Wzmocnienie odbiornika'),
        ):
            if key in params:
                acquisition.append((label,parameter_text(params[key])))
        diffusion = [('Źródło b',source['description'])]
        for key,label in (('grad_shape','Kształt gradientu'),('delta','δ'),('delta_large','Δ')):
            if key in params:
                diffusion.append((label,parameter_text(params[key])))
        for label,values,unit in (('Zakres g',dataset.gradients,'T/m'),('Zakres b',dataset.b,'s/m²')):
            if values is not None:
                finite = [v for v in values if float('-inf')<v<float('inf')]
                if finite:
                    diffusion.append((label,f'{min(finite):.6g} – {max(finite):.6g} {unit}'))
        self.experiment.setHtml(
            section('Plik i sekwencja',[('Nazwa',dataset.name),('Ścieżka',meta['source_path']),
                                       ('Sekwencja',meta['sequence'])])
            + section('Akwizycja',acquisition) + section('Kodowanie dyfuzji',diffusion))
        o = dataset.options
        delay = o['filter_delay_samples']
        delay_model = meta.get('filter_delay_model')
        delay_text = f'{delay:.6g} próbki'
        if delay and delay_model and abs(delay-delay_model['samples'])<1e-6:
            delay_text += ' (model JDF)'
        self.processing.setHtml(section('Widmo i integracja',[
            ('Rozmiar FFT',dataset.diagnostics['nfft']),('Warstwa referencyjna',o['reference_layer']+1),
            ('LB',f"{o['lb_Hz']:g} Hz"),('Gauss',f"{o['gaussian_FWHM_Hz']:g} Hz"),
            ('PH0',f"{o['phase0_deg']:.6g}°"),('PH1',f"{o['phase1_deg']:.6g}° / SW"),
            ('Pivot',f"{o['pivot_ppm']:.6g} ppm"),
            ('Korekcja opóźnienia filtra',delay_text),
            ('Korekcja bazy','Włączona' if o['baseline_enabled'] else 'Wyłączona'),
            ('Mnożnik wyświetlanych całek',f'{integral_scale:.9g}'),
        ]))
        self.parameters.setHtml(section('Pola zapisane w JDF',
                                        [(key,parameter_text(value)) for key,value in sorted(params.items())]))
