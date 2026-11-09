"""Classical 1D DOSY processing, signed integration and nonlinear diffusion fits."""
from dataclasses import dataclass
from copy import deepcopy
import json
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid
from scipy.stats import t as student_t
from scipy.signal import find_peaks
from jeol import RawData,read_jdf,validate_b,reference_layer
from storage import require_extension


def default_options(raw):
    return dict(lb_Hz=0.,gaussian_FWHM_Hz=0.,zero_fill_factor=1,
                phase0_deg=0.,phase1_deg=0.,pivot_ppm=raw.metadata['offset_ppm'],
                ppm_shift=0.,reference_layer=reference_layer(raw),crop_filter=True,filter_delay_samples=0.,
                baseline_enabled=False,baseline_preview=True,baseline_degree=1,
                baseline_per_layer=True,baseline_regions=[])


def validate_options(options,nlayers=None):
    o=deepcopy(options)
    # Older projects already contain their full manual phase ramp.
    o.setdefault('filter_delay_samples',0.)
    o.setdefault('baseline_degree',1)
    o.setdefault('baseline_per_layer',True)
    for k in ('lb_Hz','gaussian_FWHM_Hz','phase0_deg','phase1_deg','pivot_ppm','ppm_shift','filter_delay_samples'):
        if not np.isfinite(float(o[k])):raise ValueError('Nieprawidłowy parametr '+k)
        o[k]=float(o[k])
    if min(o['lb_Hz'],o['gaussian_FWHM_Hz'])<0:raise ValueError('LB i szerokość Gaussa muszą być nieujemne.')
    if not 0<=o['filter_delay_samples']<=250:raise ValueError('Opóźnienie filtra musi być w zakresie 0–250 próbek.')
    for k,allowed in [('zero_fill_factor',(1,2,4,8)),('baseline_degree',(1,2,3,4,5))]:
        if o[k] not in allowed:raise ValueError('Nieprawidłowy parametr '+k)
    ref=o['reference_layer']
    if not isinstance(ref,int) or ref<0 or (nlayers is not None and ref>=nlayers):raise ValueError('Nieprawidłowa warstwa referencyjna.')
    for key in ('crop_filter','baseline_enabled','baseline_preview','baseline_per_layer'):
        if not isinstance(o[key],bool):raise ValueError('Nieprawidłowa flaga '+key)
    if o['baseline_enabled']:
        o['baseline_per_layer']=True
    for region in o['baseline_regions']:
        if len(region)!=2 or not np.isfinite(region).all() or region[0]==region[1]:raise ValueError('Nieprawidłowy obszar bazy.')
    return o


@dataclass
class Dataset:
    raw: RawData
    ppm: np.ndarray
    phased: np.ndarray
    baseline: np.ndarray
    spectrum: np.ndarray
    options: dict
    diagnostics: dict

    @property
    def b(self):return self.raw.b
    @property
    def gradients(self):return self.raw.gradients
    @property
    def metadata(self):return self.raw.metadata
    @property
    def name(self):return self.metadata['source_name']


def baseline_fit(ppm,real,regions,degree,reference=0,per_layer=True):
    mask=np.zeros(len(ppm),bool)
    for lo,hi in regions:
        lo,hi=sorted((lo,hi));mask|=(ppm>=lo)&(ppm<=hi)
    if np.count_nonzero(mask)<degree+2:raise ValueError('Wskaż więcej punktów bazowych; potrzeba co najmniej rząd+2.')
    span=np.ptp(ppm)
    if span<=0:raise ValueError('Zerowa szerokość osi ppm.')
    x=2*(ppm-(ppm.min()+ppm.max())/2)/span
    design=np.polynomial.legendre.legvander(x[mask],degree)
    condition=float(np.linalg.cond(design))
    if not np.isfinite(condition) or condition>1e10:raise ValueError('Niestabilna baza: rozmieść obszary szerzej po widmie.')
    full=np.polynomial.legendre.legvander(x,degree)
    if per_layer:coeff=np.linalg.lstsq(design,real[:,mask].T,rcond=None)[0].T
    else:
        one=np.linalg.lstsq(design,real[reference,mask],rcond=None)[0]
        coeff=np.repeat(one[None,:],len(real),axis=0)
    return coeff@full.T,dict(condition=condition,points=int(mask.sum()),coefficients=coeff.tolist())


def phase_radians(ppm,metadata,options):
    """One correction on the full acquired SW, increasing ppm, for every layer."""
    width=metadata['x_sweep_Hz']/(metadata['x_frequency_Hz']/1e6)
    carrier=metadata['offset_ppm']+options['ppm_shift']
    filter_ramp=360*options['filter_delay_samples']*(ppm-carrier)/width
    return np.deg2rad(filter_ramp+options['phase0_deg']+options['phase1_deg']*(ppm-options['pivot_ppm'])/width)


def process_fids(raw,options=None,roi_bounds=()):
    o=validate_options(options or default_options(raw),len(raw.fid))
    sw=raw.metadata['x_sweep_Hz'];freq=raw.metadata['x_frequency_Hz']/1e6
    n=raw.fid.shape[1];time=np.arange(n)/sw
    window=np.exp(-np.pi*o['lb_Hz']*time-(np.pi*o['gaussian_FWHM_Hz']*time)**2/(4*np.log(2)))
    nfft=n*o['zero_fill_factor']
    z=np.fft.fftshift(np.fft.fft(raw.fid*window[None,:],n=nfft,axis=1),axes=1)
    f=np.fft.fftshift(np.fft.fftfreq(nfft,1/sw))
    ppm=raw.metadata['offset_ppm']+f/freq+o['ppm_shift']
    z*=np.exp(1j*phase_radians(ppm,raw.metadata,o))[None,:]
    keep=np.ones(nfft,bool)
    if o['crop_filter'] and raw.metadata.get('digital_filter') and not raw.metadata.get('already_clipped'):
        half_band=raw.metadata['x_sweep_clipped_Hz']/2
        keep=(f>=-half_band)&(f<half_band)
    ppm=ppm[keep];z=z[:,keep]
    base=np.zeros(z.shape,float);diag={'nfft':nfft,'window':'exp(-pi*LB*t) * exp(-(pi*GB*t)^2/(4*ln(2)))',
                                     'baseline':None,'warnings':[]}
    if o['baseline_regions'] and (o['baseline_preview'] or o['baseline_enabled']):
        if o['baseline_enabled']:
            for a,b in o['baseline_regions']:
                for c,d in roi_bounds:
                    if max(min(a,b),min(c,d))<min(max(a,b),max(c,d)):
                        raise ValueError('Obszar bazy przecina ROI. Zmień maskę bazy przed zatwierdzeniem.')
        base,diag['baseline']=baseline_fit(ppm,z.real,o['baseline_regions'],o['baseline_degree'],
                                          o['reference_layer'],o['baseline_per_layer'])
        if min(min(x) for x in o['baseline_regions'])>ppm.min() or max(max(x) for x in o['baseline_regions'])<ppm.max():
            diag['warnings'].append('Baza jest ekstrapolowana poza skrajne wybrane obszary; sprawdź podgląd.')
    elif o['baseline_enabled']:raise ValueError('Najpierw zaznacz obszary linii bazowej.')
    corrected=z-base if o['baseline_enabled'] else z.copy()
    if not np.isfinite(corrected).all():raise ValueError('Processing wygenerował NaN/Inf.')
    for a in (ppm,z,base,corrected):a.setflags(write=False)
    return Dataset(raw,ppm,z,base,corrected,o,diag)


def starting_phase0(dataset):
    """Coarse positive absorption phase from coherent peaks in the reference only."""
    z=dataset.phased[dataset.options['reference_layer']]
    if len(z)<32:return None
    magnitude=abs(z)
    noise=np.median(abs(np.diff(z)))
    df=abs(dataset.ppm[1]-dataset.ppm[0])*dataset.metadata['x_frequency_Hz']/1e6
    peaks,_=find_peaks(magnitude,prominence=max(.04*magnitude.max(),12*noise),
                       distance=max(3,round(3.5/df)))
    peaks=peaks[(peaks>=.04*len(z))&(peaks<.96*len(z))]
    if len(peaks)<2:return None
    amplitudes=magnitude[peaks]
    weights=np.sqrt(amplitudes)
    direction=np.sum(weights*z[peaks]/amplitudes)/weights.sum()
    if abs(direction)<.70:return None
    return -float(np.angle(direction,deg=True))


def initial_dataset(raw):
    """New JDF only: visible filter correction and a shared editable PH0 seed."""
    options=default_options(raw)
    model=raw.metadata.get('filter_delay_model')
    if model is None:return process_fids(raw,options)
    options['filter_delay_samples']=model['samples']
    dataset=process_fids(raw,options)
    phase0=starting_phase0(dataset)
    if phase0 is None:return dataset
    options['phase0_deg']=round(phase0,6)
    return process_fids(raw,options)


def signed_integrals(ppm,spectra,bounds):
    ppm=np.asarray(ppm);z=np.asarray(spectra)
    if ppm.ndim!=1 or z.ndim!=2 or z.shape[1]!=len(ppm):raise ValueError('Niezgodne wymiary integracji.')
    if not np.isfinite(ppm).all() or not np.isfinite(z).all():raise ValueError('NaN/Inf w integracji.')
    order=np.argsort(ppm);x=ppm[order];y=z[:,order]
    if np.any(np.diff(x)<=0):raise ValueError('Oś ppm musi być ściśle monotoniczna.')
    lo,hi=sorted(map(float,bounds))
    if not np.isfinite([lo,hi]).all() or lo>=hi or lo<x[0] or hi>x[-1]:raise ValueError('ROI poza widmem lub zerowej szerokości.')
    if np.count_nonzero((x>=lo)&(x<=hi))<2:raise ValueError('ROI zawiera mniej niż dwa punkty widma.')
    inside=(x>lo)&(x<hi);xx=np.r_[lo,x[inside],hi]
    yy=np.column_stack(([np.interp(lo,x,row) for row in y],y[:,inside],
                        [np.interp(hi,x,row) for row in y]))
    return trapezoid(yy,x=xx,axis=1)


def integrate(dataset,bounds):
    pre=signed_integrals(dataset.ppm,dataset.phased.real,bounds)
    values=signed_integrals(dataset.ppm,dataset.spectrum.real,bounds)
    return {'raw_integral':pre,'corrected_integral':values,
            'unit':'a.u. ppm','method':'signed real trapezoid, interpolated ppm endpoints',
            'baseline_applied':bool(dataset.options['baseline_enabled']),
            'warnings':['Całki niedodatnie: sprawdź fazę, bazę i ROI.'] if np.any(values<=0) else []}


def fit_decay(b,intensity,fit_mask=None,sigma=None):
    y=np.asarray(intensity,dtype=float)
    bad=lambda message:{'ok':False,'message':message,'model':'I0 exp(-bD)',
                        'fit_mask':np.asarray(fit_mask,bool) if fit_mask is not None else np.ones(len(y),bool)}
    if b is None:return bad('Brak zweryfikowanej listy b; uzupełnij ją ręcznie.')
    try:b=validate_b(b,len(y))
    except ValueError as e:return bad(str(e))
    mask=np.ones(len(y),bool) if fit_mask is None else np.asarray(fit_mask)
    if mask.shape!=y.shape or not np.isin(mask,[False,True]).all():return bad('Nieprawidłowa maska punktów.')
    mask=mask.astype(bool);n=int(mask.sum())
    if not np.isfinite(y).all():return bad('Całki zawierają NaN/Inf; nie usuwamy ich automatycznie.')
    if n<3 or len(np.unique(b[mask]))<3:return bad('Potrzeba co najmniej 3 aktywnych, różnych wartości b.')
    if np.max(y[mask])<=0:return bad('Całki są niedodatnie; sprawdź fazowanie i bazę.')
    scale=float(np.max(abs(y[mask])));x=b/1e9;yn=y/scale;weights=None
    if sigma is not None:
        se=np.asarray(sigma,dtype=float)
        if se.shape!=y.shape or not np.isfinite(se[mask]).all() or (se[mask]<=0).any():return bad('Wagi wymagają dodatnich, skończonych niepewności aktywnych całek.')
        weights=se[mask]/scale
    model=lambda x,a,d:a*np.exp(-x*d)
    try:
        popt,pcov=curve_fit(model,x[mask],yn[mask],p0=[max(yn[mask].max(),.1),.5],
                           bounds=([0,0],[np.inf,np.inf]),sigma=weights,absolute_sigma=weights is not None,maxfev=20000)
        prediction=model(x,*popt)*scale;residual=y-prediction
        ss=float(np.sum(residual[mask]**2));den=float(np.sum((y[mask]-np.mean(y[mask]))**2))
        D=float(popt[1]*1e-9);Dse=float(np.sqrt(max(0,pcov[1,1]))*1e-9)
        warnings=[];ci=None;ase=float(np.sqrt(max(0,pcov[0,0]))*scale)
        if n<5:warnings.append('Mało aktywnych punktów; niepewności są słabo określone.')
        if np.any(y[mask]<=0):warnings.append('Uwzględniono niedodatnie całki bez przycinania i odrzucania.')
        reliable=np.isfinite(pcov).all() and np.linalg.cond(pcov)<1e12 and D>0
        if reliable:
            width=float(student_t.ppf(.975,n-2)*Dse)
            if D-width>0:ci=[D-width,D+width]
            else:warnings.append('Przybliżony przedział D obejmuje 0; nie podano fizycznego CI.')
        else:warnings.append('Dopasowanie słabo określone lub na granicy D=0; SE/CI wymagają ostrożności.')
        if den<=0 or (1-ss/den)<.95:warnings.append('Niska jakość monoeksponenty; obejrzyj reszty.')
        grid=np.linspace(float(b.min()),float(b.max()),240)
        return {'ok':True,'D_m2_s':D,'D_standard_error_m2_s':Dse if np.isfinite(Dse) else None,
                'D_CI95_m2_s':ci,'amplitude':float(popt[0]*scale),'amplitude_SE':ase if np.isfinite(ase) else None,
                'prediction':prediction,'residuals':residual,'R2':1-ss/den if den>0 else None,
                'RMSE':float(np.sqrt(ss/n)),'dof':n-2,'n_active':n,'fit_mask':mask,
                'b_grid':grid,'curve':model(grid/1e9,*popt)*scale,'model':'I0 exp(-bD)',
                'weighted':weights is not None,'warnings':warnings,
                'uncertainty_note':'SE i CI95: lokalne przybliżenie z kowariancji, warunkowe względem modelu i processingu; bez błędu kalibracji b i korelacji ROI.'}
    except (ValueError,RuntimeError,FloatingPointError,np.linalg.LinAlgError) as exc:return bad('Nieudany fit: '+str(exc))


def analyze(dataset,entries):
    rows=[]
    for entry in deepcopy(entries):
        try:
            result=integrate(dataset,entry['bounds'])
            fit=fit_decay(dataset.b,result['corrected_integral'],entry.get('fit_mask'),entry.get('sigma') if entry.get('weighted') else None)
            rows.append(dict(entry,result=result,fit=fit,error=''))
        except Exception as exc:rows.append(dict(entry,result={},fit={'ok':False,'message':str(exc)},error=str(exc)))
    return rows


def clean_json(value):
    if isinstance(value,np.ndarray):return clean_json(value.tolist())
    if isinstance(value,dict):return {str(k):clean_json(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean_json(v) for v in value]
    if isinstance(value,np.generic):return clean_json(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def save_profile(path,options):
    require_extension(path,'.list').write_text(json.dumps({'schema':'dosy_lab_processing','version':1,'options':clean_json(validate_options(options))},indent=2,ensure_ascii=False),encoding='utf-8')


def load_profile(path):
    data=json.loads(require_extension(path,'.list').read_text(encoding='utf-8'))
    if (not isinstance(data,dict) or data.get('schema')!='dosy_lab_processing'
            or data.get('version')!=1 or not isinstance(data.get('options'),dict)):
        raise ValueError('Nieobsługiwany profil processingu.')
    try:return validate_options(data['options'])
    except (KeyError, TypeError) as exc:raise ValueError('Niekompletne ustawienia profilu processingu.') from exc
