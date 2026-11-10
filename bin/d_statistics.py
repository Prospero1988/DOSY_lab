"""Read-only snapshots and descriptive statistics; never call a fitter here."""
from copy import deepcopy
from datetime import datetime,timezone
from pathlib import Path
import csv,json,math,uuid
import numpy as np
from core import clean_json

SCIENTIFIC_NOTE = ('Statystyka obejmuje wyłącznie jawnie wybraną grupę/związek. '
    'SD opisuje rozrzut między ROI, a SE(D) niepewność pojedynczego fitu. '
    'ROI z jednego pomiaru mogą mieć skorelowane błędy i nie są niezależnymi powtórzeniami. '
    'SEM jest opisowe; nie obejmuje pełnej niepewności D. Bez automatycznego odrzucania odstających.')

def finite(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError, OverflowError):
        return None

def number(value, precision=5):
    value = finite(value)
    return f'{value:.{precision}e}' if value is not None else '—'

def quality(fit, error=''):
    fit = fit or {}
    d, se, r2 = (finite(fit.get(key)) for key in ('D_m2_s', 'D_standard_error_m2_s', 'R2'))
    if error or not fit.get('ok'):
        message = error or fit.get('message', 'Brak poprawnego dopasowania.')
        return False, 'BRAK FITU', ' '.join([message, *(fit.get('measurement_warnings') or [])])
    if d is None or d < 0:
        return False, 'BŁĘDNE D', 'D jest brakujące, ujemne lub nieskończone/NaN.'
    provided_warnings = fit.get('warnings') or []
    if isinstance(provided_warnings, str): provided_warnings = [provided_warnings]
    warnings = [str(w) for w in provided_warnings]
    warnings.extend(fit.get('measurement_warnings') or [])
    if se is None or se <= 0:
        warnings.append('Brak dodatniego, skończonego SE(D); niedostępne wagi 1/SE².')
    if r2 is None:
        warnings.append('R² niedostępne.')
    elif not 0 <= r2 <= 1:
        warnings.append('R² poza zakresem 0–1.')
    if d == 0:
        warnings.append('D na granicy 0.')
    return True, 'UWAGA' if warnings else 'OK', ' '.join(warnings) or 'Poprawny wynik; QC nie jest walidacją eksperymentu.'

def descriptive_statistics(values, uncertainties=None, weighted=False):
    """SD(ddof=1), descriptive SEM, and opt-in inverse-variance mean in SI units."""
    clean = [(i, finite(v)) for i, v in enumerate(values)]
    clean = [(i, v) for i, v in clean if v is not None]
    n = len(clean)
    out = dict(N=n, mean=None, median=None, SD=None, CV_percent=None, SEM=None,
               minimum=None, maximum=None, weighted_mean=None, weighted_SE_formal=None,
               weighted_status='Nie włączono; wynikiem domyślnym jest średnia arytmetyczna.')
    if not n:
        if weighted: out['weighted_status'] = 'Brak zaznaczonych poprawnych D.'
        return out
    a = np.array([v for _, v in clean], dtype=float)
    scale = float(np.max(np.abs(a))) or 1.
    normalized = a/scale
    mean = float(np.mean(normalized))
    out.update(mean=finite(mean*scale), median=finite(float(np.median(normalized))*scale),
               minimum=float(np.min(a)), maximum=float(np.max(a)))
    if n > 1:
        sd = float(np.std(normalized, ddof=1))
        out.update(SD=finite(sd*scale), SEM=finite(sd/math.sqrt(n)*scale),
                   CV_percent=finite(sd/abs(mean)*100) if mean else None)
    if weighted:
        uncertainties = uncertainties if uncertainties is not None else []
        se = [finite(uncertainties[i]) if i < len(uncertainties) else None for i, _ in clean]
        if any(v is None or v <= 0 for v in se):
            out['weighted_status'] = 'Nie obliczono: każde zaznaczone ROI musi mieć dodatnie, skończone SE(D).'
        else:
            smallest = min(se)
            weights = (smallest/np.array(se))**2  # proportional to inverse variance; no overflow
            out.update(weighted_mean=finite(float(np.sum(weights*normalized)/weights.sum())*scale),
                       weighted_SE_formal=finite(smallest/math.sqrt(float(weights.sum()))),
                       weighted_status='Opcjonalne 1/SE² po potwierdzeniu użytkownika. SE formalne zakłada niezależne, poprawne wariancje; korelacji nie uwzględniono.')
    return out

def _vector(value):
    try:
        a = np.asarray(value, dtype=float)
        return a if a.ndim == 1 else None
    except (ValueError, TypeError):
        return None

def diagnostic_points(snapshot, record):
    """Return measured rows without interpolating/dropping points or refitting.

    Explicit fit_mask/included_mask is respected. A malformed mask is reported,
    never guessed. Cached predictions include excluded experimental points too.
    """
    result, fit = record.get('result') or {}, record.get('fit') or {}
    b, g = _vector(snapshot.get('b')), _vector(snapshot.get('gradients'))
    y = _vector(result.get('corrected_integral'))
    if b is None or y is None or len(b) != len(y):
        return [], 'Brak zgodnych tablic rzeczywistych b i całek; punktów nie odtworzono.'
    n = len(b)
    notes = []
    mask = fit.get('fit_mask', fit.get('included_mask', result.get('fit_mask')))
    if mask is None:
        mask = np.ones(n, dtype=bool) if fit.get('ok') else None
    else:
        try:
            a = np.asarray(mask)
            if a.shape != (n,) or not np.isin(a, [False, True]).all():
                raise ValueError()
            mask = a.astype(bool)
        except (ValueError, TypeError):
            mask = None
            notes.append('Niezgodna maska dopasowania; udział punktów oznaczono jako nieznany.')
    def aligned(value):
        a = _vector(value)
        if a is None: return [None]*n
        if len(a) == n: return [finite(v) for v in a]
        if mask is not None and len(a) == int(mask.sum()):
            out = [None]*n
            for i, v in zip(np.flatnonzero(mask), a): out[int(i)] = finite(v)
            return out
        return [None]*n
    prediction = aligned(fit.get('prediction'))
    residual = aligned(fit.get('residuals'))
    snr = aligned(result.get('snr', result.get('SNR', result.get('signal_to_noise'))))
    rows = []
    for i in range(n):
        intensity, pred, resid = finite(y[i]), prediction[i], residual[i]
        if resid is None and intensity is not None and pred is not None:
            resid = finite(intensity-pred)  # identity using cached prediction, not a new fit
        rows.append({'index': i, 'g': finite(g[i]) if g is not None and len(g) == n else None,
                     'b': finite(b[i]), 'integral': intensity, 'prediction': pred, 'residual': resid,
                     'snr': snr[i], 'fit_included': bool(mask[i]) if mask is not None else None})
    missing = sum(r['b'] is None or r['integral'] is None for r in rows)
    if missing: notes.append(f'{missing} punktów ma niefinity b lub całkę; nie można ich narysować. Ich indeksy pozostają w eksporcie.')
    return rows, ' '.join(notes)


def build_snapshot(dataset,entries,experiment_id,label,**unused):
    records=[]
    for e in entries:
        fit=deepcopy(e.get('fit') or {})
        eligible,status,detail=quality(fit,e.get('error',''))
        bounds=e['bounds']
        records.append(dict(id=e['id'],name=e.get('name',e['id']),group=e.get('group',''),
            color=e.get('color','#67e8f9'),ppm=sum(bounds)/2,bounds=list(bounds),fit=fit,
            result=deepcopy(e.get('result') or {}),D=finite(fit.get('D_m2_s')),
            SE=finite(fit.get('D_standard_error_m2_s')),R2=finite(fit.get('R2')),
            eligible=eligible,qc=status,qc_detail=detail,analysis_error=e.get('error',''),
            use_background=dataset.options['baseline_enabled']))
    return dict(id=experiment_id,label=label,source_name=dataset.name,
        source_sha256=dataset.metadata['source_sha256'],updated_utc=datetime.now(timezone.utc).isoformat(),
        b=None if dataset.b is None else dataset.b.copy(),gradients=dataset.gradients.copy(),
        records=records,processing=deepcopy(dataset.options),
        experiment={'sequence':dataset.metadata['sequence']},diffusion_model=deepcopy(dataset.metadata['b_source']))


def groups(snapshot):return sorted({r.get('group','').strip() for r in snapshot['records'] if r.get('group','').strip()})


def included_records(snapshot,choices,group=None):
    group=snapshot.get('statistics_group','') if group is None else group
    return [r for r in snapshot['records'] if group and r.get('group')==group and r['eligible'] and choices.get(r['id'],True)]


def summarize(snapshot,choices,weighted=False,group=None):
    rows=included_records(snapshot,choices,group)
    return descriptive_statistics([r['D'] for r in rows],[r['SE'] for r in rows],weighted)


def export_statistics(directory,snapshot,choices,weighted=False,experiments=(),selected_id=None):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    summary=summarize(snapshot,choices,weighted)
    with (directory/'statistics_signals.csv').open('w',newline='',encoding='utf-8-sig') as h:
        w=csv.writer(h);w.writerow(['ID','name','group','ppm','D_m2_s','SE_D_m2_s','R2','QC','selected'])
        for r in snapshot['records']:w.writerow([r['id'],r['name'],r['group'],r['ppm'],r['D'],r['SE'],r['R2'],r['qc'],choices.get(r['id'],True)])
    with (directory/'statistics_summary.csv').open('w',newline='',encoding='utf-8-sig') as h:
        w=csv.writer(h);w.writerow(['group','statistic','value'])
        for k,v in summary.items():w.writerow([snapshot.get('statistics_group',''),k,v])
    with (directory/'experiment_comparison.csv').open('w',newline='',encoding='utf-8-sig') as h:
        w=csv.writer(h);w.writerow(['experiment_id','experiment','group','N','mean_D_m2_s','SD_m2_s','CV_percent','SEM_m2_s'])
        for e,c in experiments or [(snapshot,choices)]:
            for group in groups(e):
                s=summarize(e,c,group=group)
                w.writerow([e['id'],e['label'],group,s['N'],s['mean'],s['SD'],s['CV_percent'],s['SEM']])
    if selected_id:
        r=next((r for r in snapshot['records'] if r['id']==selected_id),None)
        if r:
            points,_=diagnostic_points(snapshot,r)
            with (directory/'selected_roi_points.csv').open('w',newline='',encoding='utf-8-sig') as h:
                w=csv.DictWriter(h,fieldnames=['index','g','b','integral','prediction','residual','snr','fit_included']);w.writeheader();w.writerows(points)
    payload=dict(schema='dosy_lab_statistics',version=1,experiment=snapshot,choices=choices,summary=summary)
    (directory/'statistics_snapshot.json').write_text(json.dumps(clean_json(payload),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    return summary


def import_statistics(path):
    data=json.loads(Path(path).read_text(encoding='utf-8'))
    if data.get('schema')!='dosy_lab_statistics' or data.get('version')!=1:raise ValueError('Nieobsługiwany zapis statystyki DOSY Lab.')
    e=data['experiment'];seen=set()
    for r in e['records']:
        if r['id'] in seen:raise ValueError('Powtórzone ID sygnału.')
        seen.add(r['id']);eligible,status,detail=quality(r['fit'],r.get('analysis_error',''))
        r.update(eligible=eligible,qc=status,qc_detail=detail,D=finite(r['fit'].get('D_m2_s')),
                 SE=finite(r['fit'].get('D_standard_error_m2_s')),R2=finite(r['fit'].get('R2')))
    e['id']=uuid.uuid4().hex;e['label']+=' · zapis'
    return e,{sid:bool(data.get('choices',{}).get(sid,True)) for sid in seen}
