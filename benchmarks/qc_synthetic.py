"""Seeded, heteroscedastic controls for measurement QC and spectral rejection."""
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bin'))
from core import clean_json, fit_decay
from dosy_alignment import ilt_profiles
from roi_quality import qualify_roi


def run():
    previous = {}
    source = subprocess.check_output(['git','show','aedf2b8:bin/roi_quality.py'],cwd=ROOT)
    exec(compile(source,'v2.1 baseline roi_quality','exec'),previous)
    ppm = np.linspace(-1.,1.,81)
    bounds = (ppm[0],ppm[-1])
    b = np.linspace(0.,14.,24)*1e9
    mask = np.ones(len(b),bool)
    d = 3.8e-10
    decay = np.exp(-b*d)
    shape = 1/(1+1j*(ppm-.07)/.18)
    rng = np.random.default_rng(2102)
    repeats = 80
    result = {}
    for scenario, level, gradient_phase, overlap in (
        ('clean_low_noise',.003,0.,False), ('clean_late_noise',.04,0.,False),
        ('constant_phase',.012,0.,False), ('gradient_phase',.012,55.,False),
        ('overlap',.012,0.,True),
    ):
        counts = dict(old_flags=0,new_warnings=0,new_info=0,new_ok=0,phase_flags=0,shape_flags=0,
                      old_map_reject=0,new_map_reject=0,layers_used=[])
        for _ in range(repeats):
            sigma = level*np.linspace(.5,1.5,len(b))
            # Nearest-neighbour spectral correlation with unit per-quadrature variance.
            eps = rng.normal(size=(len(b),len(ppm)))+1j*rng.normal(size=(len(b),len(ppm)))
            eps = (eps+.4*np.roll(eps,1,axis=1))/np.sqrt(1.16)
            rotation = np.exp(1j*np.deg2rad((25 if scenario=='constant_phase' else 0.)
                                           + gradient_phase*b/b.max()))
            spectra = (decay*rotation)[:,None]*shape[None,:]
            if overlap:
                spectra += .7*np.exp(-b*1.9e-9)[:,None]/(1+1j*(ppm+.25)/.14)[None,:]
            spectra += sigma[:,None]*eps
            y = trapezoid(spectra.real,x=ppm,axis=1)
            row = dict(id='R1',bounds=bounds,result={'corrected_integral':y,'raw_integral':y},fit=fit_decay(b,y))
            dataset = SimpleNamespace(b=b,ppm=ppm,spectrum=spectra,options={'baseline_enabled':False})
            old = previous['assess_roi_profile'](dataset,bounds,mask)
            diag,findings = qualify_roi(dataset,bounds,row['fit'],row['result'],None,sigma,fit_decay)
            row['fit']['measurement_findings'] = findings
            warned = any(f['severity']=='warning' for f in findings)
            counts['old_flags'] += bool(old['flagged'])
            counts['new_warnings'] += warned
            counts['new_info'] += bool(findings) and not warned
            counts['new_ok'] += not findings
            counts['phase_flags'] += 'phase' in diag['flagged']
            counts['shape_flags'] += 'shape' in diag['flagged']
            counts['layers_used'].append(len(diag['diagnostic_layers']))
            # Oracle isolated component. Phase drift and overlap are incompatible
            # with this assignment; clean cases test false rejection.
            solve = SimpleNamespace(converged=True,diffusion=np.array([1e-10,d,2e-9]),
                weights=np.array([[0.],[1.],[0.]]), intervals=((0,3),),components=(1,),
                scales=np.array([max(abs(y))]),rms=np.array([0.]))
            ref = spectra[0,:-1].real
            target = spectra[:,:-1].real
            old_bad = np.linalg.norm(decay[:,None]*ref-target)>.15*np.linalg.norm(target)
            counts['old_map_reject'] += bool(old['flagged']) or old_bad
            _, rejected = ilt_profiles(b,ppm,spectra.real,0,[row],('R1',),solve,noise_sigma=sigma)
            counts['new_map_reject'] += bool(rejected)
        counts['median_layers_used'] = float(np.median(counts.pop('layers_used')))
        result[scenario] = counts
    return dict(seed=2102,replicates=repeats,result=result)


if __name__=='__main__':
    data = clean_json(run())
    (ROOT/'benchmarks'/'qc_maps_synthetic.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,indent=2))
