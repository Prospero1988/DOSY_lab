"""Reproduce the strychnine/CDCl3 QC and map-stage comparison."""
from hashlib import sha256
import json
import sys
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
from adaptive_ilt import adaptive_local_ilt
from core import analyze, clean_json, initial_dataset, process_fids
from exports import read_project
from d_statistics import quality
from dosy_alignment import alignment_plan, ilt_profiles
from dosy_map import fitted_diffusion, roi_bands
from global_ilt import GlobalParameters, solve_global_ilt
from ilt import ILTParameters, solve_ilt
from jeol import read_jdf


def run():
    raw = read_jdf(next((ROOT/'examples'/'Strychnina').glob('*.jdf')))
    dataset = initial_dataset(raw)
    reference = dataset.spectrum[0].real
    peaks, _ = find_peaks(reference,prominence=.01*reference.max(),distance=20)
    centers = dataset.ppm[peaks[(dataset.ppm[peaks]>0)&(dataset.ppm[peaks]<10)
                                &(reference[peaks]>.01*reference.max())]]
    windows = np.split(centers,np.where(np.diff(centers)>.055)[0]+1)
    entries = [dict(id=f'R{i:02}',bounds=(round(w[0]-.024,5),round(w[-1]+.024,5)))
               for i,w in enumerate(windows,1)]
    args = [arg for arg in sys.argv[1:] if not arg.startswith('--')]
    project_hash = sha256(Path(args[0]).read_bytes()).hexdigest() if args else None
    if args:
        raw, settings, entries = read_project(args[0])
        dataset = process_fids(raw,settings,[entry['bounds'] for entry in entries])
        reference = dataset.spectrum[settings['reference_layer']].real
    reference_layer = dataset.options['reference_layer']
    rows = analyze(dataset,entries)
    ids = tuple(r['id'] for r in rows)
    diffusion = np.array([fitted_diffusion(r) for r in rows])
    intensity = np.column_stack([r['result']['corrected_integral'] for r in rows])
    grid = ILTParameters(diffusion.min()/20,diffusion.max()*20,order=0,regularization=1e-6)
    joint = solve_global_ilt(dataset.b,intensity,GlobalParameters(grid))
    local = adaptive_local_ilt(dataset,rows,5.,10.,1e-6)
    withheld = {}
    options = {}
    new = '--before' not in sys.argv
    if new:
        from roi_quality import estimate_spectral_noise
        noise = estimate_spectral_noise(dataset.phased)
        options = dict(noise_sigma=noise, diagnostics={})
    profiles, rejected = ilt_profiles(dataset.b,dataset.ppm,dataset.spectrum.real,reference_layer,rows,ids,joint,
                                      withheld=withheld,**options)
    align = alignment_plan(dataset.b,rows,5.,(),False)
    shared = alignment_plan(dataset.b,rows,5.,(),True)
    modes = {
        'original': ({},{},set()),
        'aligned': ({sid:g.representative for g in align for sid in g.ids},{},set()),
        'shared': ({sid:g.shared for g in shared if not g.incompatible for sid in g.ids},{},set()),
        'global': ({},profiles,(set(ids)-profiles.keys()) if new else set(withheld)),
        'adaptive': ({},local.profiles,(set(ids)-local.profiles.keys()) if new else set(local.withheld)),
    }
    matrices = {}
    for mode,(overrides,mode_profiles,excluded) in modes.items():
        reconstructed = new and mode in ('global','adaptive')
        display_options = dict(reconstructed_only=True) if reconstructed else {}
        bands, overlap = roi_bands(dataset.ppm,reference,rows,() if reconstructed else excluded,
                                   overrides=overrides,profiles=mode_profiles,**display_options)
        positive = {band.sid for band in bands if np.max(band.intensity)>0}
        matrices[mode] = dict(matrix=sorted({band.sid for band in bands}),visible=sorted(positive),
                              overlap_mask=overlap)
    results = []
    for i,r in enumerate(rows):
        f = r['fit']
        ind = solve_ilt(dataset.b,intensity[:,i],grid)
        sid = r['id']
        results.append(dict(id=sid, bounds=r['bounds'], spectrum_present=bool(np.any(
            reference[(dataset.ppm>=r['bounds'][0])&(dataset.ppm<=r['bounds'][1])]>0)),
            integrated=True, D_m2_s=f['D_m2_s'], R2=f['R2'], RMS=f['RMSE'],
            QC=quality(f)[1], fit_warnings=f['warnings'],
            qc_detail=quality(f)[2], diagnostic=r['result']['profile_diagnostics'],
            independent_ILT_D=[ind.diffusion[j] for j in ind.peaks], independent_ILT_rms=ind.rms,
            global_optimization=True, global_converged=joint.converged,
            global_candidate_D=[joint.diffusion[k] for k,q in zip(joint.components,joint.fractions[:,i]) if q>0],
            global_rejection=withheld.get(sid), global_assigned=sid in profiles,
            adaptive_optimization=(local.diagnostics.get(sid,{}).get('optimized',False)) if new else any(sid in g.ids for g in local.groups),
            adaptive_rejection=local.withheld.get(sid), adaptive_assigned=sid in local.profiles,
            global_profile_diagnostic=options.get('diagnostics',{}).get(sid),
            adaptive_profile_diagnostic=getattr(local,'diagnostics',{}).get(sid),
            matrix_and_visible={mode: [sid in record['matrix'],sid in record['visible']]
                                for mode,record in matrices.items()}))
    return clean_json(dict(sha256=raw.metadata['source_sha256'], options=dataset.options,
        project_sha256=project_hash,
        peak_rule=('saved project ROI definitions' if args else
                   'positive peaks >1% maximum, prominence >=1%, 0–10 ppm, merge gap 0.055 ppm, pad 0.024 ppm'),
        floor_percent=0, gain=1, reference=reference_layer+1, rows=results, modes=matrices,
        global_maxima=joint.diffusion[list(joint.components)], global_rejected=rejected,
        local_unresolved=local.unresolved))


if __name__ == '__main__':
    target = 'qc_maps_before.json' if '--before' in sys.argv else 'qc_maps_after.json'
    (ROOT/'benchmarks'/target).write_text(json.dumps(run(),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
