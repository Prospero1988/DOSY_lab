"""Compare identical ROI-map settings against archived v2.0.0 and the 2.1 baseline."""
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]


def measure(project, source):
    sys.path.insert(0, str(source))
    import storage

    # Archived implementations contain only bin. Resolve relocated project inputs
    # against the same real workspace for every version.
    storage.DATA, storage.WORK, storage.EXAMPLES = (ROOT/name for name in ('data','work','examples'))
    from dataclasses import asdict
    import numpy as np
    from scipy.integrate import trapezoid
    import adaptive_ilt
    from core import analyze, clean_json, process_fids
    from dosy_alignment import alignment_plan, ilt_profiles
    from dosy_map import roi_bands
    from exports import read_project
    from global_ilt import GlobalParameters, solve_global_ilt
    from ilt import ILTParameters, solve_ilt

    raw, options, entries = read_project(project)
    data = process_fids(raw, options, [row['bounds'] for row in entries])
    rows = analyze(data,entries)
    ids = tuple(row['id'] for row in rows)
    b, ppm, spectra = data.b, data.ppm, data.spectrum.real
    ref = data.options['reference_layer']
    y = np.column_stack([row['result']['corrected_integral'] for row in rows])
    mask = np.column_stack([row['fit']['fit_mask'] for row in rows])
    d = np.array([row['fit']['D_m2_s'] for row in rows])
    grid = ILTParameters(d.min()/20,d.max()*20,order=0,regularization=1e-6)
    joint = solve_global_ilt(b,y,GlobalParameters(grid),include=mask)
    diagnostic, withheld = {}, {}
    args = (b,ppm,spectra,ref,rows,ids,joint)
    solutions = []

    def save_local_solution(*args, **kwargs):
        result = solve_global_ilt(*args, **kwargs)
        solutions.append(dict(diffusion=result.diffusion, weights=result.weights, scales=result.scales,
                              rms=result.rms, converged=result.converged))
        return result

    adaptive_ilt.solve_global_ilt = save_local_solution
    local = adaptive_ilt.adaptive_local_ilt(data,rows,5.,10.,1e-6)
    strict = hasattr(local,'diagnostics')
    if strict:
        from roi_quality import estimate_spectral_noise
        accepted, _ = ilt_profiles(*args, withheld=withheld, diagnostics=diagnostic,
                                   noise_sigma=estimate_spectral_noise(data.phased))
    else:
        accepted, _ = ilt_profiles(*args)
    group_plan = alignment_plan(b,rows,5.,(),True)
    shared = {sid:group.shared for group in group_plan if not group.incompatible for sid in group.ids}

    def visible(profiles=None, overrides=None, reconstructed=False, exclude=()):
        extras = dict(reconstructed_only=True) if reconstructed else {}
        bands, _ = roi_bands(ppm,spectra[ref],rows,exclude,profiles=profiles,overrides=overrides,**extras)
        return {sid: [dict(D_m2_s=band.diffusion,area=trapezoid(np.maximum(band.intensity,0.),x=band.ppm))
                      for band in bands if band.sid==sid and np.any(band.intensity>0)]
                for sid in ids if any(band.sid==sid and np.any(band.intensity>0) for band in bands)}

    modes = dict(mono=visible(), shared=visible(overrides=shared),
                 global_filtered=visible(accepted,reconstructed=strict,exclude=() if strict else withheld),
                 adaptive_filtered=visible(local.profiles,reconstructed=strict))
    individual = {}
    for index, row in enumerate(rows):
        params = ILTParameters(d[index]/20,d[index]*20)
        individual[row['id']] = solve_ilt(b,y[:,index],params,include=mask[:,index])
    extra = {}
    if hasattr(local,'exploratory'):
        from dosy_alignment import exploratory_profiles, individual_map_solution
        recovered, records = exploratory_profiles(*args)
        independent, independent_records, independent_filtered = {}, {}, {}
        noise = estimate_spectral_noise(data.phased)
        for sid,result in individual.items():
            inputs = (b,ppm,spectra,ref,rows,(sid,),individual_map_solution(result))
            profiles, details = exploratory_profiles(*inputs)
            independent.update(profiles)
            independent_records.update(details)
            filtered, _ = ilt_profiles(*inputs,noise_sigma=noise)
            independent_filtered.update(filtered)
        modes.update(global_exploratory=visible(recovered,reconstructed=True),
                     adaptive_exploratory=visible(local.exploratory,reconstructed=True),
                     individual_exploratory=visible(independent,reconstructed=True),
                     individual_filtered=visible(independent_filtered,reconstructed=True))
        extra = dict(global_exploratory=records, adaptive_exploratory=local.exploratory_diagnostics,
                     individual_exploratory=independent_records)
    return clean_json(dict(ids=ids, mono_D=d, b=b, ppm=ppm, y=y, reference=ref+1,
        source_sha256=raw.metadata['source_sha256'], modes=modes, global_fit=asdict(joint),
        individual_fits={sid:asdict(result) for sid,result in individual.items()},
        adaptive_groups=[asdict(group) for group in local.groups], local_joint_fits=solutions,
        global_filtered_diagnostics=diagnostic,
        global_assigned=sorted(accepted), adaptive_assigned=sorted(local.profiles), **extra))


def run(project):
    versions = {}
    env = {**os.environ,'OPENBLAS_NUM_THREADS':'1','QT_QPA_PLATFORM':'offscreen','PYTHONIOENCODING':'utf-8'}
    for ref in ('v2.0.0','e601223','working'):
        with tempfile.TemporaryDirectory(prefix='dosy-map-comparison-') as scratch:
            directory = Path(scratch)
            if ref == 'working':
                source = ROOT/'bin'
            else:
                archive = subprocess.check_output(['git','archive','--format=zip',ref,'bin'],cwd=ROOT)
                with ZipFile(io.BytesIO(archive)) as bundle:
                    bundle.extractall(directory)
                source = directory/'bin'
            out = directory/'result.json'
            subprocess.run([sys.executable,__file__,str(project),str(source),str(out)],
                           cwd=ROOT,env=env,check=True)
            versions[ref] = json.loads(out.read_text(encoding='utf8'))
    return versions


def summarize(versions):
    import numpy as np

    current = versions['working']
    baseline = versions['e601223']
    for key in ('b','ppm','y','mono_D','global_fit','individual_fits','local_joint_fits'):
        if current[key] != baseline[key]:
            raise AssertionError(f'Visualization changed the baseline numerical results: {key}')

    def delta(a, b):
        return float(np.max(abs(np.asarray(a)-b)))

    version_results = {}
    for ref, version in versions.items():
        modes = version['modes']
        diffs = {key:delta(version[key],current[key]) for key in ('b','ppm','y','mono_D')}
        diffs.update({f'global_{key}':delta(version['global_fit'][key],current['global_fit'][key])
                      for key in ('diffusion','weights','scales','rms')})
        diffs['individual_weights'] = max(delta(fit['weights'],current['individual_fits'][sid]['weights'])
                                          for sid,fit in version['individual_fits'].items())
        diffs['local_joint_weights'] = max(delta(a['weights'],b['weights'])
                                          for a,b in zip(version['local_joint_fits'],current['local_joint_fits'],strict=True))
        version_results[ref] = dict(visible={k:len(v) for k,v in modes.items()}, max_absolute_difference=diffs,
                                    global_native_fallback=sorted(set(modes['global_filtered'])-set(version['global_assigned'])),
                                    adaptive_native_fallback=sorted(set(modes['adaptive_filtered'])-set(version['adaptive_assigned'])))
    projections = {}
    for mode in ('global_exploratory','adaptive_exploratory','individual_exploratory'):
        profiles = current['modes'][mode]
        for sid, components in profiles.items():
            limit = current[mode][sid]['model_reference_area']
            area = sum(component['area'] for component in components)
            measured = sum(component['area'] for component in current['modes']['mono'][sid])
            if area > measured + 1e-9 or area > sum(limit) + 1e-9:
                raise AssertionError(f'Unaccounted reference area: {mode} {sid}')
        projections[mode] = dict(visible=sorted(profiles), S2=current[mode]['S2'])
    return dict(source_sha256=current['source_sha256'], reference=current['reference'],
                versions=version_results, exploratory=projections)


if __name__ == '__main__':
    project = Path(sys.argv[1]).resolve()
    if len(sys.argv)==4:
        Path(sys.argv[3]).write_text(json.dumps(measure(project,Path(sys.argv[2])),ensure_ascii=False),encoding='utf8')
    else:
        results = run(project)
        print(json.dumps(summarize(results) if '--summary' in sys.argv else results,ensure_ascii=False,indent=2))
