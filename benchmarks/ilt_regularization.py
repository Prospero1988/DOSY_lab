"""Benchmarks on known attenuation truth and the unchanged strychnine project."""
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bin'))
from adaptive_ilt import local_sweep_inputs
from core import analyze, clean_json, process_fids
from exports import read_project
from global_ilt import GlobalParameters
from ilt import ILTParameters
from ilt_regularization import SweepInput, SweepOptions, scan


def metrics(candidate, b_eval, truth, d):
    prediction = candidate.solution.predict(b_eval).reshape(truth.shape)
    error = np.linalg.norm((prediction-truth)/np.max(truth,axis=0))/np.sqrt(truth.size)
    result = dict(lambda_value=candidate.regularization, maxima_m2_s=candidate.maxima,
                  prediction_rms=error, weighted_residual=candidate.residual, penalty=candidate.penalty)
    if d is not None:
        cost = abs(np.log(np.asarray(candidate.maxima)[:,None]/d))
        i,j = linear_sum_assignment(cost)
        valid = cost[i,j]<=.20
        matched = int(valid.sum())
        result.update(false_maxima=len(candidate.maxima)-matched, missed=len(d)-matched,
                      relative_bias=list(np.asarray(candidate.maxima)[i[valid]]/d[j[valid]]-1))
    return result


def describe(result, b_eval, truth, d):
    cases = [metrics(c,b_eval,truth,d) for c in result.candidates if c.valid]
    chosen = next((c for c in cases if c['lambda_value']==result.recommendation),None)
    return dict(method=result.method,recommendation=chosen, interval=result.interval, noise=result.expected,
                warnings=result.warnings, converged=len(cases), trials=len(result.candidates),
                oracle=min(cases,key=lambda c:c['prediction_rms']) if cases else None,
                grid=[dict(kind=label,lambda_value=c.regularization,maxima=c.maxima,valid=c.valid)
                      for label,c in result.grid_checks])


def synthetic():
    rng = np.random.default_rng(81207)
    b = np.linspace(0.,1.6e10,32)
    b_eval = np.linspace(.13e9,16e9,151)
    output = {}
    cases = {
        'mono': ([3.8],[1.],.004),
        'separated': ([2.7,13.],[1.,.6],.004),
        'close': ([4.,5.],[1.,.55],.004),
        'weak': ([2.7,13.],[1.,.07],.004),
        'low_SNR': ([2.7,13.],[1.,.6],.15),
        'continuous': (np.geomspace(1.8,12.,101),np.full(101,.01),.005),
        'poor_b_coverage': ([2.7,13.],[1.,.6],.01),
        'heteroscedastic': ([2.7,13.],[1.,.6],np.linspace(.002,.012,len(b))),
    }
    for key,(d,amps,sigma) in cases.items():
        d = np.asarray(d)*1e-10
        sampled = b*.015 if key=='poor_b_coverage' else b
        sigma = np.broadcast_to(sigma,b.shape)
        y = np.exp(-np.outer(sampled,d))@amps+rng.normal(0,sigma)
        grid = ILTParameters(.2e-10,50e-10,50,regularization=1e-5)
        data = SweepInput(sampled,y,np.ones_like(b,dtype=bool),grid,sigma,sigma)
        truth = (np.exp(-np.outer(b_eval,d))@amps)[:,None]
        result = scan(data,SweepOptions(1e-10,1.,25,trust_noise=True,check_grid=True))
        case = describe(result,b_eval,truth,None if key=='continuous' else d)
        heuristic = scan(data,SweepOptions(1e-10,1.,25,method='corner'))
        case['corner_lambda'] = heuristic.recommendation
        if key in ('mono','heteroscedastic'):
            case['noise_scale_check'] = {}
            for factor in (.5,2.):
                changed = SweepInput(sampled,y,data.mask,grid,sigma,sigma*factor)
                diagnostic = scan(changed,SweepOptions(1e-10,1.,25,trust_noise=True))
                case['noise_scale_check'][str(factor)] = diagnostic.recommendation
        output[key] = case

    d = np.array([2.7,13.])*1e-10
    # Different ROI amplitudes, with both D values overlapping in the middle ROI.
    amps = np.array([[1.,.9,0.],[0.,.22,.05]])
    sigma = np.linspace(.002,.008,len(b))[:,None]*np.array([[1.,1.,.15]])
    ppm = np.linspace(3.9,4.2,401)
    profiles = np.exp(-.5*((ppm-np.array([4.04,4.06])[:,None])/.025)**2)
    profiles /= trapezoid(profiles,x=ppm,axis=1)[:,None]
    spectra = np.einsum('nd,dm,dp->nmp',np.exp(-np.outer(b,d)),amps,profiles)
    y = trapezoid(spectra,x=ppm,axis=2)+rng.normal(0,sigma)
    mask = np.ones_like(y,dtype=bool)
    grid = ILTParameters(.2e-10,50e-10,50,regularization=1e-6,order=0)
    source = SweepInput(b,y,mask,GlobalParameters(grid,1e-4),sigma,sigma)
    result = scan(source,SweepOptions(1e-10,1.,25,trust_noise=True,check_grid=True))
    output['joint_overlap_different_amplitudes'] = describe(result,b_eval,np.exp(-np.outer(b_eval,d))@amps,d)
    return output


def repeated():
    rng = np.random.default_rng(20261110)
    b = np.linspace(0.,1.6e10,32)
    grid = ILTParameters(.2e-10,50e-10,50)
    settings = SweepOptions(1e-10,1.,25,trust_noise=True)
    output = {}
    for key,d,amp,noise in (('mono',[3.8],[1.],.004),
                            ('separated',[2.7,13.],[1.,.6],.004),
                            ('close',[4.,5.],[1.,.55],.004),
                            ('low_SNR',[2.7,13.],[1.,.6],.15)):
        sigma = np.full(b.shape,noise)
        d = np.asarray(d)*1e-10
        truth = np.exp(-np.outer(b,d))@amp
        recommendations = []
        for _ in range(20):
            y = truth+rng.normal(0,sigma)
            data = SweepInput(b,y,np.ones(len(b),bool),grid,sigma,sigma)
            result = scan(data,settings)
            if result.recommendation is not None:
                chosen = next(c for c in result.candidates if c.regularization==result.recommendation)
                recommendations.append(metrics(chosen,b,truth[:,None],d))
        output[key] = dict(replicates=20,selected=len(recommendations),results=recommendations)
    return output


def real(project):
    raw, options, entries = read_project(project)
    data = process_fids(raw,options,[r['bounds'] for r in entries])
    rows = analyze(data,entries)
    groups = local_sweep_inputs(data,rows,5.,10.,1e-6,(),{})
    output = {}
    for ids,source in groups.items():
        result = scan(source,SweepOptions(1e-9,1e-3,13,check_grid=True))
        output[','.join(ids)] = dict(method=result.method,recommendation=result.recommendation,
            range=result.interval,fixed_mu=result.fixed_mu,warnings=result.warnings,
            trials=[dict(lambda_value=c.regularization,residual=c.residual,maxima=c.maxima,valid=c.valid)
                    for c in result.candidates])
    return dict(source_sha256=raw.metadata['source_sha256'],groups=output)


if __name__=='__main__':
    data = dict(seed=81207,repeated_seed=20261110,synthetic=synthetic(), repeated=repeated())
    if len(sys.argv)>1:
        data['strychnine'] = real(Path(sys.argv[1]))
    print(json.dumps(clean_json(data),ensure_ascii=False,indent=2))
