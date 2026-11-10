"""Independent audit fixtures for DOSY Lab v2.0.0. The release is the system under test."""
from copy import deepcopy
from dataclasses import replace
import json
import os
from pathlib import Path
import platform
import subprocess
import tempfile
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np
import scipy
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.linalg import expm
from scipy.optimize import linear_sum_assignment, lsq_linear, minimize_scalar

import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from adaptive_ilt import adaptive_local_ilt
from core import analyze, default_options, fit_decay, initial_dataset, process_fids, signed_integrals
from dosy_alignment import alignment_plan, ilt_profiles
from dosy_map import roi_bands
from exports import export_bundle, read_integrals, read_project, save_integrals, save_project
from global_ilt import GlobalParameters, solve_global_ilt
from ilt import ILTParameters, solve_ilt
from jeol import RawData, factory_b, read_jdf

SEED = 20261110
DUNIT = 1e-10
GAMMA = 2 * np.pi * 42.577478518e6  # independent proton angular gyromagnetic ratio
B = np.r_[0., np.geomspace(.02, 35., 48)]  # dimensionless b*1e-9
TOL = .10  # matching threshold for discrete populations, in absolute log D
OUT = ROOT / "audit" / "results"


def decay(b9, d10, amplitude):
    return np.exp(-np.outer(b9, np.asarray(d10) / 10)) @ np.asarray(amplitude)


def reference_mono(b9, y, sigma=None):
    """Profile out positive I at min(b), bracket all sampled minima in log D."""
    scale = np.max(np.abs(y))
    ys = y / scale
    w = np.ones(len(y)) if sigma is None else 1 / np.asarray(sigma)
    w /= np.linalg.norm(w)
    x = (b9 - min(b9)) / 10

    def profiled(d10):
        k = w * np.exp(-x * d10)
        a = max(0., np.dot(w * ys, k) / np.dot(k, k))
        return np.linalg.norm(w * ys - a * k) ** 2, a

    grid = np.geomspace(1e-6, 1e5, 180)
    values = [profiled(d)[0] for d in grid]
    candidates = [0., grid[-1]]
    for i in range(1, len(grid) - 1):
        if values[i] <= values[i-1] and values[i] < values[i+1]:
            optimum = minimize_scalar(lambda z: profiled(np.exp(z))[0],
                                      bounds=np.log(grid[[i-1, i+1]]), method="bounded",
                                      options={"xatol": 1e-12})
            candidates.append(np.exp(optimum.x))
    best = min(candidates, key=lambda d: profiled(d)[0])
    _, a = profiled(best)
    return best, scale * a * np.exp(-x * best)


def reference_ilt(b9, y, params, sigma=None):
    """Independent bin centers and BVLS active set on the augmented design."""
    n = params.points
    ends = np.array([params.d_min, params.d_max])
    if params.spacing == "log":
        edges = np.linspace(*np.log(ends), n+1)
        d = np.exp((edges[1:] + edges[:-1]) / 2)
    else:
        edges = np.linspace(*ends, n+1)
        d = (edges[1:] + edges[:-1]) / 2
    kernel = np.exp(-b9[:, None] * 1e9 * d)
    scale = max(abs(y))
    w = np.ones(len(y)) if sigma is None else 1 / np.asarray(sigma)
    w /= np.linalg.norm(w)
    penalty = np.eye(n)
    if params.order == 2:
        penalty = (penalty[:-2] - 2 * penalty[1:-1] + penalty[2:]) * (.05 * n)**2
    design = np.vstack([kernel*w[:, None], np.sqrt(params.regularization*n)*penalty])
    target = np.r_[w*y/scale, np.zeros(len(penalty))]
    ref = lsq_linear(design, target, bounds=(0., np.inf), method="bvls", tol=1e-12)
    return d, kernel, design, target, ref


def global_certificate(b9, y, params, result, sigma=None, mask=None):
    """Independent primal/Fenchel dual and active-row KKT residual."""
    mask = np.ones_like(y, dtype=bool) if mask is None else mask
    scale = np.max(np.where(mask, np.abs(y), 0.), axis=0)
    weight = mask.astype(float) if sigma is None else np.where(mask, 1/sigma, 0.)
    weight /= np.linalg.norm(weight, axis=0)
    kernel = np.exp(-np.outer(b9 * 1e9, result.diffusion))
    target = weight * y / scale
    f = result.weights
    alpha = params.grid.points * params.grid.regularization
    beta = params.sparsity * np.sqrt(y.shape[1])
    residual = weight * (kernel @ f) - target
    row_norm = np.linalg.norm(f, axis=1)
    primal = np.sum(residual**2) + alpha*np.sum(f*f) + beta*sum(row_norm)
    v = 2*residual
    linear = kernel.T @ (weight*v)
    dual = -np.sum(v*target) - np.sum(v*v)/4
    dual -= np.sum(np.maximum(np.linalg.norm(np.maximum(-linear, 0), axis=1)-beta, 0)**2)/(4*alpha)
    active = row_norm > 1e-12
    grad = linear[active] + 2*alpha*f[active] + beta*f[active]/row_norm[active, None]
    violations = np.where(f[active] > 1e-12, abs(grad), np.maximum(-grad, 0))
    inactive = np.maximum(np.linalg.norm(np.maximum(-linear[~active], 0), axis=1)-beta, 0)
    return dict(primal=primal, dual=dual, gap=primal-dual,
                objective_difference=primal-result.objective,
                kkt_max=max(np.max(violations, initial=0.), np.max(inactive, initial=0.)),
                scale_error=np.max(np.abs(result.scales-scale)))


def match_peaks(estimated, known):
    if not len(estimated) or not len(known):
        return dict(matched=0, extra=len(estimated), missed=len(known), matched_bias_pct=[])
    errors = np.abs(np.log(np.asarray(estimated)[:, None]/np.asarray(known)))
    # Penalize pairs outside tolerance before minimizing remaining distances.
    i, j = linear_sum_assignment(np.where(errors <= TOL, errors, 10.))
    valid = errors[i,j] <= TOL
    bias = 100*(np.asarray(estimated)[i[valid]]/np.asarray(known)[j[valid]]-1)
    return dict(matched=len(bias), extra=len(estimated)-len(bias), missed=len(known)-len(bias),
                matched_bias_pct=bias.tolist())


def ilt_record(result, known):
    peaks = (result.diffusion[list(result.peaks)] / DUNIT).tolist()
    # Full bin-mass distribution, summarized as a central log-grid interval.
    quantiles = []
    if result.weights.sum() > 0:
        quantiles = (np.interp([.05,.95], np.cumsum(result.weights)/result.weights.sum(),
                               result.diffusion)/DUNIT).tolist()
    return dict(peaks=peaks, rms=result.rms, rms_normalized=result.rms/result.scale,
                width90=quantiles, rank=result.effective_rank, warnings=result.warnings,
                **(match_peaks(peaks, known) if known is not None else {}))


def benchmark_cases():
    rng = np.random.default_rng(SEED)
    grid = ILTParameters(.25*DUNIT, 40*DUNIT, points=100)
    x = B
    cases = {}
    for label, d, a in (
        ("A_single", [4], [1]), ("B_separated", [2, 16], [.6,.4]),
        ("C_close", [4,4.5], [.5,.5]), ("D_three", [1,5,20], [.25,.5,.25]),
    ):
        cases[label] = (x, decay(x,d,a), d, None)
    for ratio in (1,10,100,1000):
        d = [2,16]
        y = decay(x,d,[1.,1/ratio]) + rng.normal(0,.001,len(x))
        cases[f"E_ratio_{ratio}"] = (x,y,d,None)
    se = np.linspace(.002,.02,len(x))
    cases["F_heteroscedastic"] = (x,decay(x,[4],[1])+rng.normal(0,se),[4],se)
    cases["G_phase_baseline"] = (x,decay(x,[4],[1])*np.cos(np.linspace(0,.9,len(x)))+.08,[4],None)
    for tag, sample in (("short_span", np.linspace(0.,1.,12)), ("late_only", np.linspace(100.,180.,12)),
                        ("clustered", np.r_[np.linspace(.1,.12,7),12.,12.1])):
        d = [.4] if tag == "late_only" else [4]
        cases["H_"+tag] = (sample,decay(sample,d,[1]),d,None)
    cases["I_overlap"] = (x,decay(x,[3.8,18.8],[1,.35]),[3.8,18.8],None)
    cases["L_exchange_fast_limit"] = (x,decay(x,[7.6],[1]),[7.6],None)
    # Two-site stationary 60:40 exchange in the ideal narrow-pulse limit.
    for rate in (0., 10., 1000.):
        transfer = rate*np.array([[-.4,.6],[.4,-.6]])
        population = np.array([.6,.4])
        diffusivities = np.diag([2.,16.])
        duration = .1
        signal = np.array([(expm(duration*transfer - bb*.1*diffusivities) @ population).sum() for bb in x])
        cases[f"L_exchange_{rate:g}_per_s"] = (x,signal,[2.,16.] if rate==0 else None,None)
    continuum = np.geomspace(.5,28,5001)
    density = np.exp(-.5*(np.log(continuum/4)/.38)**2)
    cases["M_continuous"] = (x,decay(x,continuum,density/density.sum()),None,None)
    cases["N_noise_only"] = (x,rng.normal(0,.01,len(x)),[],None)
    cases["O_solvent_overlap"] = (x,decay(x,[3.8,19],[1,40]),[3.8,19],None)
    cases["P_negative_tail"] = (x,decay(x,[4],[1])-.14,[4],None)
    out = {}
    for key,(bx,y,truth,sd) in cases.items():
        fit = fit_decay(bx*1e9,y,sigma=sd)
        ref_d, ref_y = reference_mono(bx,y,sd)
        alt = {}
        for lam in (1e-8, 1e-5, 1e-2):
            result = solve_ilt(bx*1e9,y,replace(grid,regularization=lam),sigma=sd)
            alt[f"{lam:g}"] = ilt_record(result,truth)
        out[key] = dict(b9=bx, y=y, truth=truth, sigma=sd, mono=mono_record(fit), ilt=alt,
                        ref_D10=ref_d, ref_rms=np.sqrt(np.mean((y-ref_y)**2)),
                        mono_signed_prediction=fit.get("prediction"))
    return out


def mono_record(fit):
    if not fit.get("ok"):
        return dict(ok=False, message=fit["message"])
    ci = fit["D_CI95_m2_s"]
    return dict(ok=True, D10=fit["D_m2_s"]/DUNIT, se10=fit["D_standard_error_m2_s"]/DUNIT
                if fit["D_standard_error_m2_s"] is not None else None, rms=fit["RMSE"],
                ci10=[c/DUNIT for c in ci] if ci else None, R2=fit["R2"], warnings=fit["warnings"])


def monte_carlo():
    rng = np.random.default_rng(SEED+1)
    out = {}
    x = np.linspace(0,16,24)
    truth = 4.
    clean = decay(x,[truth],[1.])
    for noise in (.005,.025,.075):
        sigma = noise*np.linspace(.45,1.4,len(x))
        records = {kind:[] for kind in ("unweighted","known_sigma","observed_above_2sigma")}
        for _ in range(160):
            y = clean + rng.normal(0,sigma)
            for name, mask, sd in (
                ("unweighted", None, None), ("known_sigma", None, sigma),
                ("observed_above_2sigma", y>2*sigma, sigma),
            ):
                fit = fit_decay(x*1e9,y,fit_mask=mask,sigma=sd)
                records[name].append(mono_record(fit))
        for name, rows in records.items():
            good = [r for r in rows if r["ok"]]
            estimates = np.array([r["D10"] for r in good])
            intervals = [r["ci10"] for r in good if r["ci10"]]
            out[f"{noise:g}_{name}"] = dict(n=len(good), mean=estimates.mean(), bias_pct=100*(estimates.mean()/truth-1),
                sd=estimates.std(ddof=1), median_se=np.median([r["se10"] for r in good]),
                coverage=np.mean([lo <= truth <= hi for lo,hi in intervals]) if intervals else None,
                reported_CI=len(intervals), D10=estimates)
    return out


def waveform(u, shape):
    if shape == "SQUARE":
        return np.ones_like(u)
    if shape == "SINE":
        return np.sin(np.pi*u)
    if shape == "SQUARE_SINE":
        return np.sin(np.pi*u)**2
    ramp = np.minimum(u, 1-u)/.1
    if shape == "TRAPEZOID":
        return np.minimum(ramp, 1.)
    if shape == "S_RECTANGLE":
        return (1-np.cos(np.pi*np.minimum(ramp,1.)))/2
    raise ValueError(shape)


def integrate_sequence(shape, delta, big_delta, tau1, tau2, rf90=0.):
    """Construct signed encoding lobes from JEOL pulse ordering.

    tau2=None selects STE. Spoilers during longitudinal storage have p=0.
    RF is a gradient-free interval with instantaneous toggling in its center.
    """
    length = delta/2
    u = np.linspace(0.,1.,20001)
    f = waveform(u, shape)
    moment = length*cumulative_trapezoid(f, x=u, initial=0.)
    q = b = m1 = time = 0.
    def wait(duration):
        nonlocal b, m1, time
        if duration < 0:
            raise ValueError("Negative delay in independent timeline.")
        b += q*q*duration
        m1 += q*duration
        time += duration
    def lobe(sign):
        nonlocal q, b, m1, time
        ramp = q + sign*moment
        b += length*trapezoid(ramp**2, x=u)
        m1 += length*trapezoid(ramp,x=u)
        q = ramp[-1]
        time += length
    def bipolar(sign, tau, repeats=1):
        wait(rf90)
        for n in range(2*repeats):
            lobe(sign)
            wait(tau + (2*rf90 if n == repeats-1 else 0.))
        wait(rf90)
    if tau2 is None:
        bipolar(1, tau1)
        wait(big_delta-4*rf90-delta-2*tau1)
        bipolar(-1, tau1)
    else:
        storage = big_delta/2-4*rf90-delta-2*tau1
        bipolar(1, tau1)
        wait(storage)
        bipolar(-1, tau2, 2)
        wait(storage-tau1+tau2)
        bipolar(1, tau1)
    return dict(coefficient_s=b/delta**2, m0_s=q, integral_q_s2=m1,
                elapsed_s=time, area_factor=(2*moment[-1]/delta)**2)


def metadata_number(params, key):
    """Decode SI from record exponents without calling the production converter."""
    record = params[key]
    return record["value"] * 10.**(record["decimal_scaler"]
        -3*sum(unit["prefix_code"]*unit["power"] for unit in record["units"]))


def physical_b():
    out = {"gamma_ref":GAMMA, "waveforms":{}, "jdf":{}}
    for path in sorted((ROOT/"examples").rglob("*.jdf")):
        raw = read_jdf(path)
        m = raw.metadata
        if raw.b is None:
            out["jdf"][path.name] = dict(sha256=m["source_sha256"], shape=raw.fid.shape, b_source=m["b_source"])
            continue
        p = m["parameters"]
        model = m["b_source"]
        shape = p["grad_shape"]["value"]
        small, large = (metadata_number(p,k) for k in ("delta","diffusion_time"))
        tau1 = metadata_number(p, "tau_1" if "tau_1" in p else "tau")
        tau2 = metadata_number(p,"tau_2") if "tau_2" in p else None
        rf = metadata_number(p,"x_pulse")
        args = (small,large,tau1,tau2)
        q = integrate_sequence(shape,*args)
        finite = integrate_sequence(shape,*args,rf)
        rb = (GAMMA*raw.gradients*small)**2*q["coefficient_s"]
        derived = (GAMMA*raw.gradients*small)**2*finite["coefficient_s"]
        phase = initial_dataset(raw).options["phase0_deg"]
        out["jdf"][path.name] = dict(
            sha256=m["source_sha256"], shape=raw.fid.shape, sequence=model["sequence"], pulse=model["pulse_program"],
            g_T_m=raw.gradients, b_app=raw.b, b_q=rb, b_q_finite_rf=derived, shape_name=shape,
            x_pulse_s=rf, tau1=tau1, tau2=tau2, delta=small, Delta=large, q=q, finite_rf=finite,
            st_coef_s=metadata_number(p,"st_coef"), app_gamma=model["gamma_rad_s_T"],
            relative_b_max=np.max(abs(rb/raw.b-1)), finite_relative_b_max=np.max(abs(derived/raw.b-1)),
            phase0_seed=phase, filter_samples=m["filter_delay_model"]["samples"])
        for wave in ("SQUARE","SINE","SQUARE_SINE","TRAPEZOID","S_RECTANGLE"):
            mod = deepcopy(p)
            mod.pop("st_coef")
            mod["grad_shape"]["value"] = wave
            predicted, info = factory_b(mod,m["context"],raw.gradients)
            integrated = integrate_sequence(wave,*args)
            name = model["sequence"]+"_"+wave
            out["waveforms"][name] = dict(independent=integrated,
                 app_coefficient_s=info["coefficient_s"], ratio=integrated["coefficient_s"]/info["coefficient_s"])
        if "Sacharoza" in path.name and tau2 is None:
            mutated = deepcopy(p)
            mutated["diffusion_time"]["value"] *= 2
            inconsistent_b, info = factory_b(mutated,m["context"],raw.gradients)
            actual_q = integrate_sequence(shape,small,2*large,tau1,None)
            corrupted_attenuation = np.exp(-raw.b*1.15**2*4*DUNIT)
            scaled_fit = fit_decay(raw.b, corrupted_attenuation)
            out["gradient_scale_15pct"] = dict(fit=mono_record(scaled_fit), true_D10=4.,
                                               predicted_D10=4*1.15**2)
            out["metadata_disagreement"] = dict(accepted=info["kind"]=="factory",
                nominal_relative_change=np.max(abs(inconsistent_b/raw.b-1)),
                D_inflation=actual_q["coefficient_s"]/q["coefficient_s"])
    return out


def individual_sweeps():
    data, reference = {}, {}
    for label, d, a in (("mono",[4],[1]),("separated",[2,16],[.6,.4]),("close",[4,4.5],[.5,.5])):
        y = decay(B,d,a)
        for points in (50,100,200,400):
            for spacing in ("log","linear"):
                for lo,hi in ((.25,40),(1.,25.)):
                    for lam in (1e-8,1e-5,1e-2):
                        p = ILTParameters(lo*DUNIT,hi*DUNIT,points=points,spacing=spacing,regularization=lam)
                        s = solve_ilt(B*1e9,y,p)
                        data[f"{label}_{points}_{spacing}_{lo:g}_{hi:g}_{lam:g}"] = ilt_record(s,d)
    rng = np.random.default_rng(SEED+2)
    for spacing in ("log","linear"):
        for order in (0,2):
            y = decay(B,[2,16],[.6,.4]) + rng.normal(0,.002,len(B))
            sigma = np.linspace(.001,.005,len(B))
            p = ILTParameters(.25*DUNIT,40*DUNIT,points=50,spacing=spacing,order=order)
            s = solve_ilt(B*1e9,y,p,sigma=sigma)
            d,k,design,target,ref = reference_ilt(B,y,p,sigma)
            gradient = 2*design.T@(design@s.weights-target)
            kkt = np.where(s.weights > 1e-12, abs(gradient), np.maximum(-gradient,0))
            reference[f"{spacing}_order{order}"] = dict(success=ref.success,
                grid_relerr=np.max(abs(s.diffusion/d-1)), weight_delta=np.max(abs(ref.x-s.weights)),
                prediction_delta=np.max(abs(k@(ref.x-s.weights)))*s.scale,
                objective_delta=np.sum((design@s.weights-target)**2)-np.sum((design@ref.x-target)**2),
                kkt_inf=np.max(kkt), nonnegative=bool(np.all(s.weights>=0)))
    noise = {}
    for sd in (0.,.001,.005,.02):
        for n in (12,24,49):
            for schedule in ("linear","geometric","short"):
                b = (np.linspace(0,35,n) if schedule=="linear" else np.r_[0,np.geomspace(.02,35,n-1)]
                     if schedule=="geometric" else np.linspace(0,3,n))
                for lam in (1e-8,1e-5,1e-2):
                    rows = []
                    for i in range(12):
                        truth = [2,16] if i%2 else [4]
                        amp = [.9,.1] if i%2 else [1.]
                        y = decay(b,truth,amp)+sd*np.random.default_rng(SEED+100+i).normal(size=n)
                        p = ILTParameters(.25*DUNIT,40*DUNIT,points=50,regularization=lam)
                        rows.append(dict(truth=truth, **ilt_record(solve_ilt(b*1e9,y,p),truth)))
                    noise[f"{sd:g}_{n}_{schedule}_{lam:g}"] = rows
    return data,reference,noise


def joint_sweeps():
    rng = np.random.default_rng(SEED+3)
    grid = ILTParameters(.25*DUNIT,40*DUNIT,points=100,order=0)
    out = {}
    def solve(label, truth, y, lam, mu, sigma=None, mask=None):
        params = GlobalParameters(replace(grid,regularization=lam),sparsity=mu)
        result = solve_global_ilt(B*1e9,y,params,sigma=sigma,include=mask)
        d = result.diffusion[list(result.components)]/DUNIT
        all_truth = sorted(set(sum(truth,[])))
        matched = match_peaks(d,all_truth)
        certificate = global_certificate(B,y,params,result,sigma,mask)
        out[label] = dict(truth=truth, peaks=d, fractions=result.fractions, rms=result.rms,
            rms_normalized=result.rms/result.scales, certificate=certificate, warnings=result.warnings,
            iterations=result.iterations, converged=result.converged, stationarity=result.stationarity,
            **matched)
        return result
    for ratio in (1,10,100,1000):
        y = np.column_stack([decay(B,[4],[1]),decay(B,[4],[1/ratio])])
        se = np.full(y.shape,.0015)
        y += rng.normal(0,se)
        solve(f"J_shared_ratio_{ratio}", [[4],[4]], y,1e-6,1e-4,se)
    clean = np.column_stack([decay(B,[2],[1]),decay(B,[2,16],[.5,.5]),decay(B,[16],[.005])])
    sigma = np.column_stack([np.full(len(B),.001),np.linspace(.001,.006,len(B)),np.full(len(B),.0005)])
    y = clean+rng.normal(0,sigma)
    mask = np.ones(y.shape,bool)
    mask[3::5,1] = False
    for lam in (1e-8,1e-6,1e-3):
        for mu in (0.,1e-4,.01,.1,10.):
            solve(f"K_partial_{lam:g}_{mu:g}", [[2],[2,16],[16]],y,lam,mu,sigma,mask)
    for m in (2,4,8):
        y = np.tile(decay(B,[2,16],[.92,.08])[:,None],(1,m))
        solve(f"replicated_{m}", [[2,16]]*m,y,1e-6,1e-4)
    y = np.column_stack([decay(B,[4],[1]),decay(B,[4.35],[.01])])
    for mu in (0.,1e-4,.01,.1):
        solve(f"close_distinct_{mu:g}", [[4],[4.35]],y,1e-6,mu)
    y = np.column_stack([decay(B,[4],[1.]),rng.normal(0,.001,len(B))])
    solve("noise_only_second_roi", [[4],[]],y,1e-6,1e-4,np.full(y.shape,.001))
    y = np.column_stack([decay(B,[3],[1.]),decay(B,[3,16],[.7,.3])])
    r = solve("separable_comparison", [[3],[3,16]],y,1e-5,0.)
    independent = np.column_stack([reference_ilt(B,v,replace(grid,regularization=1e-5))[4].x for v in y.T])
    independent_prediction = np.exp(-np.outer(B*1e9,r.diffusion)) @ independent * r.scales
    out["separable_comparison"]["max_prediction_diff_vs_BVLS"] = np.max(abs(r.predict(B*1e9)-independent_prediction))
    scaled = y*np.array([1.,1000.])
    control = solve_global_ilt(B*1e9,scaled,GlobalParameters(replace(grid,regularization=1e-5),sparsity=0.))
    out["separable_comparison"]["max_mass_change_after_rescaling"] = np.max(abs(r.weights-control.weights))
    out["separable_comparison"]["max_rescaled_prediction_difference"] = np.max(
        abs(r.predict(B*1e9)-control.predict(B*1e9)/np.array([1.,1000.])))
    loose = GlobalParameters(replace(grid,regularization=1e-8),max_iterations=20)
    limited = solve_global_ilt(B*1e9,y,loose)
    out["forced_20_iteration_limit"] = dict(iterations=limited.iterations, converged=limited.converged,
        stationarity=limited.stationarity, certificate=global_certificate(B,y,loose,limited))
    return out


def synthetic_raw(b9, components, *, phase_slope=0., baseline=0., noise=0.):
    """Causal complex FIDs with analytic exponential transverse decay."""
    sw, freq, offset, n = 8000., 600e6, 4.5, 8192
    t = np.arange(n)/sw
    out = np.zeros((len(b9),n),complex)
    for ppm, rate, amplitude, width, phase in components:
        hertz = (ppm-offset)*(freq/1e6)
        fid = np.exp(2j*np.pi*hertz*t-np.pi*width*t)
        theta = np.deg2rad(phase+phase_slope*(b9-b9.min())/np.ptp(b9))
        out += (decay(b9,[rate],[amplitude])*np.exp(1j*theta))[:,None]*fid
    # A perturbation confined to the first FID point becomes a constant real offset.
    out[:,0] += baseline*(1 + 1.5*b9/np.max(b9))
    if noise:
        rng = np.random.default_rng(SEED+4)
        out += noise*(rng.normal(size=out.shape)+1j*rng.normal(size=out.shape))
    return RawData(out,np.full(len(b9),np.nan),b9*1e9,
                   dict(x_sweep_Hz=sw,x_frequency_Hz=freq,offset_ppm=offset,
                        source_name="analytic_FID",digital_filter=False))


def processing_cases():
    b = np.linspace(0,14,24)
    core = [(5.25,3.8,1,3.,0)]
    roi = (5.18,5.32)
    fixtures = [
        ("clean", core, {}, {}),
        ("common_25deg", [(5.25,3.8,1,3.,25)], {}, {}),
        ("gradient_phase_55deg",core,dict(phase_slope=55),{}),
        ("uncorrected_baseline",core,dict(baseline=8.),{}),
        ("clean_baseline_mask",core,dict(baseline=8.),
             dict(baseline_enabled=True,baseline_regions=[(0,1),(9,10)])),
        ("linear_baseline_solvent_tail",core+[(4.94,19,60,9,12)],dict(baseline=8.),
             dict(baseline_enabled=True,baseline_regions=[(4.98,5.08),(5.55,5.8)])),
        ("unresolved_two_peaks",core+[(5.258,18.8,.35,4,0)],{},{}),
        ("solvent_tail",core+[(5.01,19,90,12,14)],{},{}),
        ("negative_analyte",[(5.25,3.8,1,3.,180)],{},{}),
        ("low_SNR",[(5.25,3.8,1e-3,3.,0)],dict(noise=.04),{}),
    ]
    out = {}
    for key, components, effects, extra in fixtures:
        raw = synthetic_raw(b,components,**effects)
        options = default_options(raw)
        options.update(extra)
        entry = dict(id="A",bounds=roi)
        d = process_fids(raw,options,[roi])
        row = analyze(d,[entry])[0]
        integrals = row["result"]["corrected_integral"]
        peak = d.ppm[np.argmax(d.spectrum[0].real)]
        out[key] = dict(b9=b, y=integrals, fit=mono_record(row["fit"]), raw_readonly=not raw.fid.flags.writeable,
                        prediction=row["fit"].get("prediction"), max_peak_ppm=peak,
                        negative_count=np.count_nonzero(integrals<=0), process_warnings=d.diagnostics["warnings"])
    raw = synthetic_raw(b,core)
    raw_freq = 600.
    n = raw.fid.shape[-1]
    # Direct evaluation of the DFT at five signed integer bins, independent of FFT.
    bins = np.array([-1020,-1,0,13,461])
    fft_reference = raw.fid @ np.exp(-2j*np.pi*np.arange(n)[:,None]*bins[None,:]/n)
    no_change = process_fids(raw)
    freq_per_bin = raw.metadata["x_sweep_Hz"]/n
    target = 4.5+bins*freq_per_bin/raw_freq
    selected = np.searchsorted(no_change.ppm,target-1e-12)
    out["DFT_and_ppm"] = dict(error=np.max(abs(no_change.phased[:,selected]-fft_reference)),
                              ppm_error=np.max(abs(no_change.ppm[selected]-target)))
    out["common_linear_processing"] = {}
    for zf in (1,2,4,8):
        d = process_fids(raw,dict(default_options(raw),lb_Hz=3.,gaussian_FWHM_Hz=4.,
                                  zero_fill_factor=zf,phase0_deg=31.,phase1_deg=-80.,
                                  pivot_ppm=5.25,ppm_shift=.11))
        y = signed_integrals(d.ppm,d.spectrum.real,np.array(roi)+.11)
        out["common_linear_processing"][str(zf)] = dict(fit=mono_record(fit_decay(b*1e9,y)),
            amplitude0=y[0], attenuation_error=np.max(abs(y/y[0]-np.exp(-b*.38))))
    out["bad_fit_inputs"] = {}
    for name, xx, yy in (("nan",b,np.r_[np.nan,np.ones(len(b)-1)]),
                         ("infinite_b",np.r_[np.inf,b[1:]],np.ones(len(b))),
                         ("only_two_distinct_b",np.resize([0.,1.],len(b)),np.ones(len(b))),
                         ("zero_signal",b,np.zeros(len(b))),
                         ("constant_signal",b,np.ones(len(b)))):
        out["bad_fit_inputs"][name] = mono_record(fit_decay(xx*1e9,yy))
    return out


def local_fixture(specifications):
    ppm = np.linspace(.25,len(specifications)+.75,round((len(specifications)+.5)*240)+1)
    spectra = np.zeros((len(B),len(ppm)))
    entries = []
    rng = np.random.default_rng(SEED+5)
    for i,(rates, amps, sigma, offsets) in enumerate(specifications,1):
        bounds = (i-.31, i+.31)
        mask = (ppm>=bounds[0])&(ppm<=bounds[1])
        p = ppm[mask]
        # Each profile has a known shape in ppm and its own temporal evolution.
        shapes = np.exp(-.5*((p[None,:] - i-np.asarray(offsets)[:,None])/.036)**2)
        values = np.exp(-np.outer(B,np.asarray(rates)/10)) @ (np.asarray(amps)[:,None]*shapes)
        values += sigma*rng.normal(size=values.shape)
        spectra[:,mask] = values
        y = signed_integrals(ppm,spectra,bounds)
        entries.append(dict(id=f"R{i}", bounds=bounds, fit=fit_decay(B*1e9,y),
                            result=dict(corrected_integral=y)))
    dataset = SimpleNamespace(b=B*1e9,ppm=ppm,spectrum=spectra,options={"reference_layer":0})
    return dataset, entries


def adaptive_maps():
    d, rows = local_fixture([([3.],[1.],.003,[0]),
        ([3.],[.008],.0009,[0]),([3.25],[1],.016,[0]),([3.25],[.6],.012,[0]),
        ([9.],[.45],.004,[0]),([9.],[.2],.005,[0]),
        ([3.,9.],[.22,.15],0.,[-.032,.032]),([4.],[0.],.001,[0])])
    out = {"truth":[[3],[3],[3.25],[3.25],[9],[9],[3,9],[]], "native":{
        row["id"]:mono_record(row["fit"]) for row in rows}, "sweep":{}}
    for tol in (1.,3.,5.,8.,15.,30.):
        for balance in (1.,10.,100.,1000.):
            for lam in (1e-6,1e-4):
                res = adaptive_local_ilt(d,rows,tol,balance,lam)
                groups = [{"ids":g.ids, "D10":g.dominant/DUNIT if g.dominant is not None else None,
                    "converged":g.converged, "rms":g.rms,"amp_ratio":g.amplitude_ratio,
                    "rejected":g.rejected} for g in res.groups]
                profiles = {}
                for row in rows:
                    if row["id"] in res.profiles:
                        seq = res.profiles[row["id"]]
                        i,j = np.searchsorted(d.ppm,row["bounds"])
                        real = np.maximum(d.spectrum[0,i:j],0)
                        total = sum(v for _,v in seq)
                        area = trapezoid(total,x=d.ppm[i:j])
                        original = trapezoid(real,x=d.ppm[i:j])
                        profiles[row["id"]] = dict(D10=[g/DUNIT for g,_ in seq],
                            profile_max_error=np.max(abs(total-real)), area_ratio=area/original)
                out["sweep"][f"{tol:g}_{balance:g}_{lam:g}"] = dict(groups=groups, profiles=profiles,
                                                                  unresolved=res.unresolved)
    tol = 15.
    out["display_only"] = [
        dict(ids=g.ids, representative_D10=g.representative/DUNIT)
        for g in alignment_plan(d.b,rows,tol,(),False)]
    out["shared_fit"] = [
        dict(ids=g.ids, representative_D10=g.representative/DUNIT,
             shared_D10=g.shared/DUNIT if g.shared else None, incompatible=g.incompatible,
             diagnostics=g.diagnostics) for g in alignment_plan(d.b,rows,tol,(),True)]
    # Isolate the renderer with an oracle ILT result, including a real but subthreshold component.
    dm, rm = local_fixture([([3.,12.],[1.,.03],0.,[-.05,.05])])
    solution = SimpleNamespace(converged=True, diffusion=np.array([1.,3.,6.,12.,20.])*DUNIT,
        weights=np.array([0.,1.,0.,.03,0.])[:,None], intervals=((0,3),(3,5)),
        components=(1,3), scales=np.array([1.]), rms=np.array([0.]))
    profiles, rejected = ilt_profiles(dm.b,dm.ppm,dm.spectrum,0,rm,("R1",),solution,
                                      conserve_reference=True)
    i,j = np.searchsorted(dm.ppm,rm[0]["bounds"])
    trace = profiles["R1"][0][1]
    out["weak_component_oracle_map"] = dict(truth_D10=[3,12], oracle_weak_fraction=.03/1.03,
        plotted_D10=[v/DUNIT for v,_ in profiles["R1"]], rejection=rejected,
        entire_reference_at_strong_D=bool(np.array_equal(trace,dm.spectrum[0,i:j])))
    first = dict(rm[0],id="left",bounds=(.8,1.03))
    second = dict(rm[0],id="right",bounds=(.98,1.2))
    bands, overlap = roi_bands(dm.ppm,dm.spectrum[0],[first,second])
    common = (dm.ppm>=.98)&(dm.ppm<1.03)
    dm2, rm2 = local_fixture([([3.,12.],[1.,.55],0.,[-.024,.024])])
    solution2 = SimpleNamespace(converged=True, diffusion=np.array([1.,3.,6.,12.,20.])*DUNIT,
        weights=np.array([0.,1.,0.,.55,0.])[:,None], intervals=((0,3),(3,5)),
        components=(1,3), scales=np.array([1.]), rms=np.array([0.]))
    reconstructed, not_drawn = ilt_profiles(dm2.b,dm2.ppm,dm2.spectrum,0,rm2,("R1",),solution2,
                                            conserve_reference=True)
    a,z = np.searchsorted(dm2.ppm,rm2[0]["bounds"])
    ideal = np.array([amplitude*np.exp(-.5*((dm2.ppm[a:z]-1-offset)/.036)**2)
                      for amplitude,offset in ((1.,-.024),(.55,.024))])
    actual = np.array([v for _,v in reconstructed["R1"]])
    out["resolved_oracle_map"] = dict(plotted_D10=[rate/DUNIT for rate,_ in reconstructed["R1"]],
         rejected=not_drawn,max_profile_error=np.max(abs(actual-ideal)),
         max_reference_error=np.max(abs(actual.sum(axis=0)-dm2.spectrum[0,a:z])))
    out["overlapping_bounds"] = dict(overlap=overlap, positive_original_samples=np.count_nonzero(dm.spectrum[0,common]),
        nonzero_common_samples=sum(np.count_nonzero(b.intensity[(b.ppm>=.98)&(b.ppm<1.03)]) for b in bands))
    out["map_settings_example"] = dict(max_D10=19,weak_relative_height=.001,floor_pct=1,
        default_fwhm_D10=.03*19, user_min_fwhm_D10=.002*19,user_max_fwhm_D10=.15*19)
    return out


def real_processing_roundtrip():
    summaries = {}
    for path in sorted((ROOT/"examples").rglob("*.jdf")):
        raw = read_jdf(path)
        if raw.b is None:
            continue
        d = initial_dataset(raw)
        bounds = [(5.24,5.34),(4.06,4.14),(3.52,3.60)] if "Sacharoza" in path.name else [(8.00,8.24),(5.78,6.03)]
        entries = [dict(id=f"R{i+1}",name=f"R{i+1}",group="audit",bounds=roi,fit_mask=[True]*len(raw.b),
                        sigma=[None]*len(raw.b),weighted=False) for i,roi in enumerate(bounds)]
        results = analyze(d,entries)
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            save_project(folder/"repeat.proj",d,entries)
            reload, opts, rois = read_project(folder/"repeat.proj")
            repeated = analyze(process_fids(reload,opts,bounds),rois)
            reconst = [
                np.max(abs(old["result"]["corrected_integral"]-new["result"]["corrected_integral"]))
                for old,new in zip(results,repeated)
            ]
            save_integrals(folder/"repeat.ints",entries)
            templates = read_integrals(folder/"repeat.ints")
            export_bundle(folder/"bundle",d,results)
            archive = np.load(folder/"bundle"/"complex_data.npz")
            summary = dict(ppm_shape=len(d.ppm), max_roundtrip_integral_error=max(reconst),
                    integral_template_preserved=all(tuple(t["bounds"])==tuple(e["bounds"])
                        and t["id"]==e["id"] and t["group"]==e["group"] for t,e in zip(templates,entries)),
                    npz_keys=archive.files, exports=sorted(f.name for f in (folder/"bundle").iterdir()),
                    fits=[dict(bounds=r["bounds"],fit=mono_record(r["fit"])) for r in results])
            summary["max_npz_difference"] = max(float(np.max(abs(archive[key]-value))) for key,value in
                (("fid", raw.fid), ("b", raw.b), ("ppm", d.ppm), ("spectrum", d.spectrum), ("phased", d.phased)))
            archive.close()
            summaries[path.name] = summary
    return summaries


def plain(value):
    if isinstance(value,np.ndarray):
        return plain(value.tolist())
    if isinstance(value,np.generic):
        return plain(value.item())
    if isinstance(value,dict):
        return {k:plain(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [plain(v) for v in value]
    if isinstance(value,float) and not np.isfinite(value):
        return None
    return value


def main():
    import openpyxl
    import PySide6
    import pyqtgraph
    source = subprocess.check_output(["git","rev-parse","v2.0.0^{commit}"],cwd=ROOT,text=True).strip()
    dirty = subprocess.check_output(["git","diff",source,"--name-only","--","bin"],cwd=ROOT,text=True).strip()
    if dirty:
        raise RuntimeError(f"Released source differs from v2.0.0: {dirty}")
    sections = [
        ("physical_b",physical_b), ("cases",benchmark_cases), ("monte_carlo",monte_carlo),
        ("individual",individual_sweeps), ("joint",joint_sweeps),
        ("processing",processing_cases), ("adaptive_maps",adaptive_maps), ("real_roundtrip",real_processing_roundtrip),
    ]
    config = np.show_config(mode="dicts")["Build Dependencies"]
    results = dict(source=source, tag="v2.0.0", python=sys.version, platform=platform.platform(),
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,pyside=PySide6.__version__,
                      pyqtgraph=pyqtgraph.__version__,openpyxl=openpyxl.__version__,
                      blas=config["blas"]["version"]), seed=SEED, log_match_tol=TOL)
    OUT.mkdir(exist_ok=True)
    for name, function in sections:
        print(name,flush=True)
        results[name] = function()
        (OUT/"benchmark.json").write_text(json.dumps(plain(results),ensure_ascii=False,indent=2,allow_nan=False)+"\n",
                                           encoding="utf-8")
    print("Complete:", OUT/"benchmark.json",flush=True)


if __name__ == "__main__":
    main()
