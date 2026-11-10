"""Diffusion grouping and optional spectral profiles for the experimental map."""
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
from scipy.optimize import minimize_scalar, nnls
from scipy.integrate import trapezoid

from d_statistics import finite
from global_ilt import component_intervals
from ilt import D_UNIT, _kernel


MAX_PROFILE_RELATIVE_RESIDUAL = .15


@dataclass(frozen=True)
class DiffusionGroup:
    ids: tuple[str, ...]
    representative: float
    shared: float | None = None
    # Each diagnostic is (independent RMS, shared RMS), in normalized intensity units.
    diagnostics: tuple[tuple[float, float], ...] = ()
    incompatible: tuple[str, ...] = ()


def eligible_diffusion(row):
    fit = row.get('fit') or {}
    d, r2, se = [finite(fit.get(key)) for key in ('D_m2_s', 'R2', 'D_standard_error_m2_s')]
    if (not fit.get('ok') or row.get('error') or d is None or d <= 0 or r2 is None or not .95 <= r2 <= 1
            or fit.get('n_active', 0) < 5 or (se is not None and se > .25*d)):
        return None
    return d


def cluster_diffusion(rows, tolerance, excluded=()):
    """Deterministic ordered complete-link grouping for a 1D relative distance."""
    candidates = [(d, row['id'], row) for row in rows
                  if row['id'] not in excluded and (d := eligible_diffusion(row)) is not None]
    candidates.sort(key=lambda item: item[:2])
    groups = []
    current = []
    for candidate in candidates:
        d = candidate[0]
        if current and 200*(1-current[0][0]/d)/(1+current[0][0]/d) > tolerance:
            groups.append(current)
            current = []
        current.append(candidate)
    groups.append(current)
    output = []
    for group in groups:
        if len(group) < 2:
            continue
        d = np.array([entry[0] for entry in group])
        se = np.array([entry[2]['fit'].get('D_standard_error_m2_s') for entry in group], dtype=float)
        representative = np.average(d, weights=(se.min()/se)**2) if np.isfinite(se).all() and np.all(se > 0) else np.median(d)
        output.append(DiffusionGroup(tuple(sorted(entry[1] for entry in group)), float(representative)))
    return tuple(output)


def _attenuation_column(b, row):
    fit = row['fit']
    mask = np.asarray(fit['fit_mask'], dtype=bool)
    x = (b[mask]-b[mask].min())*D_UNIT
    y = np.asarray(row['result']['corrected_integral'])[mask]
    scale = np.max(np.abs(y))
    weight = np.ones(len(y))
    if fit.get('weighted'):
        sigma = np.asarray(row['sigma'], dtype=float)[mask]
        weight = sigma.min()/sigma
    weight /= np.linalg.norm(weight)
    return x, y/scale*weight, weight, np.linalg.norm(np.asarray(fit['residuals'])[mask]/scale*weight)


def fit_shared_group(b, group, by_id):
    """Profile out independent amplitudes and fit one shared positive D."""
    columns = [_attenuation_column(b, by_id[sid]) for sid in group.ids]

    def residuals(rate):
        residual = []
        for x, y, weight, _ in columns:
            model = np.exp(-x*rate)*weight
            amplitude = max(0., np.dot(model, y)/np.dot(model, model))
            residual.append(np.linalg.norm(y-model*amplitude))
        return np.array(residual)

    center = group.representative/D_UNIT
    result = minimize_scalar(lambda rate: np.linalg.norm(residuals(rate))**2,
                             bounds=(center/2, center*2), method='bounded',
                             options={'xatol': center*1e-8})
    if not result.success or not center*.5001 < result.x < center*1.999:
        return DiffusionGroup(group.ids, group.representative, incompatible=group.ids)
    diagnostics = tuple((float(col[3]), float(rms)) for col, rms in zip(columns, residuals(result.x)))
    incompatible = tuple(sid for sid, (original, shared) in zip(group.ids, diagnostics)
                         if shared > max(1.5*original, original+.005))
    return DiffusionGroup(group.ids, group.representative, float(result.x*D_UNIT), diagnostics, incompatible)


def alignment_plan(b, rows, tolerance, excluded, shared=False):
    groups = cluster_diffusion(rows, tolerance, excluded)
    if shared and b is not None:
        by_id = {row['id']: row for row in rows}
        groups = tuple(fit_shared_group(b, group, by_id) for group in groups)
    return groups


def profile_residual(target, predicted, noise_energy=0.):
    energy = np.sum(np.abs(target)**2)
    if not energy:
        return dict(relative_residual=None, excess_residual=None, relative_noise=None)
    error = np.sum(np.abs(target-predicted)**2)
    return dict(relative_residual=float(np.sqrt(error/energy)),
                excess_residual=float(np.sqrt(max(0.,error-noise_energy)/energy)),
                relative_noise=float(np.sqrt(noise_energy/energy)))


def peak_attenuation_rms(b, row, diffusion, scale, fallback):
    """Validate the peak used as the map coordinate, not its wider ILT interval."""
    fit = row['fit']
    y = row.get('result',{}).get('corrected_integral')
    if y is None:
        return None, fallback
    mask = np.asarray(fit['fit_mask'], dtype=bool)
    x = b[mask]-b[mask].min()
    values = np.asarray(y)[mask]/scale
    weights = np.ones(len(x))
    if fit.get('weighted'):
        sigma = np.asarray(row['sigma'],dtype=float)[mask]
        weights = sigma.min()/sigma
    weights /= np.linalg.norm(weights)
    values *= weights
    basis = np.exp(-x*diffusion)*weights
    amplitude = max(0.,np.dot(basis,values)/np.dot(basis,basis))
    native = np.linalg.norm(np.asarray(fit['residuals'])[mask]/scale*weights)
    return float(np.linalg.norm(values-amplitude*basis)), float(native)


def ilt_profiles(b, ppm, spectra, reference_layer, rows, ids, solution, *, conserve_reference=False,
                 withheld=None, diagnostics=None, noise_sigma=None):
    """Estimate spectral profiles and retain a per-ROI account of assignment."""
    profiles, rejected = {}, []
    withheld = {} if withheld is None else withheld
    diagnostics = {} if diagnostics is None else diagnostics
    if not solution.converged:
        diagnostics.update({sid:dict(reason='solver_not_converged', optimized=True) for sid in ids})
        return profiles, tuple(ids)
    by_id = {row['id']: row for row in rows}
    kernel = _kernel(b, solution.diffusion)
    for column, sid in enumerate(ids):
        if sid not in by_id:
            continue
        entry = by_id[sid]
        mask = np.asarray(entry['fit']['fit_mask'], dtype=bool)
        start, stop = np.searchsorted(ppm, entry['bounds'])
        prediction = np.array([kernel[:, lo:hi] @ solution.weights[lo:hi,column]
                               for lo,hi in solution.intervals]).reshape(-1,len(b)).T
        masses = prediction[reference_layer]
        keep = (masses > 0) & (masses >= .05*masses.sum())
        peaks = [peak for peak, enabled in zip(solution.components, keep) if enabled]
        omitted = float(masses[~keep].sum()/masses.sum()) if masses.sum() > 0 else 0.
        record = dict(optimized=True, omitted_fraction=omitted, component_min_fraction=.05,
                      candidate_D=tuple(float(solution.diffusion[k]) for k in peaks),
                      spectral_limit=MAX_PROFILE_RELATIVE_RESIDUAL)
        if omitted:
            record['omitted_components'] = tuple(
                (float(solution.diffusion[peak]), float(mass/masses.sum()))
                for peak,mass,retained in zip(solution.components,masses,keep) if mass>0 and not retained)
        diagnostics[sid] = record

        def reject(reason):
            record['reason'] = reason
            record['unexplained_positive_fraction'] = (1. if np.any(spectra[reference_layer,start:stop]>0) else None)
            rejected.append(sid)
            withheld[sid] = record

        if not peaks or stop-start < 2:
            reject('no_component' if not peaks else 'insufficient_profile')
            continue
        if len(peaks) == 1 and omitted > 0:
            reject('discarded_component')
            continue
        if any(peak in (0,len(solution.diffusion)-1) for peak in peaks):
            reject('boundary_peak')
            continue
        scale = solution.scales[column]
        original_rms = entry['fit']['RMSE']/scale
        integral_limit = max(1.5*original_rms,original_rms+.005)
        record.update(integral_rms=float(solution.rms[column]/scale), integral_limit=integral_limit)
        if record['integral_rms'] > integral_limit:
            reject('integral_mismatch')
            continue
        basis = prediction[mask][:,keep]/masses[keep]
        target = spectra[mask,start:stop]
        trace = spectra[reference_layer,start:stop]
        if noise_sigma is not None:
            sigma = np.asarray(noise_sigma)
            variance = sigma[mask]**2
            ref_noise = sigma[reference_layer]*np.sqrt(len(trace))
            if ref_noise > 0:
                record['reference_snr_estimate'] = float(np.linalg.norm(trace)/ref_noise)
                if record['reference_snr_estimate'] < 5:
                    reject('insufficient_spectral_snr')
                    continue
        else:
            variance = np.zeros(mask.sum())
        if len(peaks) == 1:
            peak = solution.diffusion[peaks[0]]
            rms, native = peak_attenuation_rms(b,entry,peak,scale,original_rms)
            record.update(peak_rms=rms, peak_limit=max(1.5*native,native+.005))
            if rms is not None and rms > record['peak_limit']:
                reject('peak_mismatch')
                continue
            concerns = entry['fit'].get('measurement_findings',())
            if any(finding['severity']=='warning' and finding['kind'] in ('shape','phase')
                   for finding in concerns):
                reject('profile_dependent_diffusion')
                continue
            if noise_sigma is not None:
                variance += basis[:,0]**2*sigma[reference_layer]**2
                # Reference residual is exactly zero when it is included.
                variance[np.flatnonzero(mask)==reference_layer] = 0.
            record.update(profile_residual(target,basis*trace,len(trace)*variance.sum()))
            if record['excess_residual'] is None or record['excess_residual'] > MAX_PROFILE_RELATIVE_RESIDUAL:
                reject('spectral_mismatch')
                continue
            coeff = trace[None,:]
        else:
            # Fit subthreshold positive components as nuisance terms. Removing
            # their bases before NNLS would leak their signal into retained peaks.
            modeled = masses > 0
            basis = prediction[mask][:,modeled]/masses[modeled]
            target_scale = np.max(np.abs(target))
            condition = np.linalg.cond(basis)
            record['basis_condition'] = float(condition)
            if mask.sum() < 2*np.count_nonzero(modeled)+2 or condition > 100 or not target_scale:
                reject('ill_conditioned_basis')
                continue
            try:
                coeff = np.array([nnls(basis, signal)[0] for signal in (target/target_scale).T]).T*target_scale
            except (RuntimeError, ValueError, np.linalg.LinAlgError):
                reject('spectral_solver_failed')
                continue

            def assess(values):
                record.update(profile_residual(target,basis@values,len(trace)*variance.sum()))
                areas = trapezoid(values,x=ppm[start:stop],axis=1)
                if areas.sum() <= 0:
                    return False
                shares = np.max(abs(areas/areas.sum()-masses[modeled]/masses[modeled].sum()))
                record['share_discrepancy'] = float(shares)
                return (record['excess_residual'] is not None
                        and record['excess_residual'] <= MAX_PROFILE_RELATIVE_RESIDUAL and shares <= .2)

            if not assess(coeff):
                reject('spectral_mismatch')
                continue
            reference = np.maximum(trace,0.)
            if conserve_reference and omitted == 0:
                total = coeff.sum(axis=0)
                if not reference.any() or reference[total == 0].sum() > .01*reference.sum():
                    reject('unexplained_reference')
                    continue
                adjusted = coeff*np.divide(reference,total,out=np.zeros_like(total),where=total>0)
                if not assess(adjusted):
                    reject('spectral_mismatch')
                    continue
                coeff = adjusted
            coeff = coeff[keep[modeled]]
        positive = np.maximum(trace,0.)
        represented = np.maximum(coeff,0.).sum(axis=0)
        area = trapezoid(positive,x=ppm[start:stop])
        remaining = trapezoid(np.maximum(positive-represented,0.),x=ppm[start:stop])
        record.update(reason='assigned', unexplained_positive_fraction=float(remaining/area) if area>0 else None)
        profiles[sid] = tuple((solution.diffusion[peak], values) for peak,values in zip(peaks,coeff))
    return profiles, tuple(rejected)


def individual_map_solution(result):
    """Adapt the solved individual distribution without changing its maxima."""
    peaks = tuple(sorted(result.peaks))
    return SimpleNamespace(converged=True, diffusion=result.diffusion, weights=result.weights[:,None],
                           components=peaks, intervals=component_intervals(result.weights, peaks),
                           scales=np.array([result.scale]), rms=np.array([result.rms]))


def exploratory_profiles(b, ppm, spectra, reference_layer, rows, ids, solution):
    """Project all positive intervals of a converged solution onto measured spectral bins."""
    profiles = {}
    diagnostics = {sid: dict(reason='solver_not_converged', optimized=True, solved=False) for sid in ids}
    if not solution.converged:
        return profiles, diagnostics
    kernel = _kernel(b, solution.diffusion)
    columns = dict(zip(ids, range(len(ids))))
    for row in rows:
        sid = row['id']
        if sid not in columns:
            continue
        column = columns[sid]
        record = diagnostics[sid] = dict(reason='no_component', optimized=True, solved=True)
        weights = solution.weights[:,column]
        prediction = np.column_stack([kernel[:,lo:hi] @ weights[lo:hi] for lo,hi in solution.intervals]
                                     or [np.zeros(len(b))])
        ref_masses = prediction[reference_layer]
        keep = ref_masses > 0
        if not keep.any():
            continue
        masses = ref_masses[keep]*solution.scales[column]
        intervals = [bounds for bounds, active in zip(solution.intervals,keep) if active]
        peaks = [lo + np.argmax(weights[lo:hi]) for lo,hi in intervals]
        rates = solution.diffusion[peaks]
        record.update(candidate_D=tuple(map(float,rates)), model_reference_area=tuple(map(float,masses)),
                      boundary_peak=any(k in (0,len(weights)-1) for k in peaks))
        start, stop = np.searchsorted(ppm, row['bounds'])
        mask = np.asarray(row['fit']['fit_mask'], dtype=bool)
        basis = prediction[mask][:,keep]/ref_masses[keep]
        target = spectra[mask,start:stop]
        trace = np.maximum(spectra[reference_layer,start:stop],0.)
        target_scale = np.max(abs(target),initial=0.)
        if stop-start < 2 or not target_scale or not trace.any():
            record['reason'] = 'insufficient_profile'
            continue
        try:
            singular = np.linalg.svd(basis,compute_uv=False)
            tolerance = singular[0]*max(basis.shape)*np.finfo(basis.dtype).eps
            if np.count_nonzero(singular > tolerance) < len(peaks):
                record['reason'] = 'ill_conditioned_basis'
                continue
            record['basis_condition'] = float(singular[0]/singular[-1])
            if len(peaks) == 1:
                coeff = np.maximum(0., (basis.T @ target)/np.sum(basis**2))
            else:
                coeff = np.array([nnls(basis, signal)[0] for signal in (target/target_scale).T]).T*target_scale
        except (RuntimeError, ValueError, np.linalg.LinAlgError):
            record['reason'] = 'spectral_solver_failed'
            continue
        axis = ppm[start:stop]
        area = trapezoid(coeff,x=axis,axis=1)
        # Cap component area by ILT's reference prediction and summed bin intensity
        # by the positive measured reference. Never scale amplitudes upward.
        coeff *= np.divide(np.minimum(area,masses),area,out=np.zeros_like(area),where=area>0)[:,None]
        total = coeff.sum(axis=0)
        coeff *= np.divide(np.minimum(total,trace),total,out=np.zeros_like(total),where=total>0)
        area = trapezoid(coeff,x=axis,axis=1)
        record.update(profile_residual(target,basis@coeff),
                      mapped_reference_area=tuple(map(float,area)),
                      unexplained_positive_fraction=float(np.clip(1-area.sum()/trapezoid(trace,x=axis),0.,1.)))
        active = area > 0
        record['reason'] = 'assigned' if active.any() else 'no_positive_reconstructed_profile'
        if active.any():
            profiles[sid] = tuple(zip(rates[active],coeff[active]))
    return profiles, diagnostics
