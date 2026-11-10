"""Diffusion grouping and optional spectral profiles for the experimental map."""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar, nnls
from scipy.integrate import trapezoid

from d_statistics import finite
from ilt import D_UNIT, _kernel


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


def ilt_profiles(b, ppm, spectra, reference_layer, rows, ids, solution):
    """Estimate nonnegative ppm profiles using component decay bases fixed by ILT."""
    profiles, rejected = {}, []
    if not solution.converged:
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
        scale = solution.scales[column]
        original_rms = entry['fit']['RMSE']/scale
        if (not peaks or stop-start < 2 or any(peak in (0, len(solution.diffusion)-1) for peak in peaks)
                or solution.rms[column]/scale > max(1.5*original_rms, original_rms+.005)):
            rejected.append(sid)
            continue
        if len(peaks) == 1:
            profiles[sid] = ((solution.diffusion[peaks[0]], spectra[reference_layer,start:stop]),)
            continue
        basis = prediction[mask][:,keep]/masses[keep]
        target = spectra[mask,start:stop]
        target_scale = np.max(np.abs(target))
        if (mask.sum() < 2*len(peaks)+2 or np.linalg.cond(basis) > 100 or not target_scale):
            rejected.append(sid)
            continue
        target = target/target_scale
        try:
            coeff = np.array([nnls(basis, signal)[0] for signal in target.T]).T
        except (RuntimeError, ValueError, np.linalg.LinAlgError):
            rejected.append(sid)
            continue
        weights = trapezoid(coeff, x=ppm[start:stop], axis=1)
        if (np.linalg.norm(basis@coeff-target)/np.linalg.norm(target) > .15 or weights.sum() <= 0
                or np.max(np.abs(weights/weights.sum()-masses[keep]/masses[keep].sum())) > .2):
            rejected.append(sid)
            continue
        profiles[sid] = tuple((solution.diffusion[peak], trace*target_scale) for peak, trace in zip(peaks, coeff))
    return profiles, tuple(rejected)
