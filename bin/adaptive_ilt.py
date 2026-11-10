"""Local, quality-screened joint inversion of DOSY ROI groups."""
from dataclasses import dataclass, field, replace
from types import SimpleNamespace

import numpy as np

from dosy_alignment import cluster_diffusion, fit_shared_group, ilt_profiles
from global_ilt import GlobalParameters, SolveCancelled, identify_components, solve_global_ilt
from ilt import ILTParameters, MIN_DISTINCT_B, solve_ilt


@dataclass(frozen=True)
class LocalGroup:
    ids: tuple[str, ...]
    d_range: tuple[float, float]
    dominant: float | None
    amplitude_ratio: float
    rms: float | None
    converged: bool
    rejected: tuple[str, ...]
    independent: bool = False


@dataclass(frozen=True)
class LocalResult:
    groups: tuple[LocalGroup, ...]
    profiles: dict
    unresolved: tuple[str, ...]
    membership: dict  # group IDs -> (ROI, D) components matched from a separate mixed ROI
    withheld: dict = field(default_factory=dict)


def amplitude_weights(amplitude, snr, limit):
    """Protect against relatively noisy weak traces without grouping by amplitude."""
    strong = amplitude >= amplitude.max()/limit
    weights = np.ones(len(amplitude))
    weights[~strong] = np.clip(snr[~strong]/np.median(snr[strong]), .25, 1.)
    return weights


def reference_quality(row, reference):
    fit = row.get('fit') or {}
    if not fit.get('ok') or row.get('error') or fit['D_m2_s'] <= 0:
        return None
    amplitude = float(row['result']['corrected_integral'][reference])
    if not np.isfinite(amplitude) or amplitude <= 0:
        return None
    sigma = np.asarray(row.get('sigma'), dtype=float)
    noise = fit['RMSE']
    if sigma.shape == np.shape(fit['fit_mask']) and np.isfinite(sigma[reference]) and sigma[reference] > 0:
        noise = sigma[reference]
    snr = amplitude/max(noise, amplitude*1e-12)
    return (amplitude, snr) if snr >= 5 else None


def compatible_groups(b, rows, tolerance):
    """Remove incompatible members, then look for groups among removed members."""
    by_id = {row['id']: row for row in rows}
    remaining = list(rows)
    while len(remaining) >= 2:
        candidates = cluster_diffusion(remaining, tolerance)
        if not candidates:
            break
        peeled = []
        for group in candidates:
            group = fit_shared_group(b, group, by_id)
            if not group.incompatible:
                yield group
                continue
            if not group.diagnostics:
                continue
            allowed = set(group.ids)-set(group.incompatible)
            if not allowed:
                # Retry the rest without the member with the worst standardized deterioration.
                worst = max(zip(group.ids, group.diagnostics), key=lambda pair:
                            pair[1][1]/max(1.5*pair[1][0], pair[1][0]+.005))[0]
                allowed = set(group.ids)-{worst}
            if len(allowed) >= 2:
                yield from compatible_groups(b, [by_id[sid] for sid in sorted(allowed)], tolerance)
            peeled.extend(by_id[sid] for sid in group.ids if sid not in allowed)
        remaining = peeled


def fit_inputs(rows):
    y = np.column_stack([row['result']['corrected_integral'] for row in rows])
    mask = np.column_stack([row['fit']['fit_mask'] for row in rows]).astype(bool)
    sigma = np.ones_like(y)
    weighted = False
    for i, row in enumerate(rows):
        if 'sigma' in row:
            values = np.asarray(row['sigma'], dtype=float)
            if values.shape != mask[:,i].shape:
                continue
            active = values[mask[:,i]]
            if np.isfinite(active).all() and np.all(active > 0):
                sigma[:,i] = values
                weighted = True
    return y, mask, sigma if weighted else None


def grid_for(rows, regularization):
    d = np.array([row['fit']['D_m2_s'] for row in rows])
    return ILTParameters(d.min()/20, d.max()*20, regularization=regularization, order=0)


def separate_mixed_roi(dataset, row, regularization):
    """Only retain a split when two penalties support similar interior maxima."""
    y, mask, sigma = fit_inputs([row])
    y, mask = y[:,0], mask[:,0]
    sigma = sigma[:,0] if sigma is not None else None
    grid = grid_for([row], regularization)
    probe = 2*regularization if regularization <= 50. else regularization/2
    results = [solve_ilt(dataset.b, y, replace(grid, regularization=lam), include=mask, sigma=sigma)
               for lam in (regularization, probe)]
    fitted, check = results
    weights = fitted.weights[:, None]
    peaks, intervals = identify_components(weights)
    check_peaks, _ = identify_components(check.weights[:, None])
    if (len(peaks) < 2 or len(peaks) != len(check_peaks) or min(peaks) == 0 or max(peaks) == len(weights)-1
            or min(check_peaks) == 0 or max(check_peaks) == len(weights)-1
            or np.max(np.abs(np.log(fitted.diffusion[list(peaks)]/check.diffusion[list(check_peaks)]))) > .15
            or fitted.rms >= .8*row['fit']['RMSE']):
        return None
    return SimpleNamespace(diffusion=fitted.diffusion, weights=weights, scales=[fitted.scale],
                           rms=[fitted.rms], intervals=intervals, components=peaks, converged=True)


def adaptive_local_ilt(dataset, rows, tolerance, balance, regularization, excluded=(), cancelled=None):
    if not (1 <= tolerance <= 30 and 1 <= balance <= 1000 and 1e-12 <= regularization <= 100):
        raise ValueError('Nieprawidłowe parametry Adaptive Local ILT.')
    if dataset.b is None:
        return LocalResult((), {}, tuple(row['id'] for row in rows), {})
    reference = dataset.options['reference_layer']
    by_id = {row['id']: row for row in rows}
    quality = {sid: info for sid, row in by_id.items() if sid not in excluded
               and (info := reference_quality(row, reference)) is not None
               and len(np.unique(dataset.b[row['fit']['fit_mask']])) >= MIN_DISTINCT_B}
    profiles, groups, attempted, withheld = {}, [], set(), {}

    def check_cancelled():
        if cancelled is not None and cancelled():
            raise SolveCancelled('Global ILT anulowano.')

    def project(ids, solution):
        return ilt_profiles(dataset.b, dataset.ppm, dataset.spectrum.real, reference, rows, ids,
                            solution, conserve_reference=True, withheld=withheld)

    for group in compatible_groups(dataset.b, [by_id[sid] for sid in quality], tolerance):
        check_cancelled()
        ids = group.ids
        attempted.update(ids)
        members = [by_id[sid] for sid in ids]
        amplitude, snr = np.array([quality[sid] for sid in ids]).T
        d = [row['fit']['D_m2_s'] for row in members]
        y, mask, sigma = fit_inputs(members)
        try:
            solution = solve_global_ilt(dataset.b, y, GlobalParameters(grid_for(members, regularization)),
                                        include=mask, sigma=sigma, cancelled=cancelled,
                                        roi_weights=amplitude_weights(amplitude, snr, balance))
            recovered, rejected = project(ids, solution)
        except SolveCancelled:
            raise
        except (ValueError, RuntimeError, np.linalg.LinAlgError):
            solution, recovered, rejected = None, {}, ids
        dominant = (solution.diffusion[np.argmax(solution.aggregate)] if solution is not None and solution.components else None)
        rms = (float(np.sqrt(np.mean((solution.rms/solution.scales)**2))) if solution is not None else None)
        groups.append(LocalGroup(ids, (min(d), max(d)), dominant, float(amplitude.max()/amplitude.min()),
                                 rms, solution is not None and solution.converged, rejected))
        profiles.update(recovered)

    # An ungrouped ROI may have useful independent evidence for a mixture.
    # Each ROI is reconstructed at most once, with no duplication across local solves.
    for sid in sorted(set(quality)-attempted):
        check_cancelled()
        row = by_id[sid]
        if row['fit']['RMSE'] < .003*quality[sid][0]:
            continue
        try:
            solution = separate_mixed_roi(dataset, row, regularization)
            if solution is None:
                continue
            recovered, _ = project((sid,), solution)
        except (ValueError, RuntimeError, np.linalg.LinAlgError):
            continue
        if sid in recovered and len(recovered[sid]) >= 2:
            d = row['fit']['D_m2_s']
            groups.append(LocalGroup((sid,), (d,d), float(solution.diffusion[np.argmax(solution.weights)]),
                                     1., solution.rms[0]/solution.scales[0], True, (), independent=True))
            profiles.update(recovered)
    membership = {}
    accepted = [g for g in groups if not g.independent and g.converged and g.dominant is not None
                and any(sid in profiles for sid in g.ids)]
    for group in groups:
        if not group.independent:
            continue
        sid = group.ids[0]
        for d, _ in profiles[sid]:
            matches = [g for g in accepted if 200*abs(d-g.dominant)/(d+g.dominant) <= tolerance]
            if len(matches) == 1:
                membership.setdefault(matches[0].ids, []).append((sid, d))
    return LocalResult(tuple(groups), profiles, tuple(sorted(set(by_id)-set(profiles)-set(excluded))),
                       {key: tuple(value) for key, value in membership.items()}, withheld)
