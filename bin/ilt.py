"""Nonnegative Tikhonov inversion of a single measured DOSY attenuation."""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import svdvals
from scipy.optimize import nnls
from scipy.signal import find_peaks


MIN_DISTINCT_B = 6
GRID_SIZES = (50, 100, 200, 400)
D_UNIT = 1e-10


@dataclass(frozen=True)
class ILTParameters:
    d_min: float
    d_max: float
    points: int = 100
    spacing: str = 'log'
    regularization: float = 1e-5
    order: int = 2


@dataclass(frozen=True)
class ILTResult:
    diffusion: np.ndarray
    weights: np.ndarray  # per-bin contributions at b=0, in units of scale
    scale: float
    peaks: tuple[int, ...]
    rms: float  # in original integral units
    weighted_rms: float  # in normalized intensity units, relative weights
    penalty_norm: float
    effective_rank: int
    warnings: tuple[str, ...]

    @property
    def density(self):
        """Unit-area density in the selected, normalized grid coordinate u."""
        total = self.weights.sum()
        return (self.weights / total) * len(self.weights) if total else self.weights.copy()

    def predict(self, b):
        return (_kernel(b, self.diffusion) @ self.weights) * self.scale


def _kernel(b, diffusion):
    with np.errstate(over='ignore', under='ignore'):
        return np.exp(-np.multiply.outer(np.asarray(b) * D_UNIT, diffusion / D_UNIT))


def diffusion_grid(p):
    if not (np.isfinite([p.d_min, p.d_max]).all() and 0 < p.d_min < p.d_max
            and np.isfinite(p.d_max / D_UNIT)):
        raise ValueError('Zakres D musi spełniać 0 < D min < D max.')
    if p.points not in GRID_SIZES or p.spacing not in ('log', 'linear') or p.order not in (0, 2):
        raise ValueError('Nieprawidłowa siatka lub tryb regularyzacji ILT.')
    if not (np.isfinite(p.regularization) and 1e-12 <= p.regularization <= 100):
        raise ValueError('λ musi mieścić się w zakresie 10⁻¹²–100.')
    n = p.points
    u = (np.arange(n) + .5) / n
    if p.spacing == 'log':
        diffusion = np.exp(np.log(p.d_min) + u * (np.log(p.d_max) - np.log(p.d_min)))
    else:
        diffusion = p.d_min + u * (p.d_max - p.d_min)
    if np.any(np.diff(diffusion) <= 0):
        raise ValueError('Zakres D jest za wąski dla wybranej siatki.')
    return diffusion


def _grid_and_penalty(p):
    diffusion = diffusion_grid(p)
    n = p.points
    identity = np.eye(n)
    # Each cell has width 1/n in u. The amplitude density equals n times the bin mass.
    # The 0.05-u length scale keeps the curvature penalty usable on a similar λ scale.
    penalty = np.sqrt(n) * (identity if p.order == 0 else np.diff(identity, n=2, axis=0) * (.05*n)**2)
    return diffusion, penalty


def solve_ilt(b, intensity, parameters, *, include=None, sigma=None):
    diffusion, penalty = _grid_and_penalty(parameters)
    b = np.asarray(b, dtype=float)
    y = np.asarray(intensity, dtype=float)
    if b.ndim != 1 or b.shape != y.shape or not np.isfinite(b).all() or np.any(b < 0):
        raise ValueError('ILT wymaga zweryfikowanych, skończonych wartości b zgodnych z całkami.')
    if not np.isfinite(y).all():
        raise ValueError('Całki zawierają NaN/Inf; nie usuwamy ich automatycznie.')
    mask = np.ones(b.shape, dtype=bool) if include is None else np.asarray(include)
    if mask.shape != b.shape or not np.isin(mask, [True, False]).all():
        raise ValueError('Nieprawidłowa maska punktów.')
    mask = mask.astype(bool)
    active_b, active_y = b[mask], y[mask]
    if len(np.unique(active_b)) < MIN_DISTINCT_B:
        raise ValueError('ILT wymaga co najmniej 6 różnych aktywnych wartości b.')
    if active_y.max() <= 0:
        raise ValueError('Całki są niedodatnie; sprawdź fazowanie i bazę.')
    scale = float(np.max(np.abs(active_y)))
    row_weights = np.ones(active_b.size)
    warnings = []
    if sigma is not None:
        uncertainty = np.asarray(sigma, dtype=float)
        if (uncertainty.shape != b.shape or not np.isfinite(uncertainty[mask]).all()
                or np.any(uncertainty[mask] <= 0)):
            raise ValueError('Wagi wymagają dodatnich, skończonych niepewności aktywnych całek.')
        active_sigma = uncertainty[mask]
        row_weights = active_sigma.min() / active_sigma
        row_weights /= np.sqrt(np.mean(row_weights**2))
        if np.median(active_sigma) > .1 * scale:
            warnings.append('Duża niepewność całek względem amplitudy sygnału.')
    if np.any(active_y <= 0):
        warnings.append('Uwzględniono niedodatnie całki bez przycinania i odrzucania.')
    kernel = _kernel(active_b, diffusion)
    row_scale = row_weights / np.sqrt(active_b.size)
    design = kernel * row_scale[:, None]
    observed = active_y / scale
    augmented = np.vstack((design, np.sqrt(parameters.regularization) * penalty))
    target = np.r_[observed * row_scale, np.zeros(penalty.shape[0])]
    try:
        weights, _ = nnls(augmented, target, maxiter=8 * parameters.points)
    except (RuntimeError, np.linalg.LinAlgError, ValueError) as exc:
        raise RuntimeError('Solver ILT nie uzyskał rozwiązania.') from exc
    if not np.isfinite(weights).all():
        raise RuntimeError('Solver ILT nie uzyskał rozwiązania.')
    residual = observed - kernel @ weights
    singular = svdvals(design)
    rank = int(np.count_nonzero(singular > singular[0] * 1e-6))
    if rank < 3:
        warnings.append('Dane słabo rozróżniają współczynniki D w tym zakresie.')
    if _kernel([active_b.min()], diffusion)[0, -1] < .05:
        warnings.append('Duże D są silnie stłumione już przy najniższym aktywnym b.')
    peaks = ()
    if weights.max() > 0:
        dominant = int(np.argmax(weights))
        candidates = find_peaks(weights, prominence=.05 * weights[dominant])[0]
        peaks = tuple(sorted({dominant, *map(int, candidates)}, key=lambda i: weights[i], reverse=True))
        if dominant in (0, len(weights)-1):
            warnings.append('Dominujące maksimum leży przy granicy siatki D.')
    else:
        warnings.append('Rozwiązanie ILT nie zawiera dodatnich udziałów.')
    return ILTResult(diffusion, weights, scale, peaks, float(np.sqrt(np.mean(residual**2)) * scale),
                     float(np.linalg.norm(residual * row_scale)), float(np.linalg.norm(penalty @ weights)),
                     rank, tuple(warnings))
