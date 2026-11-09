"""Joint nonnegative DOSY inversion with diffusion-row sparsity."""
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.linalg import eigvalsh
from scipy.signal import find_peaks

from ilt import ILTParameters, MIN_DISTINCT_B, diffusion_grid, _kernel


@dataclass(frozen=True)
class GlobalParameters:
    grid: ILTParameters  # zeroth-order Tikhonov for this solver
    sparsity: float = 1e-4
    max_iterations: int = 10000
    tolerance: float = 2e-6


@dataclass(frozen=True)
class GlobalResult:
    diffusion: np.ndarray
    weights: np.ndarray  # normalized bin masses, rows D, columns ROI
    scales: np.ndarray
    components: tuple[int, ...]
    intervals: tuple[tuple[int, int], ...]
    contributions: np.ndarray  # original intensity units extrapolated to b=0
    rms: np.ndarray  # original intensity units, included points only
    objective: float
    stationarity: float
    iterations: int
    converged: bool
    warnings: tuple[str, ...]

    @property
    def aggregate(self):
        return np.linalg.norm(self.weights, axis=1)

    @property
    def fractions(self):
        totals = self.contributions.sum(axis=0)
        return np.divide(self.contributions, totals, out=np.zeros_like(self.contributions), where=totals > 0)

    def predict(self, b):
        return (_kernel(b, self.diffusion) @ self.weights) * self.scales


class SolveCancelled(RuntimeError):
    """An interactive parameter change superseded this solve."""


def prepare_attenuation(b, intensity, include=None, sigma=None):
    b, y = np.asarray(b, dtype=float), np.asarray(intensity, dtype=float)
    if (b.ndim != 1 or y.ndim != 2 or y.shape[0] != b.size or y.shape[1] < 2
            or not np.isfinite(b).all() or np.any(b < 0)):
        raise ValueError('Global ILT wymaga wspólnej listy b i co najmniej dwóch całek.')
    if not np.isfinite(y).all():
        raise ValueError('Całki zawierają NaN/Inf; nie usuwamy ich automatycznie.')
    mask = np.ones(y.shape, dtype=bool) if include is None else np.asarray(include)
    if mask.shape != y.shape or not np.isin(mask, [True, False]).all():
        raise ValueError('Nieprawidłowa maska punktów.')
    mask = mask.astype(bool)
    if any(np.unique(b[active]).size < MIN_DISTINCT_B for active in mask.T):
        raise ValueError('Każde ROI wymaga co najmniej 6 różnych aktywnych wartości b.')
    if np.any(np.max(y, axis=0, where=mask, initial=-np.inf) <= 0):
        raise ValueError('Całki są niedodatnie; sprawdź fazowanie i bazę.')
    scales = np.max(np.abs(y), axis=0, where=mask, initial=0.)
    relative = mask.astype(float)
    if sigma is not None:
        sigma = np.asarray(sigma, dtype=float)
        if sigma.shape != y.shape or not np.isfinite(sigma[mask]).all() or np.any(sigma[mask] <= 0):
            raise ValueError('Wagi wymagają dodatnich, skończonych niepewności aktywnych całek.')
        minima = np.min(sigma, axis=0, where=mask, initial=np.inf)
        np.divide(minima, sigma, out=relative, where=mask)
    relative /= np.linalg.norm(relative, axis=0)
    normalized = np.divide(y, scales, out=np.zeros_like(y), where=mask)
    return b, normalized, mask, scales, relative**2


def nonnegative_group_prox(values, threshold):
    """Project positive, then shrink each diffusion row by its L2 norm."""
    positive = np.maximum(values, 0.)
    norms = np.linalg.norm(positive, axis=1)
    shrink = np.divide(np.maximum(norms-threshold, 0.), norms, out=np.zeros_like(norms), where=norms > 0)
    return positive * shrink[:, None]


def identify_components(weights):
    aggregate = np.linalg.norm(weights, axis=1)
    maximum = aggregate.max()
    if maximum <= 0:
        return (), ()
    detected = find_peaks(np.r_[0., aggregate, 0.], prominence=.05*maximum)[0] - 1
    peaks = tuple(map(int, detected)) or (int(np.argmax(aggregate)),)
    cuts = [0]
    for left, right in zip(peaks, peaks[1:]):
        valley = aggregate[left+1:right]
        if valley.size:
            minima = np.flatnonzero(valley == valley.min())
            cuts.append(left + 2 + int(minima[len(minima)//2]))
        else:
            cuts.append(right)
    cuts.append(len(weights))
    return peaks, tuple(zip(cuts, cuts[1:]))


def solve_global_ilt(b, intensity, parameters, *, include=None, sigma=None,
                     cancelled: Callable[[], bool] | None = None):
    """Accelerated proximal gradient on one coupled convex objective."""
    p = parameters
    if p.grid.order != 0 or not (np.isfinite(p.sparsity) and 0 <= p.sparsity <= 100
            and isinstance(p.max_iterations, (int, np.integer)) and 1 <= p.max_iterations <= 100000
            and 1e-10 <= p.tolerance <= .01):
        raise ValueError('Nieprawidłowe parametry solvera Global ILT.')
    diffusion = diffusion_grid(p.grid)
    b, y, mask, scales, w2 = prepare_attenuation(b, intensity, include, sigma)
    kernel = _kernel(b, diffusion)
    ridge = p.grid.regularization * p.grid.points  # L = sqrt(P) I
    group = p.sparsity * np.sqrt(y.shape[1])
    rhs = kernel.T @ (w2*y)
    n = len(b)
    gram = kernel @ kernel.T
    lipschitz = 0.
    for column in w2.T:
        w = np.sqrt(column)
        largest = eigvalsh(gram * w[:, None] * w[None, :], subset_by_index=(n-1, n-1))[0]
        lipschitz = max(lipschitz, float(largest))
    step = 1 / (2 * (lipschitz + ridge))
    gradient_scale = max(np.linalg.norm(2*rhs), np.finfo(float).eps)

    def iterate(point):
        gradient = 2 * (kernel.T @ (w2 * (kernel @ point)) - rhs + ridge * point)
        return nonnegative_group_prox(point-step*gradient, step*group)

    weights = np.zeros_like(rhs)
    accelerated, momentum, stationarity = weights, 1., np.inf
    for iteration in range(1, p.max_iterations+1):
        if cancelled is not None and cancelled():
            raise SolveCancelled('Global ILT anulowano.')
        candidate = iterate(accelerated)
        change = candidate-weights
        if np.sum((accelerated-candidate)*change) > 0:
            momentum = 1.
            accelerated = candidate
        else:
            next_momentum = (1 + np.sqrt(1 + 4*momentum**2))/2
            accelerated = candidate + ((momentum-1)/next_momentum)*change
            momentum = next_momentum
        weights = candidate
        if iteration % 20 == 0 or iteration == p.max_iterations:
            stationarity = float(np.linalg.norm(weights-iterate(weights))/(step*gradient_scale))
            if stationarity <= p.tolerance:
                break

    residual = np.where(mask, kernel @ weights-y, 0.)
    rms = np.sqrt(np.sum(residual**2, axis=0)/mask.sum(axis=0))*scales
    aggregate = np.linalg.norm(weights, axis=1)
    objective = float(np.sum(w2*residual**2) + ridge*np.sum(weights**2) + group*aggregate.sum())
    peaks, intervals = identify_components(weights)
    contributions = np.array([weights[a:z].sum(axis=0)*scales for a,z in intervals]).reshape(-1, len(scales))
    warnings = []
    if np.any(np.asarray(intensity)[mask] <= 0):
        warnings.append('Uwzględniono niedodatnie całki bez przycinania i odrzucania.')
    if any(peak in (0, len(diffusion)-1) for peak in peaks):
        warnings.append('Maksimum przy granicy siatki D.')
    if not peaks:
        warnings.append('Brak dodatnich udziałów. Sprawdź zakres D, λ i μ.')
    elif np.any(weights.sum(axis=0) == 0):
        warnings.append('Co najmniej jeden ROI nie ma dodatnich udziałów.')
    if any(np.exp(-b[active].min()*diffusion[-1]) < .05 for active in mask.T):
        warnings.append('Duże D są silnie stłumione już przy najniższym aktywnym b.')
    return GlobalResult(diffusion, weights, scales, peaks, intervals, contributions, rms, objective,
                        stationarity, iteration, stationarity <= p.tolerance, tuple(warnings))
