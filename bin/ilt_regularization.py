"""Regularization scans on the original ILT objectives, separate from active fits."""
from dataclasses import dataclass, replace

import numpy as np
from scipy.signal import find_peaks
from scipy.stats import chi2

from global_ilt import GlobalParameters, SolveCancelled, component_intervals, solve_global_ilt
from ilt import GRID_SIZES, solve_ilt


@dataclass(frozen=True)
class SweepInput:
    b: np.ndarray
    y: np.ndarray
    mask: np.ndarray
    parameters: object
    weights_sigma: np.ndarray | None = None
    noise_sigma: np.ndarray | None = None
    roi_weights: np.ndarray | None = None
    systematic: bool = False

    @property
    def current(self):
        p = self.parameters
        return (p.grid if isinstance(p,GlobalParameters) else p).regularization


@dataclass(frozen=True)
class SweepOptions:
    minimum: float
    maximum: float
    count: int = 13
    method: str = 'auto'
    trust_noise: bool = False
    check_grid: bool = False


@dataclass
class Candidate:
    regularization: float
    solution: object = None
    residual: float = np.nan
    rms: float = np.nan
    penalty: float = np.nan
    sparsity_cost: float = 0.
    objective_residual: float = np.nan
    maxima: tuple = ()
    fractions: tuple = ()
    masses: tuple = ()
    boundary: bool = False
    status: str = ''

    @property
    def valid(self):
        return self.solution is not None and np.isfinite(self.residual)

    @property
    def eligible(self):
        return (self.valid and bool(self.maxima) and not self.boundary
                and getattr(self.solution,'effective_rank',3)>=3)


@dataclass
class SweepResult:
    candidates: tuple
    tracks: tuple
    recommendation: float | None
    interval: tuple
    method: str
    warnings: tuple
    expected: float | None
    noise_band: tuple
    curvature: np.ndarray
    grid_checks: tuple
    fixed_mu: float | None
    corner: float | None


def supplied_noise(entries):
    """Keep missing point errors as NaN, never substitute unit weights for noise."""
    columns = []
    for entry in entries:
        mask = np.asarray(entry['fit']['fit_mask'])
        sigma = np.asarray(entry.get('sigma'),dtype=float)
        columns.append(sigma if sigma.shape == mask.shape else np.full(mask.shape,np.nan))
    return np.column_stack(columns)


def has_systematic_warning(entries):
    return any(f['severity']=='warning' and f['kind'] in ('shape','phase','baseline')
               for row in entries for f in row['fit'].get('measurement_findings',()))


def objective_data(source):
    y = np.asarray(source.y,dtype=float)
    mask = np.asarray(source.mask,dtype=bool).reshape(y.shape)
    y, mask = y.reshape(len(y),-1), mask.reshape(len(y),-1)
    if not mask.any(axis=0).all() or not np.isfinite(y).all():
        raise ValueError('Brak skończonych całek lub aktywnych punktów do skanu λ.')
    scale = np.max(abs(y),axis=0,where=mask,initial=0.)
    if np.any(scale<=0):
        raise ValueError('Brak dodatniej skali całek do skanu λ.')
    weights = mask.astype(float)
    if source.weights_sigma is not None:
        sigma = np.asarray(source.weights_sigma,dtype=float).reshape(y.shape)
        if not (np.isfinite(sigma[mask]).all() and np.all(sigma[mask]>0)):
            raise ValueError('Wagi wymagają dodatnich, skończonych niepewności aktywnych całek.')
        minimum = np.min(sigma,axis=0,where=mask,initial=np.inf)
        np.divide(minimum,sigma,out=weights,where=mask)
    weights /= np.linalg.norm(weights,axis=0)
    if source.roi_weights is not None:
        weights *= np.sqrt(source.roi_weights)
    return y, mask, scale, weights


def expected_noise(source, mask, scale, weights, trust):
    if not trust or source.noise_sigma is None or source.systematic:
        return None, ()
    sigma = np.asarray(source.noise_sigma,dtype=float).reshape(mask.shape)
    if not (np.isfinite(sigma[mask]).all() and np.all(sigma[mask]>0)):
        return None, ()
    terms = ((weights*np.where(mask,sigma,0.))/scale)[mask]**2
    mean, squares = terms.sum(), np.dot(terms,terms)
    if not (mean>0 and squares>0 and np.isfinite([mean,squares]).all()):
        return None, ()
    factor, dof = squares/mean, mean**2/squares
    return float(np.sqrt(mean)), tuple(np.sqrt(factor*chi2.ppf([.025,.975],dof)))


def measure(source, grid, scale, weights, mask, y, cancelled):
    lam = grid.regularization
    record = Candidate(lam)
    joint = isinstance(source.parameters,GlobalParameters)
    try:
        if joint:
            p = replace(source.parameters,grid=grid)
            r = solve_global_ilt(source.b,source.y,p,include=source.mask,sigma=source.weights_sigma,
                                 roi_weights=source.roi_weights,cancelled=cancelled)
            record.solution = r
            if not r.converged:
                record.status = 'Brak zbieżności'
                return record
            indices, intervals = r.components, r.intervals
            penalty = np.sqrt(grid.points)*np.linalg.norm(r.weights)
            cost = p.sparsity*np.sqrt(r.weights.shape[1])*r.aggregate.sum()
            mass = r.contributions
            normalized_mass = np.array([r.weights[a:z].sum() for a,z in intervals])
        else:
            r = solve_ilt(source.b,source.y,grid,include=source.mask,sigma=source.weights_sigma)
            record.solution = r
            indices = tuple(sorted(r.peaks))
            intervals = component_intervals(r.weights,indices)
            penalty = r.penalty_norm
            cost = 0.
            normalized_mass = np.array([r.weights[a:z].sum() for a,z in intervals])
            mass = (normalized_mass*r.scale)[:,None]
    except SolveCancelled:
        raise
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        record.status = str(exc)
        return record
    prediction = r.predict(source.b).reshape(y.shape)
    normalized = np.where(mask,prediction/scale-y/scale,0.)
    record.residual = float(np.linalg.norm(weights*normalized))
    record.rms = float(np.linalg.norm(normalized)/np.sqrt(mask.sum()))
    record.penalty = float(penalty)
    record.sparsity_cost = float(cost)
    record.objective_residual = float(np.hypot(record.residual,np.sqrt(cost)))
    record.maxima = tuple(r.diffusion[list(indices)])
    if len(mass):
        # The joint aggregate matches the normalized columns used by the solver.
        record.fractions = tuple(normalized_mass/normalized_mass.sum())
        record.masses = tuple(tuple(row) for row in mass)
    record.boundary = any(i in (0,len(r.diffusion)-1) for i in indices)
    record.status = 'Zbieżny' if indices else 'Brak dodatnich udziałów'
    return record


def peak_tracks(candidates):
    """Match only mutually unique neighbors within 12 log-percent or 1.5 cells."""
    tracks, previous = [], []
    for index,c in enumerate(candidates):
        rates = np.asarray(c.maxima)
        current = []
        if len(rates):
            step = np.max(np.diff(np.log(c.solution.diffusion)))
            gate = max(.12,1.5*step)
            old = np.array([tracks[i][-1][1] for i in previous])
            close = abs(np.log(old[:,None]/rates)) <= gate
            conflicts = (close.sum(axis=1)>1) | np.any(close & (close.sum(axis=0)>1),axis=1)
            for i in np.flatnonzero(conflicts):
                point = tracks[previous[i]][-1]
                tracks[previous[i]][-1] = (*point[:3],True)
            for j,d in enumerate(rates):
                possible = np.flatnonzero(close[:,j])
                if len(possible)==1 and close[possible[0]].sum()==1:
                    track = previous[possible[0]]
                else:
                    track = len(tracks)
                    tracks.append([])
                current.append(track)
                ambiguous = bool(len(possible) and (len(possible)>1 or close[possible[0]].sum()>1))
                boundary = d in (c.solution.diffusion[0],c.solution.diffusion[-1])
                tracks[track].append((index,d,c.fractions[j],ambiguous or boundary))
        previous = current
    results = []
    for track in tracks:
        count = len(track)/len(candidates)
        span = np.ptp(np.log([p[1] for p in track]))
        if any(p[3] for p in track) or len(track)==1:
            status = 'Nierozstrzygnięte'
        elif count>=.8 and span<=.12:
            status = 'Stabilne w badanym zakresie'
        elif count>=.5 and span<=.30:
            status = 'Umiarkowanie wrażliwe'
        else:
            status = 'Silnie zależne od λ'
        results.append((status,tuple(track)))
    return tuple(results)


def find_corner(candidates):
    """Curvature of the log-log curve, parameterized by log λ."""
    curvature = np.full(len(candidates),np.nan)
    if not all(c.valid and c.penalty>0 and c.objective_residual>0 for c in candidates):
        return None, (), curvature
    t = np.log([c.regularization for c in candidates])
    x = np.log([c.objective_residual for c in candidates])
    y = np.log([c.penalty for c in candidates])
    dx, dy = np.gradient(x,t), np.gradient(y,t)
    speed = np.hypot(dx,dy)
    cross = dx*np.gradient(dy,t)-dy*np.gradient(dx,t)
    # Regularization must actually trade loss for penalty on both sides.
    good = (dx>0)&(dy<0)&(speed>.02)
    curvature[good] = cross[good]/speed[good]**3
    curvature[:2] = curvature[-2:] = np.nan
    maxima = find_peaks(np.r_[-np.inf,np.nan_to_num(curvature,nan=-np.inf),-np.inf],height=.1)[0]-1
    maxima = [i for i in maxima if candidates[i].eligible
              and x[-1]-x[i]>=np.log(1.15) and y[0]-y[i]>=np.log(1.15)
              and x[i]-x[0]>=np.log(1.05) and y[i]-y[-1]>=np.log(1.05)]
    if not maxima:
        return None, (), curvature
    best = max(maxima,key=lambda i:curvature[i])
    rival = [i for i in maxima if curvature[i]>=.75*curvature[best]]
    if len(rival)>1:
        return None, (min(rival),max(rival)), curvature
    return best, (best-1,best+1), curvature


def scan(source, options, cancelled=None):
    if not (1e-12 <= options.minimum < options.maximum <= 100 and 7 <= options.count <= 41
            and options.method in ('auto','corner','noise')):
        raise ValueError('Zakres skanu λ musi mieścić się w 10⁻¹²–100 i zawierać 7–41 prób.')
    y,mask,scale,weights = objective_data(source)
    params = source.parameters
    grid = params.grid if isinstance(params,GlobalParameters) else params
    mu = params.sparsity if isinstance(params,GlobalParameters) else None
    lambdas = np.geomspace(options.minimum,options.maximum,options.count)
    candidates = []
    for lam in lambdas:
        if cancelled is not None and cancelled():
            raise SolveCancelled('Anulowano skan λ.')
        candidates.append(measure(source,replace(grid,regularization=lam),scale,weights,mask,y,cancelled))
    tracks = peak_tracks(candidates)
    corner, region, curvature = find_corner(candidates)
    expected, band = expected_noise(source,mask,scale,weights,options.trust_noise)
    warnings = []
    if mu is not None:
        warnings.append('Dobór warunkowy przy stałym μ. Zmiana μ wymaga nowego skanu.')
    if source.systematic:
        warnings.append('QC wskazuje możliwy efekt systematyczny. Rozbieżności nie przypisano do szumu.')
    if expected is None and options.method != 'corner':
        warnings.append('Brak potwierdzonych, kompletnych σ(I) do zasady rozbieżności.')
    if any(not c.valid for c in candidates):
        warnings.append('Część prób nie osiągnęła zbieżności lub nie dała poprawnego rozwiązania.')
    if any(c.boundary for c in candidates):
        warnings.append('Część maksimów wypada przy granicy siatki D.')
    # Reject a nearly constant prediction relative to the observed attenuation.
    square = weights**2
    normalized = np.where(mask,y/scale,0.)
    constant = np.sum(square*normalized,axis=0)/square.sum(axis=0)
    reference = np.linalg.norm(weights*(normalized-constant))
    usable = [i for i,c in enumerate(candidates) if c.eligible and c.residual < .75*reference]
    chosen = None
    interval = ()
    method = 'Brak jednoznacznej rekomendacji'
    if expected is not None and options.method != 'corner':
        compatible = [i for i in usable if band[0]<=candidates[i].residual<=band[1]]
        crossings = [i for i,(a,b) in enumerate(zip(candidates,candidates[1:]))
                     if a.valid and b.valid and (a.residual<=expected) != (b.residual<=expected)]
        if compatible:
            interval = tuple(float(lambdas[i]) for i in (compatible[0],compatible[-1]))
        if len(crossings)==1 and candidates[crossings[0]].residual<=expected:
            i = crossings[0]
            near = [j for j in (i,i+1) if j in compatible]
            if near:
                def score(j):
                    unstable = sum(point[2] for status,track in tracks
                                   if status in ('Nierozstrzygnięte','Silnie zależne od λ')
                                   for point in track if point[0]==j)
                    return unstable,abs(np.log(candidates[j].residual/expected))

                chosen = min(near,key=score)
                method = 'Zasada rozbieżności'
                if corner is not None and corner not in compatible:
                    warnings.append('Załamanie krzywej wypada poza zakresem zgodnym z szumem.')
        if chosen is None:
            warnings.append('Skan nie wyznaczył jednoznacznego przejścia przez poziom szumu.')
            if corner is not None:
                warnings.append('Załamanie krzywej jest dostępne jako osobna wskazówka.')
    elif options.method != 'noise':
        if region:
            interval = tuple(float(lambdas[i]) for i in region)
        if corner in usable:
            chosen = corner
            method = 'Heurystyka krzywej warunkowej' if mu else 'Heurystyka L-curve'
        else:
            warnings.append('Brak wyraźnego wewnętrznego załamania krzywej lub za mały kontrast zaniku.')
    if chosen is not None:
        # A recommendation is still local to the grid. Flag short or ambiguous tracks.
        if any(status not in ('Stabilne w badanym zakresie','Umiarkowanie wrażliwe')
               and any(p[0]==chosen for p in track) for status,track in tracks):
            warnings.append('Część maksimów przy sugerowanym λ jest niestabilna lub niejednoznaczna.')
    checks = []
    if options.check_grid:
        lam = lambdas[chosen] if chosen is not None else source.current
        alternative = next((size for size in GRID_SIZES if size>grid.points),GRID_SIZES[-2])
        for label,test in (
                ('Szerszy zakres D',replace(grid,d_min=grid.d_min/2,d_max=grid.d_max*2,regularization=lam)),
                ('Inna rozdzielczość siatki D',replace(grid,points=alternative,regularization=lam))):
            if cancelled is not None and cancelled():
                raise SolveCancelled('Anulowano skan λ.')
            checks.append((label,measure(source,test,scale,weights,mask,y,cancelled)))
    return SweepResult(tuple(candidates),tracks,None if chosen is None else float(lambdas[chosen]),interval,
                       method,tuple(warnings),expected,band,curvature,tuple(checks),mu,
                       None if corner is None else float(lambdas[corner]))
