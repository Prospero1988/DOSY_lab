"""Noise-aware, read-only checks of complex profiles, phase and baseline."""
import numpy as np

MIN_PROFILE_SNR = 5.
MIN_DISTINCT_B = 8
FAMILY_ALPHA = .05
CI_MULTIPLIER = 1.96


def estimate_spectral_noise(phased):
    """Estimate per-layer per-quadrature noise in quiet, locally detrended blocks."""
    if phased.shape[1] < 128:
        return None
    count = min(32, phased.shape[1]//32)
    width = phased.shape[1]//count
    block = phased[:, :count*width].reshape(len(phased), count, width)
    t = np.arange(width)/(width-1)
    quarter = width//4
    estimates = []
    for quadrature in (block.real, block.imag):
        left = np.median(quadrature[:, :, :quarter], axis=2)
        right = np.median(quadrature[:, :, -quarter:], axis=2)
        detrended = quadrature-(left[:, :, None]+(right-left)[:, :, None]*(t-.125)/.75)
        mad = np.median(abs(detrended-np.median(detrended,axis=2,keepdims=True)),axis=2)
        estimates.append(1.4826*mad)
    scales = np.maximum(*estimates)
    median = np.median(scales,axis=1)
    rank = np.median(scales/np.maximum(median[:,None],np.finfo(float).tiny),axis=0)
    quiet = np.argsort(rank,kind='stable')[:count//2]
    return np.median(scales[:,quiet],axis=1)


def layer_contrast(residual, predictor, signs):
    """Wild-sign diagnostic at the layer level, retaining spectral correlations."""
    contributions = predictor[:, None]*residual
    effect = np.linalg.norm(contributions.sum(axis=0))
    if effect <= 128*np.finfo(float).eps:
        return 0., 1.
    gram = (contributions@contributions.conj().T).real
    simulated_squared = np.sum((signs@gram)*signs, axis=1)
    return float(effect), float((1+np.count_nonzero(simulated_squared >= effect**2))/(1+len(signs)))


MESSAGES = {
    'shape':'Zespolony profil ROI zmienia się z b. D opisuje całe ROI; sprawdź nakładanie sygnałów lub ogon rozpuszczalnika.',
    'phase':'Faza ROI zmienia się z b. Sprawdź część urojoną i stabilność fazy między warstwami.',
    'baseline':'Wpływ korekcji bazy zmienia się z b. Sprawdź, czy maska obejmuje wyłącznie tło.',
    'insufficient':'Za mało informacji do oceny zmian profilu między warstwami.',
}


def assess_roi_profile(dataset, bounds, mask, noise_sigma=None):
    lo, hi = sorted(bounds)
    columns = (dataset.ppm >= lo) & (dataset.ppm <= hi)
    layer_indices = np.flatnonzero(mask)
    spectra = dataset.spectrum[mask][:,columns]
    norm = np.linalg.norm(spectra)
    diagnostic = dict(status='insufficient', warnings=[], flagged=[], checks={},
                      min_snr=MIN_PROFILE_SNR, min_distinct_b=MIN_DISTINCT_B,
                      active_layers=layer_indices+1, noise_source='unavailable')
    if spectra.shape[1] < 2 or not np.isfinite(norm) or norm == 0:
        return diagnostic
    spectra = spectra/norm
    _, _, vh = np.linalg.svd(spectra, full_matrices=False)
    template = vh[0]
    amplitude = abs(spectra@template.conj())
    usable = amplitude > 128*np.finfo(float).eps
    sigma = None
    if noise_sigma is not None:
        # Conservative coherent-projection noise scale, allowing correlated ppm bins.
        sigma = np.asarray(noise_sigma)[mask]*np.sqrt(spectra.shape[1])/norm
        known = np.isfinite(sigma)&(sigma>=0)
        usable &= known & (amplitude >= MIN_PROFILE_SNR*sigma)
        snr = np.divide(amplitude,sigma,out=np.full_like(amplitude,np.nan),where=known&(sigma>0))
        diagnostic.update(snr_lower_estimate=snr, noise_source='quiet_spectral_blocks',
                          excluded_low_information_layers=layer_indices[~usable]+1)
    b = dataset.b[layer_indices[usable]]
    diagnostic['diagnostic_layers'] = layer_indices[usable]+1
    diagnostic['diagnostic_b_span'] = float(np.ptp(b)) if len(b) else 0.
    if np.unique(b).size < MIN_DISTINCT_B:
        return diagnostic

    spectra = spectra[usable]
    weights = np.ones(len(b))
    if sigma is not None and np.all(sigma[usable]>0):
        weights = np.min(sigma[usable])/sigma[usable]
    weighted = spectra*weights[:,None]
    _, singular, vh = np.linalg.svd(weighted,full_matrices=False)
    template = vh[0]
    coefficients = spectra@template.conj()
    amplitude = abs(coefficients)
    weighted_amplitude = amplitude*weights
    gradient_order = (b-b.min())/np.ptp(b)
    predictor = weighted_amplitude*(gradient_order-np.average(gradient_order,weights=weighted_amplitude**2))
    leverage = np.linalg.norm(predictor)
    if leverage <= np.sqrt(np.finfo(float).eps)*np.linalg.norm(weighted_amplitude):
        return diagnostic
    predictor /= leverage
    signs = np.random.default_rng(21).choice((-1.,1.),size=(511,len(b)))

    rotation = np.divide(coefficients.conj(),amplitude,out=np.ones_like(coefficients),where=amplitude>0)
    aligned = spectra*rotation[:,None]
    phase_reference = np.angle(np.vdot(weighted_amplitude,coefficients*weights))
    projected = coefficients*np.exp(-1j*phase_reference)
    effects = {
        'shape': layer_contrast((aligned-amplitude[:,None]*template)*weights[:,None],predictor,signs),
        'phase': layer_contrast((projected.imag*weights)[:,None],predictor,signs),
    }
    if dataset.options['baseline_enabled']:
        correction = dataset.baseline[mask][:,columns][usable]*weights[:,None]/norm
        correction -= weighted_amplitude[:,None]*(weighted_amplitude@correction)/np.dot(weighted_amplitude,weighted_amplitude)
        effects['baseline'] = layer_contrast(correction,predictor,signs)
    threshold = FAMILY_ALPHA/len(effects)
    checks = {key:dict(effect=effect,p_heuristic=p,p_threshold=threshold) for key,(effect,p) in effects.items()}
    flagged = [key for key,check in checks.items() if check['p_heuristic']<=threshold]
    diagnostic.update(status='flagged' if flagged else 'no_detected_contrast',checks=checks,flagged=flagged,
                      relative_phase_deg=np.angle(projected,deg=True),
                      explained_energy=float(singular[0]**2/np.linalg.norm(weighted)**2),
                      warnings=[MESSAGES[key] for key in flagged])
    if sigma is not None:
        diagnostic['phase_noise_deg'] = np.degrees(np.arctan2(sigma[usable],amplitude))
    if 'phase' in flagged:
        diagnostic['phase_amplitude'] = amplitude
    if 'shape' in flagged:
        power = abs(template)**2
        split = min(max(1,np.searchsorted(np.cumsum(power),.5)+1),len(power)-1)
        halves = (slice(None,split),slice(split,None))
        if all(power[h].sum()>np.finfo(float).eps for h in halves):
            diagnostic['shape_amplitudes'] = [np.real(aligned[:,h]@template[h].conj())/power[h].sum() for h in halves]
    return diagnostic


def compare_diffusion(first, second):
    """Conservative screening scale for a difference without assuming fit independence."""
    if not first.get('ok') or not second.get('ok'):
        return dict(status='unavailable')
    d1,d2 = first['D_m2_s'],second['D_m2_s']
    s1,s2 = first['D_standard_error_m2_s'],second['D_standard_error_m2_s']
    bound = CI_MULTIPLIER*(s1+s2) if s1 is not None and s2 is not None else None
    change = d2-d1
    return dict(status='compared', D_first=d1, D_second=d2, D_change=change,
                screening_bound=bound, exceeds_bound=bool(bound is not None and abs(change)>bound
                                                           and not np.isclose(d1,d2,atol=0,rtol=1e-7)))


def qualify_roi(dataset, bounds, fit, result, sigma, noise, fit_decay):
    mask = fit['fit_mask']
    b = dataset.b
    diagnostic = assess_roi_profile(dataset,bounds,mask,noise)
    findings = []
    if diagnostic['status'] == 'insufficient':
        findings.append(dict(kind='information', severity='info', message=MESSAGES['insufficient']))
    indices = np.asarray(diagnostic.get('diagnostic_layers',()),dtype=int)-1
    shape = diagnostic.pop('shape_amplitudes',None)
    phase = diagnostic.pop('phase_amplitude',None)
    comparisons = {}
    if shape is not None:
        comparisons['shape'] = compare_diffusion(*(fit_decay(b[indices],y) for y in shape))
    if phase is not None:
        comparisons['phase'] = compare_diffusion(fit_decay(b[indices],result['corrected_integral'][indices]),
                                                 fit_decay(b[indices],phase))
    if dataset.options['baseline_enabled']:
        before = result['raw_integral'][mask]
        scale = np.max(abs(before),initial=0.)
        if scale>0:
            delta = result['corrected_integral'][mask]-before
            diagnostic['baseline_relative_rms'] = float(np.linalg.norm(delta/scale)/np.linalg.norm(before/scale))
            comparisons['baseline'] = compare_diffusion(fit_decay(b,result['raw_integral'],mask,sigma),fit)
    diagnostic['D_comparisons'] = comparisons
    for key,message in MESSAGES.items():
        if key == 'insufficient':
            continue
        matters = comparisons.get(key,{}).get('exceeds_bound',False)
        if key in diagnostic['flagged'] or matters:
            findings.append(dict(kind=key,severity='warning' if matters else 'info',message=message))
    return diagnostic, findings
