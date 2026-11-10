"""Read-only checks of layer-dependent complex profiles, phase and baseline."""
import numpy as np


def layer_contrast(residual, predictor, signs):
    """Wild-sign diagnostic at the layer level, retaining spectral correlations."""
    contributions = predictor[:, None]*residual
    effect = np.linalg.norm(contributions.sum(axis=0))
    if effect <= 128*np.finfo(float).eps:
        return 0., 1.
    gram = (contributions@contributions.conj().T).real
    simulated_squared = np.sum((signs@gram)*signs, axis=1)
    return float(effect), float((1+np.count_nonzero(simulated_squared >= effect**2))/(1+len(signs)))


def assess_roi_profile(dataset, bounds, mask):
    lo, hi = sorted(bounds)
    columns = (dataset.ppm >= lo) & (dataset.ppm <= hi)
    b = dataset.b[mask]
    if np.unique(b).size < 8 or np.count_nonzero(columns) < 2:
        return dict(status='insufficient', warnings=[])
    spectra = dataset.spectrum[mask][:, columns]
    norm = np.linalg.norm(spectra)
    if not np.isfinite(norm) or norm == 0:
        return dict(status='insufficient', warnings=[])
    spectra = spectra/norm
    # Estimate the complex template from all active layers, independently of b.
    _, singular, vh = np.linalg.svd(spectra, full_matrices=False)
    template = vh[0]
    coefficients = spectra@template.conj()
    amplitude = abs(coefficients)
    if not np.any(amplitude):
        return dict(status='insufficient', warnings=[])
    gradient_order = (b-b.min())/np.ptp(b)
    predictor = amplitude*(gradient_order-np.average(gradient_order, weights=amplitude**2))
    leverage = np.linalg.norm(predictor)
    if leverage <= np.sqrt(np.finfo(float).eps)*np.linalg.norm(amplitude):
        return dict(status='insufficient', warnings=[])
    predictor /= leverage
    signs = np.random.default_rng(21).choice((-1.,1.), size=(511,len(b)))

    # Rotate only the diagnostic traces. The integral and processing are untouched.
    rotation = np.divide(coefficients.conj(), amplitude, out=np.ones_like(coefficients), where=amplitude>0)
    aligned = spectra*rotation[:,None]
    shape = aligned-amplitude[:,None]*template
    shape_effect, shape_p = layer_contrast(shape, predictor, signs)
    phase_reference = np.angle(np.vdot(amplitude, coefficients))
    projected = coefficients*np.exp(-1j*phase_reference)
    phase_effect, phase_p = layer_contrast(projected.imag[:,None], predictor, signs)
    relative_phase = np.angle(projected, deg=True)
    checks = dict(shape=dict(effect=shape_effect, p_heuristic=shape_p),
                  phase=dict(effect=phase_effect, p_heuristic=phase_p))
    if dataset.options['baseline_enabled']:
        correction = dataset.baseline[mask][:,columns]/norm
        # Remove the component proportional to the measured attenuation.
        correction = correction-amplitude[:,None]*(amplitude@correction)/np.dot(amplitude,amplitude)
        effect, probability = layer_contrast(correction,predictor,signs)
        checks['baseline'] = dict(effect=effect, p_heuristic=probability)

    messages = {
        'shape':'Zespolony profil ROI zmienia się z b. D opisuje całe ROI; sprawdź nakładanie sygnałów lub ogon rozpuszczalnika.',
        'phase':'Faza ROI zmienia się z b. Sprawdź część urojoną i stabilność fazy między warstwami.',
        'baseline':'Wpływ korekcji bazy zmienia się z b. Sprawdź, czy maska obejmuje wyłącznie tło.',
    }
    flagged = [key for key,check in checks.items() if check['p_heuristic'] <= .05/len(checks)]
    return dict(status='flagged' if flagged else 'no_detected_contrast', checks=checks, flagged=flagged,
                active_layers=np.flatnonzero(mask)+1, relative_phase_deg=relative_phase,
                explained_energy=float(singular[0]**2), warnings=[messages[k] for k in flagged])
