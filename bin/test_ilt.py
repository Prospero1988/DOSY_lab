from dataclasses import replace
import unittest
from unittest.mock import patch

import numpy as np

from ilt import D_UNIT, ILTParameters, solve_ilt


class InverseLaplaceTests(unittest.TestCase):
    b = np.linspace(0, 15e9, 48)
    parameters = ILTParameters(.4 * D_UNIT, 40 * D_UNIT)

    def decay(self, diffusion, amplitudes):
        return np.exp(-np.outer(self.b, np.array(diffusion)*D_UNIT)) @ amplitudes

    def test_mono_peak_and_grid_density(self):
        y = self.decay([4.], [2.5])
        for spacing, n, order in [('log', 50, 0), ('log', 200, 2), ('linear', 400, 2)]:
            with self.subTest(spacing=spacing, points=n, order=order):
                r = solve_ilt(self.b, y, replace(self.parameters, spacing=spacing, points=n, order=order))
                self.assertAlmostEqual(np.mean(r.density), 1.)
                self.assertEqual(len(r.diffusion), n)
                # A uniform linear grid spreads the same normalized curvature scale over
                # a wider physical interval near this D. Allow for regularization bias.
                np.testing.assert_allclose(r.diffusion[r.peaks[0]]/D_UNIT, 4., rtol=.12 if spacing == 'linear' else .10)
                self.assertLess(np.sqrt(np.mean((y-r.predict(self.b))**2))/y.max(), .015)

    def test_well_separated_biexponential_with_nonzero_minimum_b(self):
        b = self.b + 2e7
        truth = np.array([2., 17.])
        y = np.exp(-np.outer(b, truth*D_UNIT)) @ [.6, .4]
        r = solve_ilt(b, y, self.parameters)
        # Favorable high SNR, D ratio 8.5, and well-sampled initial and late decay.
        self.assertEqual(len(r.peaks), 2)
        np.testing.assert_allclose(np.sort(r.diffusion[list(r.peaks)])/D_UNIT, truth, rtol=.15)
        self.assertLess(r.rms, .004)

    def test_noise_and_regularization_tradeoff(self):
        y = self.decay([2., 17.], [.6, .4]) + np.random.default_rng(25).normal(0, .008, len(self.b))
        low = solve_ilt(self.b, y, replace(self.parameters, regularization=1e-8))
        high = solve_ilt(self.b, y, replace(self.parameters, regularization=.02))
        self.assertLess(low.rms, high.rms)
        self.assertGreater(low.penalty_norm, high.penalty_norm)
        self.assertGreater(np.linalg.norm(low.density-high.density)/np.sqrt(len(low.density)), .1)

    def test_close_components_below_noise_need_not_resolve(self):
        clean = self.decay([4., 4.5], [.5, .5])
        noise = .004
        self.assertLess(np.sqrt(np.mean((clean-self.decay([4.25], [1.]))**2)), noise)
        y = clean + np.random.default_rng(13).normal(0, noise, len(self.b))
        r = solve_ilt(self.b, y, self.parameters)
        self.assertEqual(len(r.peaks), 1)
        np.testing.assert_allclose(r.diffusion[r.peaks[0]]/D_UNIT, 4.25, rtol=.12)

    def test_small_signed_signal_with_mask_and_weights(self):
        y = self.decay([5.], [1e-120])
        y[-2] = -1e-122
        y[-1] = 1e-111  # an excluded point is visible but cannot set the scale
        original = y.copy()
        mask = np.arange(len(y)) < len(y)-1
        sigma = np.full(len(y), 1e-122)
        sigma[-1] = np.nan
        r = solve_ilt(self.b, y, self.parameters, include=mask, sigma=sigma)
        np.testing.assert_array_equal(y, original)
        self.assertEqual(r.scale, 1e-120)
        np.testing.assert_allclose(r.rms, np.sqrt(np.mean(((y-r.predict(self.b))[mask]/r.scale)**2))*r.scale, atol=0.)
        self.assertTrue(any('niedodatnie' in warning for warning in r.warnings))

    def test_validation_and_solver_failure(self):
        y = self.decay([4.], [1.])
        bad = y.copy(); bad[0] = np.nan
        for b, signal, kwargs in [(self.b, bad, {}), (self.b, -y, {}), (self.b, y, {'sigma': y*0}),
                                  (self.b, y, {'include': np.arange(len(y)) < 5}), (self.b*0, y, {}),
                                  (-self.b, y, {}), (self.b, y, {'include': [True]})]:
            with self.assertRaises(ValueError):
                solve_ilt(b, signal, self.parameters, **kwargs)
        with self.assertRaises(ValueError):
            solve_ilt(self.b, y, replace(self.parameters, d_min=self.parameters.d_max))
        with patch('ilt.nnls', side_effect=RuntimeError('iteration limit')):
            with self.assertRaisesRegex(RuntimeError, 'Solver ILT'):
                solve_ilt(self.b, y, self.parameters)
