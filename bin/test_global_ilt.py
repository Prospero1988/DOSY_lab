from dataclasses import replace
import unittest

import numpy as np

from global_ilt import GlobalParameters, SolveCancelled, identify_components, solve_global_ilt
from ilt import D_UNIT, ILTParameters, solve_ilt


class GlobalInverseLaplaceTests(unittest.TestCase):
    # Dense early b coverage retains fast components, while late b constrains slow ones.
    b = np.r_[0., np.geomspace(5e6, 35e9, 64)]
    p = GlobalParameters(ILTParameters(.3*D_UNIT, 40*D_UNIT, regularization=1e-6, order=0), tolerance=1e-6)
    partial = np.array([[1., .6, 0., 0.], [0., .4, 1., .55], [0., 0., 0., .45]])

    def signal(self, diffusion, amplitudes):
        return np.exp(-np.outer(self.b, np.array(diffusion)*D_UNIT)) @ amplitudes

    def fit(self, y, parameters=None, **kwargs):
        return solve_global_ilt(self.b, y, parameters or self.p, **kwargs)

    @staticmethod
    def peaks(result):
        return result.diffusion[list(result.components)]/D_UNIT

    def test_shared_diffusion_with_different_amplitudes_and_noise(self):
        y = self.signal([4.], [[.4, 1., 2., 4.]])
        sigma = np.broadcast_to([.004, .012, .009, .015], y.shape)
        y += sigma*np.random.default_rng(28).normal(size=y.shape)
        result = self.fit(y, sigma=sigma)
        self.assertTrue(result.converged)
        np.testing.assert_allclose(self.peaks(result), [4.], rtol=.08)
        np.testing.assert_allclose(result.fractions, np.ones((1, 4)))
        np.testing.assert_allclose(result.contributions.sum(axis=0), [.4, 1., 2., 4.], rtol=.04)

    def test_two_shared_populations_with_different_fractions(self):
        fractions = np.array([[.8, .2, .5], [.2, .8, .5]])
        y = self.signal([2., 17.], fractions)
        y += .001*np.random.default_rng(22).normal(size=y.shape)
        result = self.fit(y)
        # The log cells are about 5% apart; allow two cell spacings for regularization and noise.
        np.testing.assert_allclose(self.peaks(result), [2., 17.], rtol=.11)
        np.testing.assert_allclose(result.fractions, fractions, atol=.035)
        self.assertLess(result.rms.max(), .003)

    def test_partially_shared_components_do_not_force_roi_membership(self):
        result = self.fit(self.signal([1.5, 7., 24.], self.partial))
        np.testing.assert_allclose(self.peaks(result), [1.5, 7., 24.], rtol=.12)
        np.testing.assert_allclose(result.fractions, self.partial, atol=.065)
        self.assertLess(result.fractions[self.partial == 0].max(), .025)
        np.testing.assert_allclose(result.contributions.sum(axis=0), result.weights.sum(axis=0)*result.scales)
        self.assertEqual(result.intervals[0][0], 0)
        self.assertEqual(result.intervals[-1][1], self.p.grid.points)

    def test_noise_comparison_uses_identical_data_and_reports_tradeoff(self):
        clean = self.signal([4.], [[.4, 1., 2.]])
        for seed in (1, 2, 3, 4):
            with self.subTest(seed=seed):
                y = clean + .01*np.random.default_rng(seed).normal(size=clean.shape)
                joint = self.fit(y)
                independent = [solve_ilt(self.b, col, self.p.grid) for col in y.T]
                d_ind = [r.diffusion[r.peaks[0]]/D_UNIT for r in independent]
                rms_ind = np.array([r.rms for r in independent])
                np.testing.assert_allclose(self.peaks(joint), [4.], rtol=.08)
                np.testing.assert_allclose(d_ind, [4.]*3, rtol=.08)
                self.assertLess(np.max(abs(joint.rms-rms_ind)), .0004)
                # Compare recovery against the clean simulated curves, independently of residuals.
                for prediction in (joint.predict(self.b), np.column_stack([r.predict(self.b) for r in independent])):
                    self.assertLess(np.sqrt(np.mean((prediction-clean)**2)), .006)

    def test_close_diffusion_remains_unresolved_at_tested_noise(self):
        y = self.signal([4., 4.5], [[.6,.2,.5], [.4,.8,.5]])
        mono = np.exp(-np.outer(self.b, np.array([4.2,4.4,4.25])*D_UNIT))
        self.assertLess(np.sqrt(np.mean((y-mono)**2)), .003)
        y += .003*np.random.default_rng(22).normal(size=y.shape)
        result = self.fit(y)
        np.testing.assert_allclose(self.peaks(result), [4.25], rtol=.10)

    def test_penalty_sweep_reports_unstable_fast_maximum(self):
        y = self.signal([1.5,7.,24.], self.partial)
        y += .004*np.random.default_rng(14).normal(size=y.shape)
        sweep = [(1e-8, 1e-6), (1e-6, 1e-4), (1e-4, .003), (1e-3, .02)]
        solutions = [self.fit(y, replace(self.p, grid=replace(self.p.grid, regularization=lam), sparsity=mu))
                     for lam,mu in sweep]
        self.assertEqual([len(r.components) for r in solutions], [2, 3, 2, 1])
        np.testing.assert_allclose(self.peaks(solutions[1]), [1.5,7.,24.], rtol=.12)
        self.assertGreater(solutions[-1].rms.mean(), 3*solutions[1].rms.mean())
        self.assertGreater(solutions[1].fractions[-1,-1], .25)
        self.assertEqual(solutions[-1].fractions[0,-1], 1.)

    def test_scale_invariance_with_weights_masks_and_signed_input(self):
        y = self.signal([1.5,7.,24.], self.partial)
        y[-2, 2] = -.01
        mask = np.ones(y.shape, dtype=bool)
        mask[4, 1] = False
        sigma = np.broadcast_to(np.linspace(.01,.02,len(y))[:,None], y.shape).copy()
        sigma[4, 1] = np.nan
        original = y.copy()
        result = self.fit(y, include=mask, sigma=sigma)
        factors = np.array([1.,137.,1.,1.])
        scaled = self.fit(y*factors, include=mask, sigma=sigma*factors)
        np.testing.assert_array_equal(original, y)
        self.assertTrue(any('niedodatnie' in warning for warning in result.warnings))
        self.assertEqual(result.components, scaled.components)
        np.testing.assert_allclose(result.weights, scaled.weights, atol=1e-8)
        np.testing.assert_allclose(result.predict(self.b), scaled.predict(self.b)/factors, atol=1e-8)

    def test_mu_zero_agrees_with_separable_nnls_and_limits_are_reported(self):
        y = self.signal([2.,17.], [[.8,.2], [.2,.8]])
        p = replace(self.p, sparsity=0., tolerance=1e-8, max_iterations=40000)
        joint = self.fit(y, p)
        separate = np.column_stack([solve_ilt(self.b, col, p.grid).predict(self.b) for col in y.T])
        np.testing.assert_allclose(joint.predict(self.b), separate, atol=.00015, rtol=0.)
        limited = self.fit(y, replace(p, max_iterations=1))
        self.assertFalse(limited.converged)
        self.assertEqual(limited.iterations, 1)
        self.assertGreater(limited.objective, joint.objective)
        with self.assertRaises(SolveCancelled):
            self.fit(y, cancelled=lambda: True)
        with self.assertRaises(ValueError):
            self.fit(y, include=np.arange(len(y))[:,None] < 5)
        with self.assertRaises(ValueError):
            self.fit(-y)
        self.assertEqual(len(identify_components(np.ones((5, 2)))[0]), 1)
