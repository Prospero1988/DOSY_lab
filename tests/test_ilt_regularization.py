from dataclasses import replace
from types import SimpleNamespace
import unittest

import numpy as np

from global_ilt import GlobalParameters, SolveCancelled
from ilt import ILTParameters
from ilt_regularization import Candidate, SweepInput, SweepOptions, objective_data, peak_tracks, scan


class RegularizationTests(unittest.TestCase):
    def fixture(self):
        b = np.linspace(0.,1.6e10,32)
        d = np.array([2.7,13.])*1e-10
        sigma = np.linspace(.002,.009,len(b))
        y = np.exp(-np.outer(b,d))@[1.,.6] + np.random.default_rng(0).normal(0,sigma)
        grid = ILTParameters(.2e-10,50e-10,50,regularization=1e-5)
        return SweepInput(b,y,np.ones(len(b),bool),grid,sigma,sigma)

    def test_weighted_discrepancy_and_scale_invariance(self):
        source = self.fixture()
        saved = source.y.copy()
        settings = SweepOptions(1e-10,1.,21,trust_noise=True,check_grid=True)
        result = scan(source,settings)
        self.assertEqual(result.method,'Zasada rozbieżności')
        selected = next(c for c in result.candidates if c.regularization==result.recommendation)
        self.assertLess(result.noise_band[0],selected.residual)
        self.assertLess(selected.residual,result.noise_band[1])
        self.assertAlmostEqual(selected.residual,selected.solution.weighted_rms)
        np.testing.assert_allclose(selected.maxima,[2.7e-10,13e-10],rtol=.15,atol=0.)
        self.assertEqual(len(result.grid_checks),2)
        scaled = replace(source,y=1e5*source.y,weights_sigma=1e5*source.weights_sigma,
                         noise_sigma=1e5*source.noise_sigma)
        check = scan(scaled,replace(settings,check_grid=False))
        self.assertAlmostEqual(result.expected,check.expected)
        self.assertEqual(result.recommendation,check.recommendation)
        np.testing.assert_array_equal(source.y,saved)

    def test_explicit_noise_and_systematic_qc(self):
        source = self.fixture()
        settings = SweepOptions(1e-9,1.,13,method='noise',trust_noise=True)
        for data,options in ((source,replace(settings,trust_noise=False)),
                             (replace(source,noise_sigma=None),settings),
                             (replace(source,systematic=True),settings),
                             (replace(source,noise_sigma=np.r_[np.nan,source.noise_sigma[1:]]),settings)):
            with self.subTest(data=data.systematic,trust=options.trust_noise):
                r = scan(data,options)
                self.assertIsNone(r.expected)
                self.assertIsNone(r.recommendation)

    def test_conditional_joint_curve_accounts_for_fixed_mu_and_roi_weights(self):
        s = self.fixture()
        grid = replace(s.parameters,order=0)
        data = SweepInput(s.b,np.column_stack((s.y,10*s.y)),np.column_stack((s.mask,s.mask)),
                          GlobalParameters(grid,sparsity=3e-4),
                          np.column_stack((s.noise_sigma,10*s.noise_sigma)),
                          np.column_stack((s.noise_sigma,10*s.noise_sigma)),np.array([1.,.4]))
        _,mask,scale,weights = objective_data(data)
        expected = np.linalg.norm((weights*data.noise_sigma/scale)[mask])
        r = scan(data,SweepOptions(1e-8,1e-2,9,trust_noise=True))
        self.assertEqual(r.fixed_mu,3e-4)
        self.assertAlmostEqual(r.expected,expected)
        valid = [c for c in r.candidates if c.valid]
        self.assertGreaterEqual(len(valid),7)
        for c in valid:
            self.assertGreater(c.sparsity_cost,0)
            self.assertAlmostEqual(c.objective_residual**2,c.residual**2+c.sparsity_cost)
            self.assertAlmostEqual(c.solution.objective,c.objective_residual**2+c.regularization*c.penalty**2)
            np.testing.assert_allclose(np.sum(c.masses,axis=0),c.solution.predict([0.])[0])

    def test_flat_or_failed_scan_does_not_force_selection(self):
        source = self.fixture()
        flat = replace(source,y=np.ones_like(source.y),noise_sigma=None)
        r = scan(flat,SweepOptions(1e-10,1.,9))
        self.assertIsNone(r.recommendation)
        self.assertEqual(r.interval,())
        with self.assertRaises(SolveCancelled):
            scan(flat,SweepOptions(1e-10,1.,9),cancelled=lambda: True)
        with self.assertRaises(ValueError):
            scan(flat,SweepOptions(1e-10,1.,3))

    def test_peak_merge_marks_both_tracks_ambiguous(self):
        d = np.geomspace(1e-10,1e-9,100)
        sol = SimpleNamespace(diffusion=d)
        c = [Candidate(1e-7,sol,maxima=tuple(d[[20,24]]),fractions=(.5,.5)),
             Candidate(1e-6,sol,maxima=(d[22],),fractions=(1.,))]
        tracks = peak_tracks(c)
        self.assertEqual(len(tracks),3)
        self.assertTrue(all(status=='Nierozstrzygnięte' and any(p[3] for p in points)
                            for status,points in tracks))
