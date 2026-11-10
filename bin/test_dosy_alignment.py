import unittest

import numpy as np
from scipy.integrate import trapezoid

from core import fit_decay
from dosy_alignment import alignment_plan, cluster_diffusion, ilt_profiles
from global_ilt import GlobalParameters, solve_global_ilt
from ilt import D_UNIT, ILTParameters


class DiffusionAlignmentTests(unittest.TestCase):
    b = np.r_[0., np.geomspace(5e6, 30e9, 64)]

    def row(self, sid, rate, amplitude=1., noise=0., seed=0):
        y = amplitude*(np.exp(-self.b*rate*D_UNIT) + noise*np.random.default_rng(seed).normal(size=len(self.b)))
        return dict(id=sid, fit=fit_decay(self.b, y), result={'corrected_integral': y})

    def test_tolerance_is_deterministic_and_does_not_chain(self):
        rows = [self.row('A',1.), self.row('B',1.04), self.row('C',1.08), self.row('D',3.), self.row('E',3.03)]
        self.assertEqual([g.ids for g in cluster_diffusion(rows, 5.)], [('A','B'), ('D','E')])
        self.assertEqual(cluster_diffusion(rows, 5.), cluster_diffusion(rows[::-1], 5.))
        self.assertEqual([g.ids for g in cluster_diffusion(rows, 10.)], [('A','B','C'), ('D','E')])
        self.assertEqual([g.ids for g in cluster_diffusion(rows, 5., {'B'})], [('D','E')])
        rows[4]['fit']['R2'] = .9
        self.assertEqual([g.ids for g in cluster_diffusion(rows, 5.)], [('A','B')])
        rows[1]['fit']['D_standard_error_m2_s'] = None
        group, = cluster_diffusion(rows, 5.)
        self.assertAlmostEqual(group.representative/D_UNIT, 1.02, places=5)

    def test_shared_fit_improves_average_recovery_with_unequal_amplitudes(self):
        error_independent, error_joint = [], []
        for seed in range(8):
            rows = [self.row(str(i), 3., 10**i, noise, seed*13+i) for i,noise in enumerate((.009,.012,.01))]
            group, = alignment_plan(self.b, rows, 15., set(), shared=True)
            self.assertFalse(group.incompatible)
            self.assertEqual(len(group.diagnostics), 3)
            error_independent.extend(abs(row['fit']['D_m2_s']/D_UNIT-3.) for row in rows)
            error_joint.append(abs(group.shared/D_UNIT-3.))
        self.assertLess(np.mean(error_joint), np.mean(error_independent))

    def test_close_distinct_decays_can_reject_shared_d(self):
        rows = [self.row('A',4., noise=.001, seed=2), self.row('B',4.03, amplitude=10, noise=.001, seed=3),
                self.row('C',4.55, amplitude=.2, noise=.004, seed=4)]
        before = [row['fit']['D_m2_s'] for row in rows]
        group, = alignment_plan(self.b, rows, 20., set(), shared=True)
        self.assertIn('C', group.incompatible)
        self.assertGreater(group.diagnostics[2][1], group.diagnostics[2][0])
        self.assertEqual([row['fit']['D_m2_s'] for row in rows], before)

    def test_ilt_reconstructs_distinct_profiles_inside_one_roi(self):
        ppm = np.linspace(2., 5., 501)
        rates = np.array([2., 16.])*D_UNIT
        slow = np.exp(-((ppm-2.75)/.07)**2) + .8*np.exp(-((ppm-4.25)/.06)**2)
        fast = .9*np.exp(-((ppm-2.83)/.075)**2)
        spectra = np.exp(-np.outer(self.b,rates)) @ [slow,fast]
        rows = []
        for sid,bounds in [('mixed',[2.4,3.2]),('slow',[4.,4.55])]:
            start, stop = np.searchsorted(ppm, bounds)
            y = trapezoid(spectra[:,start:stop], x=ppm[start:stop], axis=1)
            rows.append(dict(id=sid, bounds=bounds, fit=fit_decay(self.b, y), result={'corrected_integral': y}))
        y = np.column_stack([row['result']['corrected_integral'] for row in rows])
        p = GlobalParameters(ILTParameters(.3*D_UNIT, 40*D_UNIT, regularization=1e-6, order=0), tolerance=1e-6)
        result = solve_global_ilt(self.b, y, p)
        profiles, rejected = ilt_profiles(self.b, ppm, spectra, 0, rows, ('mixed','slow'), result)
        self.assertFalse(rejected)
        self.assertEqual(len(profiles['mixed']), 2)
        self.assertEqual(len(profiles['slow']), 1)
        start, stop = np.searchsorted(ppm, rows[0]['bounds'])
        grid = ppm[start:stop]
        for (d, trace), expected_d, center in zip(profiles['mixed'], rates, [2.75,2.83]):
            self.assertLess(abs(grid[trace.argmax()]-center), .015)
            self.assertLess(abs(d/expected_d-1), .15)
        reconstructed = sum(trace for _,trace in profiles['mixed'])
        self.assertLess(np.linalg.norm(reconstructed-spectra[0,start:stop])/np.linalg.norm(reconstructed), .04)
