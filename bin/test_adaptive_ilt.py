from types import SimpleNamespace
import unittest

import numpy as np
from scipy.integrate import trapezoid

from adaptive_ilt import adaptive_local_ilt, amplitude_weights, fit_inputs, grid_for, local_sweep_inputs
from core import fit_decay
from global_ilt import GlobalParameters, solve_global_ilt
from ilt import D_UNIT


class AdaptiveILTTests(unittest.TestCase):
    def fixture(self, settings):
        b = np.r_[0., np.geomspace(5e6, 35e9, 64)]
        ppm = np.linspace(0., len(settings)+1, (len(settings)+1)*120+1)
        spectra = np.zeros((len(b), len(ppm)))
        rows = []
        for i, (rates, amps, noise) in enumerate(settings, 1):
            bounds = [i-.42, i+.42]
            start, stop = np.searchsorted(ppm, bounds)
            x = ppm[start:stop]
            centers = i + np.arange(len(rates))*.08
            profiles = np.exp(-((x[None,:]-centers[:,None])/.08)**2)*np.array(amps)[:,None]
            values = np.exp(-np.outer(b, np.array(rates)*D_UNIT)) @ profiles
            values += noise*max(amps)*np.random.default_rng(i).normal(size=values.shape)
            spectra[:,start:stop] = values
            intensity = trapezoid(values, x=x, axis=1)
            rows.append(dict(id=str(i), bounds=bounds, fit=fit_decay(b,intensity), result={'corrected_integral':intensity}))
        return SimpleNamespace(b=b, ppm=ppm, spectrum=spectra, options={'reference_layer':0}), rows

    def test_amplitude_balance_changes_weights_without_splitting_reliable_weak_roi(self):
        amplitude = np.array([100., 1., .3])
        snr = np.array([100., 120., 20.])
        np.testing.assert_array_equal(amplitude_weights(amplitude, snr, 10.), [1., 1., .25])
        np.testing.assert_array_equal(amplitude_weights(amplitude, snr, 1000.), [1., 1., 1.])

    def test_local_preserves_close_groups_that_whole_global_merges(self):
        dataset, rows = self.fixture([([3.],[1.],0.), ([3.],[.01],.01), ([3.25],[1.],0.),
                                      ([3.25],[.5],.003), ([9.],[2.],.001), ([9.],[1.],.001), ([9.],[1e-4],3.)])
        rows[0]['sigma'] = [None]  # an incomplete optional uncertainty vector falls back to RMS
        baseline = [row['fit']['D_m2_s'] for row in rows]
        local = adaptive_local_ilt(dataset, rows, 5., 10., 1e-6)
        self.assertEqual([g.ids for g in local.groups], [('1','2'), ('3','4'), ('5','6')])
        self.assertTrue(all(g.converged and g.rms < .008 for g in local.groups))
        self.assertEqual(local.unresolved, ('7',))
        self.assertGreater(local.groups[0].amplitude_ratio, 90.)
        recovered = [g.dominant/D_UNIT for g in local.groups]
        np.testing.assert_allclose(recovered, [3., 3.25, 9.], rtol=.05)
        np.testing.assert_allclose(np.array(baseline[:6])/D_UNIT, [3.,3.,3.25,3.25,9.,9.], rtol=.02)
        y, mask, sigma = fit_inputs(rows[:-1])
        whole = solve_global_ilt(dataset.b, y, GlobalParameters(grid_for(rows[:-1],1e-6)), include=mask, sigma=sigma)
        self.assertTrue(whole.converged)
        self.assertEqual(len(whole.components), 2)
        self.assertTrue(np.all(np.isfinite(whole.rms)))
        self.assertEqual([row['fit']['D_m2_s'] for row in rows], baseline)

    def test_overlap_contributes_to_two_local_groups_without_duplication(self):
        dataset, rows = self.fixture([([2.],[1.],0.), ([2.],[.5],0.), ([16.],[1.],0.),
                                      ([16.],[.7],0.), ([2.,16.],[1.,.9],0.)])
        result = adaptive_local_ilt(dataset, rows, 10., 10., 1e-6)
        self.assertEqual(result.unresolved, ())
        self.assertEqual(len(result.profiles['5']), 2)
        self.assertTrue(result.groups[-1].independent)
        self.assertEqual([entry[0] for entry in result.membership[('1','2')]], ['5'])
        self.assertEqual([entry[0] for entry in result.membership[('3','4')]], ['5'])
        start, stop = np.searchsorted(dataset.ppm, rows[4]['bounds'])
        total = sum(y for _,y in result.profiles['5'])
        np.testing.assert_allclose(total, dataset.spectrum[0,start:stop], atol=1e-8)
        np.testing.assert_allclose([d/D_UNIT for d,_ in result.profiles['5']], [2.,16.], rtol=.10)

    def test_group_specific_lambda_preserves_membership_and_mono_fits(self):
        dataset,rows = self.fixture([([3.],[1.],0.),([3.],[.5],0.),([9.],[1.],0.),([9.],[.8],0.)])
        original = [row['fit']['D_m2_s'] for row in rows]
        baseline = adaptive_local_ilt(dataset,rows,5.,10.,1e-6)
        override = {('1','2'):2e-5}
        inputs = local_sweep_inputs(dataset,rows,5.,10.,1e-6,(),override)
        result = adaptive_local_ilt(dataset,rows,5.,10.,1e-6,group_regularizations=override)
        self.assertEqual([g.ids for g in result.groups],[g.ids for g in baseline.groups])
        self.assertEqual([g.regularization for g in result.groups],[2e-5,1e-6])
        self.assertEqual([source.current for source in inputs.values()],[2e-5,1e-6])
        self.assertEqual([row['fit']['D_m2_s'] for row in rows],original)
