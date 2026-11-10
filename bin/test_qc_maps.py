"""Known-truth controls for the 2.1 QC and map separation."""
from types import SimpleNamespace
import unittest

import numpy as np
from scipy.integrate import trapezoid

from core import fit_decay
from d_statistics import quality
from dosy_alignment import ilt_profiles, profile_residual
from dosy_map import roi_bands
from roi_quality import qualify_roi


class QCMapTests(unittest.TestCase):
    def test_low_information_layers_are_diagnostic_only(self):
        b = np.linspace(0,14,24)*1e9
        x = np.linspace(-1,1,81)
        decay = np.exp(-b*3.8e-10)
        profile = 1/(1+1j*(x-.07)/.18)
        rng = np.random.default_rng(2102)
        noise = .04*np.linspace(.5,1.5,len(b))
        for phase in (0.,25.):
            spectrum = decay[:,None]*profile*np.exp(1j*np.deg2rad(phase))
            spectrum += noise[:,None]*(rng.normal(size=spectrum.shape)+1j*rng.normal(size=spectrum.shape))
            y = trapezoid(spectrum.real,x=x,axis=1)
            original = y.copy()
            fit = fit_decay(b,y)
            fit_mask = fit['fit_mask'].copy()
            data = SimpleNamespace(b=b,ppm=x,spectrum=spectrum,options={'baseline_enabled':False})
            diag, findings = qualify_roi(data,(-1.,1.),fit,dict(corrected_integral=y),None,noise,fit_decay)
            self.assertEqual(diag['status'],'insufficient')
            self.assertLess(len(diag['diagnostic_layers']),8)
            self.assertTrue(all(f['severity']=='info' for f in findings))
            np.testing.assert_array_equal(y,original)
            np.testing.assert_array_equal(fit['fit_mask'],fit_mask)
            self.assertEqual(quality(dict(fit,measurement_findings=findings))[1],'INFO')

    def test_residual_noise_energy_is_separate_from_assignment(self):
        target = np.ones((8,10))
        predicted = target.copy()
        predicted[:2] += .5
        noiseless = profile_residual(target,predicted)
        adjusted = profile_residual(target,predicted,noise_energy=5.)
        self.assertGreater(noiseless['excess_residual'],.15)
        self.assertEqual(noiseless['relative_residual'],adjusted['relative_residual'])
        self.assertEqual(adjusted['excess_residual'],0.)

    def test_single_maximum_must_explain_roi_even_with_exact_interval(self):
        b = np.linspace(0,9,22)*1e9
        ppm = np.linspace(1,2,121)
        profile = np.exp(-.5*((ppm-1.5)/.07)**2)
        spectra = np.exp(-b*1.9e-9)[:,None]*profile
        y = trapezoid(spectra[:,:-1],x=ppm[:-1],axis=1)
        row = dict(id='R1',bounds=(1.,2.),result={'corrected_integral':y},fit=fit_decay(b,y))
        # The interval is exact for this ROI, but its maximum across the joint
        # result is at a different D, with no weight there in this column.
        solution = SimpleNamespace(converged=True,diffusion=np.array([1e-10,3.8e-10,1.9e-9,3e-9]),
            weights=np.array([[0.],[0.],[1.],[0.]]), intervals=((0,4),),components=(1,),
            scales=np.array([max(y)]),rms=np.array([0.]))
        details = {}
        profiles,rejected = ilt_profiles(b,ppm,spectra,0,[row],('R1',),solution,diagnostics=details)
        self.assertFalse(profiles)
        self.assertEqual(rejected,('R1',))
        self.assertEqual(details['R1']['reason'],'peak_mismatch')
        self.assertGreater(details['R1']['peak_rms'],details['R1']['peak_limit'])

    def test_unresolved_overlap_cannot_be_attributed_to_an_assigned_roi(self):
        ppm = np.arange(9.)
        reference = np.ones_like(ppm)
        fit = {'ok': True, 'D_m2_s': 2e-10}
        rows = [dict(id='A', bounds=(1., 6.), fit=fit), dict(id='B', bounds=(4., 8.), fit=fit)]
        bands, overlap = roi_bands(ppm,reference,rows,profiles={'A': ((2e-10,np.ones(5)),)},
                                   reconstructed_only=True)
        self.assertTrue(overlap)
        self.assertEqual([b.sid for b in bands],['A'])
        np.testing.assert_array_equal(bands[0].intensity,[1.,1.,1.,0.,0.])

    def test_partial_multicomponent_map_does_not_renormalize_missing_mass(self):
        b = np.r_[0.,np.geomspace(.02,35.,48)]*1e9
        ppm = np.linspace(0,1,201)
        d = np.array([2.,6.,19.])*1e-10
        a = np.array([1.,.5,.045])
        shape = np.exp(-.5*((ppm[None,:]-np.array([.2,.5,.8])[:,None])/.035)**2)*a[:,None]
        spectra = np.exp(-np.outer(b,d))@shape
        y = trapezoid(spectra[:,:-1],x=ppm[:-1],axis=1)
        row = dict(id='R1',bounds=(0.,1.),result={'corrected_integral':y},fit=fit_decay(b,y))
        sol = SimpleNamespace(converged=True,diffusion=np.r_[1e-10,d,3e-9],
            weights=np.r_[0.,a,0.][:,None],intervals=((0,2),(2,3),(3,5)),components=(1,2,3),
            scales=np.array([max(y)]),rms=np.array([0.]))
        diagnostics = {}
        profiles,rejected = ilt_profiles(b,ppm,spectra,0,[row],('R1',),sol,
                                         diagnostics=diagnostics,conserve_reference=True)
        self.assertFalse(rejected)
        self.assertEqual(len(profiles['R1']),2)
        self.assertAlmostEqual(diagnostics['R1']['omitted_fraction'],a[-1]/sum(a))
        self.assertGreater(diagnostics['R1']['unexplained_positive_fraction'],.02)
        self.assertLess(sum(v[160] for _,v in profiles['R1']),spectra[0,160]*.5)
