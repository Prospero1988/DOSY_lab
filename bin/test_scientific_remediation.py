"""Regression fixtures derived from scientific audit 42895f1 and benchmark.json."""
from copy import deepcopy
import json
from pathlib import Path
import struct
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np

from core import analyze, clean_json, default_options, fit_decay, integrate, process_fids
from d_statistics import quality
from dosy_alignment import ilt_profiles
from dosy_map import roi_bands
from exports import project_payload, read_project
from ilt import ILTParameters, sensitivity_trials, solve_ilt
from jeol import RawData, factory_b, manual_b, read_jdf, resolve_b_sources

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples'
SUCROSE = EXAMPLES / 'DOSY_sacharoza' / 'Sacharoza_50mM_D2O_2D_DOSY-1-1.jdf'
UNIT = 1e-10


def audit_fid(components, phase_slope=0., baseline=0., noise=0.):
    # audit/run_audit.py synthetic_raw and processing_cases at commit 42895f1
    b9 = np.linspace(0, 14, 24)
    t = np.arange(8192)/8000.
    fid = np.zeros((len(b9), len(t)), complex)
    for ppm, rate, amplitude, width, phase in components:
        trace = np.exp(2j*np.pi*(ppm-4.5)*600*t-np.pi*width*t)
        rotation = np.exp(1j*np.deg2rad(phase+phase_slope*b9/b9.max()))
        fid += (amplitude*np.exp(-b9*rate/10)*rotation)[:, None]*trace
    fid[:,0] += baseline*(1+1.5*b9/b9.max())
    rng = np.random.default_rng(20261114)
    if noise:
        fid += noise*(rng.normal(size=fid.shape)+1j*rng.normal(size=fid.shape))
    return RawData(fid,np.full(len(b9),np.nan),b9*1e9,
                   dict(x_sweep_Hz=8000.,x_frequency_Hz=600e6,offset_ppm=4.5,
                        source_name='analytic_FID',digital_filter=False))


class ScientificRemediationTests(unittest.TestCase):
    def test_b1_conflicting_timing_withholds_b_but_loads_and_recovers_manual_project(self):
        good = read_jdf(SUCROSE)
        data = bytearray(SUCROSE.read_bytes())
        start, length = struct.unpack_from('>2I', data, 1212)
        record = next(i for i in range(start+16,start+length,64)
                      if data[i+36:i+64].strip() == b'diffusion_time')
        # Real example is little endian, with timing stored as double in ms.
        struct.pack_into('<d', data, record+16, 200.)
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            bad_path = root/'inconsistent.jdf'
            bad_path.write_bytes(data)
            bad = read_jdf(bad_path)
            self.assertIsNone(bad.b)
            self.assertEqual(bad.metadata['b_source']['kind'],'conflict')
            self.assertIn('diffusion_time',bad.metadata['b_source']['description'])
            np.testing.assert_array_equal(bad.fid,good.fid)
            self.assertFalse(bad.fid.flags.writeable)
            packet = project_payload(process_fids(good),[])
            packet.update(source_path=str(bad_path),source_sha256=bad.metadata['source_sha256'])
            saved = root/'previous.proj'
            saved.write_text(json.dumps(clean_json(packet)),encoding='utf-8')
            self.assertIsNone(read_project(saved)[0].b)
            overridden = manual_b(bad,good.b,'independently supplied values for the fixture')
            self.assertEqual(overridden.metadata['b_history'][-1]['kind'],'conflict')
            np.testing.assert_array_equal(overridden.b,good.b)

    def test_b1_axis_conflict_and_nominal_b_match_reference_measurements(self):
        expected = {
            'Sacharoza_50mM_D2O_2D_DOSY-1-1.jdf': (.05,.4,.07966308112375117),
            'Sacharoza_50mM_D2O_2D_DOSY_CC-1-1.jdf': (.05,.4,.07913616224750233),
            'Strychnina_CDCl3_bpp_ste_led_diffusion-1-1.jdf': (.05,.2,.06346308112375117),
        }
        gamma = 2*np.pi*42.577478518e6
        for name,(low, high, coeff) in expected.items():
            raw = read_jdf(next(EXAMPLES.rglob(name)))
            np.testing.assert_allclose(raw.gradients[[0,-1]],[low,high],rtol=1e-10)
            np.testing.assert_allclose(raw.b,(gamma*raw.gradients*.002)**2*coeff,rtol=1e-9)
            accepted, source = resolve_b_sources(raw.b,raw.b,raw.metadata['b_source'])
            self.assertTrue(source['cross_checked'])
            np.testing.assert_array_equal(accepted,raw.b)
            conflicting, source = resolve_b_sources(raw.b*1.1,raw.b,raw.metadata['b_source'])
            self.assertIsNone(conflicting)
            self.assertEqual(source['kind'],'conflict')
            params = deepcopy(raw.metadata['parameters'])
            params['large_t']['value'] *= 2
            _, source = factory_b(params,raw.metadata['context'],raw.gradients)
            self.assertEqual(source['kind'],'conflict')

    def test_m1_late_layers_and_heteroscedastic_reference(self):
        b = np.linspace(100,180,12)*1e9
        y = np.exp(-b*.4*UNIT)
        fit = fit_decay(b,y)
        self.assertTrue(fit['ok'])
        self.assertAlmostEqual(fit['D_m2_s']/UNIT,.4,places=8)
        self.assertAlmostEqual(fit['amplitude'],1.,places=8)
        self.assertEqual(fit['amplitude_reference_b_s_m2'],b[0])
        np.testing.assert_allclose(fit['prediction'],y,rtol=1e-10)
        far = np.linspace(1e13,1.001e13,12)
        fit = fit_decay(far,np.exp(-(far-far[0])*UNIT))
        self.assertAlmostEqual(fit['D_m2_s']/UNIT,1.,places=8)
        self.assertIsNone(fit['amplitude'])
        self.assertAlmostEqual(fit['amplitude_at_reference_b'],1.,places=8)
        # F_heteroscedastic: four E fixtures consume these seeded draws first.
        b9 = np.r_[0.,np.geomspace(.02,35.,48)]
        rng = np.random.default_rng(20261110)
        rng.normal(0,.001,(4,len(b9)))
        sigma = np.linspace(.002,.02,len(b9))
        y = np.exp(-b9*.4)+rng.normal(0,sigma)
        fit = fit_decay(b9*1e9,y,sigma=sigma)
        self.assertAlmostEqual(fit['D_m2_s']/UNIT,3.9842724077070604,places=5)
        scaled = fit_decay(b9*1e9,y*1e-150,sigma=sigma*1e-150)
        self.assertAlmostEqual(fit['D_m2_s']/UNIT,scaled['D_m2_s']/UNIT,places=8)
        self.assertAlmostEqual(fit['R2'],scaled['R2'],places=12)

    def test_p1_and_p2_flags_do_not_change_measurements_or_fit(self):
        analyte = [(5.25,3.8,1.,3.,0.)]
        cases = (
            (analyte, {}, {}, set()),
            ([(5.25,3.8,1.,3.,25.)], {}, {}, set()),
            (analyte, {'noise':.04}, {}, set()),
            (analyte, {'phase_slope':55.}, {}, {'phase'}),
            (analyte+[(5.258,18.8,.35,4.,0.)], {}, {}, {'shape','phase'}),
            (analyte+[(5.01,19.,90.,12.,14.)], {}, {}, {'shape','phase'}),
            (analyte, {'baseline':8.}, {'baseline_enabled':True,'baseline_regions':[(0,1),(9,10)]}, {'baseline'}),
            (analyte+[(4.94,19.,60.,9.,12.)], {'baseline':8.},
             {'baseline_enabled':True,'baseline_regions':[(4.98,5.08),(5.55,5.8)]}, {'shape','phase'}),
        )
        bounds = (5.18,5.32)
        for components, effects, processing, flags in cases:
            with self.subTest(components=components,effects=effects,processing=processing):
                raw = audit_fid(components,**effects)
                fid = raw.fid.copy()
                dataset = process_fids(raw,dict(default_options(raw),**processing))
                existing = integrate(dataset,bounds)
                original_fit = fit_decay(dataset.b,existing['corrected_integral'])
                row, = analyze(dataset,[dict(id='R1',bounds=bounds)])
                np.testing.assert_array_equal(row['result']['corrected_integral'],existing['corrected_integral'])
                np.testing.assert_array_equal(row['result']['raw_integral'],existing['raw_integral'])
                self.assertEqual(row['fit'].get('D_m2_s'),original_fit.get('D_m2_s'))
                np.testing.assert_array_equal(raw.fid,fid)
                diagnostic = row['result']['profile_diagnostics']
                self.assertEqual(set(diagnostic['flagged']),flags)
                if processing.get('baseline_enabled'):
                    expected = np.linalg.norm(existing['corrected_integral']-existing['raw_integral'])
                    expected /= np.linalg.norm(existing['raw_integral'])
                    self.assertAlmostEqual(diagnostic['baseline_relative_rms'],expected)
                if flags and row['fit']['ok']:
                    self.assertEqual(quality(row['fit'])[1],'UWAGA')
                if 'phase_slope' in effects:
                    self.assertEqual(row['fit']['interpretation'],'apparent')

    def test_a2_weak_oracle_remains_unassigned_and_resolved_oracle_preserves_profiles(self):
        b = np.r_[0.,np.geomspace(.02,35.,48)]*1e9
        ppm = np.linspace(.65,1.35,169)
        bounds = (.69,1.31)
        start,stop = np.searchsorted(ppm,bounds)
        for fraction, expected in ((.03,0),(.55,2),(0.,1)):
            with self.subTest(fraction=fraction):
                shapes = np.exp(-.5*((ppm[None,:]-np.array([.95,1.05])[:,None])/.036)**2)
                shapes *= np.array([1.,fraction])[:,None]
                trace = np.exp(-np.outer(b,np.array([3.,12.])*UNIT))@shapes
                solution = SimpleNamespace(converged=True,diffusion=np.array([1.,3.,6.,12.,20.])*UNIT,
                    weights=np.array([0.,1.,0.,fraction,0.])[:,None],
                    intervals=((0,3),(3,5)),components=(1,3),scales=np.array([1.]),rms=np.array([0.]))
                rows = [dict(id='R1',bounds=bounds,fit=dict(ok=True,D_m2_s=3*UNIT,
                             fit_mask=np.ones(len(b),bool),RMSE=0.))]
                withheld = {}
                profiles,rejected = ilt_profiles(b,ppm,trace,0,rows,('R1',),solution,
                                                 conserve_reference=True,withheld=withheld)
                self.assertEqual(len(profiles.get('R1',())),expected)
                if expected:
                    self.assertFalse(rejected)
                    self.assertFalse(withheld)
                    np.testing.assert_allclose(sum(z for _,z in profiles['R1']),trace[0,start:stop],atol=1e-12)
                else:
                    self.assertEqual(rejected,('R1',))
                    self.assertAlmostEqual(withheld['R1']['omitted_fraction'],.03/1.03)
                    self.assertTrue(roi_bands(ppm,trace[0],rows)[0])
                    self.assertFalse(roi_bands(ppm,trace[0],rows,excluded=withheld,profiles=profiles)[0])

    def test_i1_sensitivity_compares_si_grid_peaks_without_changing_solution(self):
        b = np.r_[0.,np.geomspace(.02,35.,48)]*1e9
        y = np.exp(-np.outer(b,np.array([2.,16.])*UNIT))@np.array([.6,.4])
        p = ILTParameters(.25*UNIT,40*UNIT,points=100,regularization=1e-5)
        baseline = solve_ilt(b,y,p)
        weights = baseline.weights.copy()
        trials = sensitivity_trials(b,y,p,baseline)
        self.assertEqual([trial['kind'] for trial in trials],['current','regularization','regularization','spacing'])
        self.assertEqual(trials[-1]['parameters'].spacing,'linear')
        np.testing.assert_allclose(np.asarray(trials[0]['peaks'])/UNIT,[2.054233,16.456558],rtol=1e-6)
        np.testing.assert_allclose(np.asarray(trials[-1]['peaks'])/UNIT,[2.83375,22.70875],rtol=1e-6)
        self.assertTrue(all(trial['mass_interval'] for trial in trials))
        np.testing.assert_array_equal(baseline.weights,weights)
