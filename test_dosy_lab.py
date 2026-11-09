"""Synthetic and actual-data tests of the classical DOSY application."""
from pathlib import Path
from copy import deepcopy
import csv,hashlib,json,struct,tempfile,unittest
import numpy as np
from jeol import RawData,read_jdf,manual_b,validate_b,factory_b,parameters,GAMMA_H,infer_filter_delay
from core import (default_options,process_fids,signed_integrals,integrate,fit_decay,analyze,
                  baseline_fit,save_profile,load_profile,initial_dataset)
from d_statistics import build_snapshot,descriptive_statistics,summarize,diagnostic_points,export_statistics,import_statistics
from exports import export_bundle,save_project,read_project
from examples.generate_demo import demo_bytes, encode_fixture

ROOT=Path(__file__).parent
DATA=ROOT/'data/DOSY_sacharoza'
STE=DATA/'Sacharoza_50mM_D2O_2D_DOSY-1-1.jdf'
DSTE=DATA/'Sacharoza_50mM_D2O_2D_DOSY_CC-1-1.jdf'
real_sucrose = unittest.skipUnless(STE.is_file() and DSTE.is_file(), 'Lokalne JDF sacharozy są opcjonalne.')


def synthetic(nlayers=16,npoints=2048):
    sw=1000.;freq=600e6;b=np.linspace(0,4e9,nlayers);t=np.arange(npoints)/sw
    trace=np.exp(-np.pi*2*t)*(np.exp(2j*np.pi*100*t)+.5*np.exp(-2j*np.pi*150*t))
    fid=np.exp(-b[:,None]*3.86e-10)*trace
    return RawData(fid.astype(complex),np.linspace(.05,.4,nlayers),b,
        dict(source_name='synthetic',source_path='synthetic.jdf',source_sha256='synthetic',sequence='synthetic',
             x_sweep_Hz=sw,x_sweep_clipped_Hz=sw,x_frequency_Hz=freq,offset_ppm=4.5,
             b_source={'kind':'manual','description':'known synthetic b'},parameters={}))


def fake_jdf(edge=4,little=True,dtype_code=0,sequence='bpp_ste_led_diffusion.jxp',shape='S_RECTANGLE'):
    """Encode a predictable FID with a padded final tile and known metadata units."""
    nx,ng=3*edge,2*edge;real=np.arange(nx*ng).reshape(ng,nx)+1.;imag=real*.2
    g=np.linspace(.05,.4,ng)
    values=dict(x_domain='Proton',x_sweep=1000.,x_freq=600e6,x_offset=4.5,delta=2.,delta_large=100.,
                tau=2.,tau_1=2.,tau_2=2.,grad_shape=shape)
    coef=.81*(.1-(2 if 'dste' in sequence else 1)*((6344*np.pi**2-207)/(19440*np.pi**2))*.002-(.002 if 'dste' in sequence else .001))
    values['st_coef']=coef*1000
    context=('experiment = "'+sequence+'";\ngrad_1 =? delta / 2;\ng => y_acq {'+', '.join(f'{v:.12g}[T/m]' for v in g)+'};\nbegin\nend pulse;')
    fid=(real-1j*imag)[:,:nx-1]
    blob=encode_fixture(fid,g,values,context,edge=edge,little=little,dtype_code=dtype_code)
    return blob,fid,g


class ReaderTests(unittest.TestCase):
    def test_generated_demo_fits_two_diffusion_coefficients(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'DOSY_demo.jdf'
            path.write_bytes(demo_bytes())
            dataset=initial_dataset(read_jdf(path))
        for bounds,expected in (([3.60,3.68],3.8e-10),([4.66,4.74],8e-10)):
            fit=fit_decay(dataset.b,integrate(dataset,bounds)['corrected_integral'])
            self.assertTrue(fit['ok'])
            self.assertAlmostEqual(fit['D_m2_s']/expected,1,delta=.06)

    def test_tiles_precision_endianness_axis_and_models(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'test.jdf'
            for edge in (4,32):
                for little in (True,False):
                    for dtype in (0,1):
                        blob,z,g=fake_jdf(edge,little,dtype);p.write_bytes(blob);r=read_jdf(p)
                        np.testing.assert_allclose(r.fid,z,rtol=1e-7);np.testing.assert_allclose(r.gradients,g,rtol=1e-12)
                        self.assertEqual(r.metadata['sequence'],'BPP-STE-LED');self.assertIsNotNone(r.b)
                        self.assertFalse(r.fid.flags.writeable)

    def test_unknown_sequence_needs_manual_b_and_does_not_guess(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'unknown.jdf';blob,_,_=fake_jdf(sequence='unknown_diffusion.jxp');p.write_bytes(blob)
            r=read_jdf(p);self.assertIsNone(r.b)
            self.assertFalse(fit_decay(r.b,np.ones(len(r.fid)))['ok'])
            new=manual_b(r,np.linspace(0,3e9,len(r.fid)),'independent experiment table')
            self.assertEqual(new.metadata['b_source']['kind'],'manual');self.assertIsNone(r.b)
            with self.assertRaises(ValueError):manual_b(r,[1,2],'source')
            with self.assertRaises(ValueError):manual_b(r,np.ones(len(r.fid)),'')

    def test_corrupt_file_and_gradient_mismatch_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.jdf';blob,_,_=fake_jdf();p.write_bytes(blob[:-1])
            with self.assertRaises(ValueError):read_jdf(p)
            blob,_,_=fake_jdf();at=struct.unpack_from('>I',blob,1224)[0];struct.pack_into('>d',blob,at,90.);p.write_bytes(blob)
            with self.assertRaises(ValueError):read_jdf(p)

    @real_sucrose
    def test_real_data_header_and_independent_tile_coordinates(self):
        for path,family in ((STE,'BPP-STE-LED'),(DSTE,'BPP-DSTE-LED')):
            original=hashlib.sha256(path.read_bytes()).hexdigest();r=read_jdf(path)
            self.assertEqual(r.fid.shape,(16,9567));self.assertEqual(r.metadata['sequence'],family)
            self.assertEqual(len(r.b),16);self.assertAlmostEqual(r.gradients[0],.05);self.assertAlmostEqual(r.gradients[-1],.4)
            self.assertTrue(all(np.diff(r.b)>0));self.assertGreater(r.metadata['x_frequency_Hz'],600e6)
            blob=path.read_bytes();ds=struct.unpack_from('>I',blob,1284)[0];nx=9568;edge=4
            for y,x in ((0,0),(1,17),(15,9566)):
                i=((y//edge)*(nx//edge)+x//edge)*edge**2+(y%edge)*edge+x%edge
                re=struct.unpack_from('<d',blob,ds+8*i)[0];im=struct.unpack_from('<d',blob,ds+8*(nx*16+i))[0]
                self.assertEqual(r.fid[y,x],re-1j*im)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),original)

    def test_wrong_factory_coefficient_is_not_silently_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'test.jdf'
            path.write_bytes(fake_jdf()[0])
            r=read_jdf(path)
        p=deepcopy(r.metadata['parameters']);p['st_coef']['value']*=2
        b,info=factory_b(p,r.metadata['context'],r.gradients)
        self.assertIsNone(b);self.assertIn('st_coef',info['description'])


class ProcessingTests(unittest.TestCase):
    def test_filter_delay_requires_supported_raw_metadata(self):
        p={key:{'value':value} for key,value in dict(digital_filter='TRUE',clipped='FALSE',
            digital_filter_status='Np',orders='2 54 73',factors='8 2',filter_factor=16).items()}
        self.assertEqual(infer_filter_delay(p,0,2048)['samples'],19.65625)
        self.assertIsNone(infer_filter_delay(p,1,2048))
        for key,value in [('orders','3 54 73'),('factors','4 2'),('clipped','TRUE'),
                          ('digital_filter','FALSE'),('digital_filter_status','unknown')]:
            with self.subTest(key=key):
                changed=dict(p,**{key:{'value':value}})
                self.assertIsNone(infer_filter_delay(changed,0,2048))

    @real_sucrose
    def test_open_sucrose_suppresses_ripple_without_baseline_or_changing_raw(self):
        for path in (STE,DSTE):
            with self.subTest(path=path):
                raw=read_jdf(path);saved=raw.fid.copy()
                before=process_fids(raw);after=initial_dataset(raw);x=after.ppm
                self.assertEqual(after.options['filter_delay_samples'],19.6875)
                self.assertEqual(after.options['phase1_deg'],0)
                self.assertFalse(after.options['baseline_enabled'])
                for lo,hi in [(2.65,2.95),(5.65,6.35)]:
                    mask=(x>=lo)&(x<=hi)
                    self.assertLess(np.std(after.phased[0,mask].real),.1*np.std(before.phased[0,mask].real))
                fit=fit_decay(after.b,integrate(after,[3.60,3.68])['corrected_integral'])
                self.assertTrue(fit['ok']);self.assertGreater(fit['R2'],.995)
                np.testing.assert_array_equal(raw.fid,saved)
                np.testing.assert_array_equal(after.b,before.b)
                with tempfile.TemporaryDirectory() as tmp:
                    project=Path(tmp)/'project.json'
                    save_project(project,after,[])
                    _,opts,_=read_project(project)
                    self.assertEqual(opts,after.options)
                    # A v1 project with a manual 7200° PH1 must not receive a second ramp.
                    old=json.loads(project.read_text(encoding='utf-8'))
                    old['processing'].pop('filter_delay_samples')
                    old['processing']['phase1_deg']=7200.
                    project.write_text(json.dumps(old),encoding='utf-8')
                    _,opts,_=read_project(project)
                    self.assertEqual(opts['filter_delay_samples'],0)
                    self.assertEqual(opts['phase1_deg'],7200)

    def test_known_FID_delay_uses_one_ramp_without_changing_amplitude_ratios(self):
        r=synthetic(npoints=4096);delay=20
        delayed=np.zeros_like(r.fid);delayed[:,delay:]=r.fid[:,:-delay]
        target=np.zeros_like(r.fid);target[:,:-delay]=r.fid[:,:-delay]
        raw=RawData(delayed,r.gradients.copy(),r.b.copy(),dict(r.metadata))
        saved=raw.fid.copy();o=default_options(raw);o['filter_delay_samples']=delay
        corrected=process_fids(raw,o)
        expected=np.fft.fftshift(np.fft.fft(target,axis=1),axes=1)
        np.testing.assert_allclose(corrected.phased,expected,atol=1e-10)
        areas=integrate(corrected,[4.60,4.73])['corrected_integral']
        np.testing.assert_allclose(areas/areas[0],np.exp(-r.b*3.86e-10),atol=1e-12)
        np.testing.assert_array_equal(raw.fid,saved)

    @unittest.skipUnless((ROOT/'data/Strychnina_CDCl3_bpp_ste_led_diffusion-1-1.jdf').exists(),'Optional strychnine reference JDF unavailable')
    def test_strychnine_manual_common_phase_reduces_periodic_blank_baseline(self):
        # Empirical diagnostic profile for this JDF, never an application default.
        raw=read_jdf(ROOT/'data/Strychnina_CDCl3_bpp_ste_led_diffusion-1-1.jdf')
        original=raw.fid.copy();o=default_options(raw);before=process_fids(raw,o)
        o['phase1_deg']=7200.;trial=process_fids(raw,o)
        o['phase0_deg']=-float(np.angle(signed_integrals(trial.ppm,trial.phased,[7.03,7.5])[0],deg=True))
        after=process_fids(raw,o);mask=(after.ppm>=10)&(after.ppm<=12.5)
        self.assertLess(np.std(after.spectrum[0,mask].real),.3*np.std(before.spectrum[0,mask].real))
        self.assertFalse(after.options['baseline_enabled']);np.testing.assert_array_equal(raw.fid,original)

    def test_one_common_phase_restores_absorption_and_all_integral_ratios(self):
        n=4096;sw=1200.;f=np.fft.fftshift(np.fft.fftfreq(n,1/sw));ppm=4.5+f/600
        ideal=10/(1+1j*(f-100)/2)+5/(1+1j*(f+180)/3)
        correction=35+240*(ppm-4.6)/(sw/600)
        amplitudes=np.exp(-np.linspace(0,4e9,16)*3.86e-10)
        distorted=amplitudes[:,None]*ideal*np.exp(-1j*np.deg2rad(correction))[None,:]
        original=synthetic(npoints=n)
        meta=dict(original.metadata,x_sweep_Hz=sw,x_sweep_clipped_Hz=sw)
        fid=np.fft.ifft(np.fft.ifftshift(distorted,axes=1),axis=1)
        raw=RawData(fid,np.linspace(.05,.4,16),np.linspace(0,4e9,16),meta)
        saved=raw.fid.copy();o=default_options(raw);o.update(phase0_deg=35,phase1_deg=240,pivot_ppm=4.6)
        d=process_fids(raw,o)
        np.testing.assert_allclose(d.spectrum,amplitudes[:,None]*ideal,atol=2e-13)
        self.assertGreater(d.spectrum.real.min(),0)
        areas=integrate(d,[4.55,4.80])['corrected_integral']
        np.testing.assert_allclose(areas/areas[0],amplitudes/amplitudes[0],atol=1e-13)
        np.testing.assert_array_equal(saved,raw.fid)

    def test_lowest_b_is_reference_without_reordering_layers(self):
        r=synthetic();order=np.r_[3,0,2,1,np.arange(4,16)]
        raw=RawData(r.fid[order].copy(),r.gradients[order].copy(),r.b[order].copy(),dict(r.metadata))
        self.assertEqual(default_options(raw)['reference_layer'],1)
        np.testing.assert_array_equal(raw.b,r.b[order])
        missing_b=RawData(r.fid.copy(),np.r_[-.4,.1,-.05,np.full(13,np.nan)],None,dict(r.metadata))
        self.assertEqual(default_options(missing_b)['reference_layer'],2)

    def test_fft_axis_and_common_phase_pivot(self):
        r=synthetic();before=r.fid.copy();o=default_options(r);a=process_fids(r,o)
        i=np.argmax(abs(a.phased[0]));self.assertAlmostEqual(a.ppm[i],4.5+100/600,delta=1/600)
        o.update(phase0_deg=37,phase1_deg=122,pivot_ppm=4.4)
        z=process_fids(r,o);phase=37+122*(a.ppm-4.4)/(1000/600)
        np.testing.assert_allclose(z.phased,a.phased*np.exp(1j*np.deg2rad(phase)),atol=1e-10)
        np.testing.assert_array_equal(before,r.fid)
        # Cropping the filter passband must not redefine PH1 over a smaller SW.
        meta=dict(r.metadata,digital_filter=True,already_clipped=False,x_sweep_clipped_Hz=600.)
        cropped=process_fids(RawData(r.fid.copy(),r.gradients.copy(),r.b.copy(),meta),o)
        keep=abs(a.ppm-4.5)<=.5
        np.testing.assert_array_equal(cropped.ppm,z.ppm[keep])
        np.testing.assert_allclose(cropped.phased,z.phased[:,keep],atol=1e-10)

    def test_exp_gaussian_and_combination_have_documented_formula(self):
        r=synthetic();o=default_options(r);o.update(lb_Hz=1.3,gaussian_FWHM_Hz=2.1,crop_filter=False)
        d=process_fids(r,o);time=np.arange(r.fid.shape[1])/1000
        window=np.exp(-np.pi*1.3*time-(np.pi*2.1*time)**2/(4*np.log(2)))
        expected=np.fft.fftshift(np.fft.fft(r.fid*window,axis=1),axes=1)
        np.testing.assert_allclose(d.spectrum,expected)

    def test_signed_integration_descending_axis_and_endpoint_interpolation(self):
        x=np.linspace(0,10,101);y=np.array([2*x+1,-2*x-1])
        forward=signed_integrals(x,y,[2.05,7.08]);reverse=signed_integrals(x[::-1],y[:,::-1],[7.08,2.05])
        np.testing.assert_allclose(forward,reverse)
        expected=7.08**2+7.08-(2.05**2+2.05)
        np.testing.assert_allclose(forward,[expected,-expected],atol=1e-12)

    def test_baseline_polynomial_mask_stability_and_undo(self):
        x=np.linspace(2,8,2001);base=.3+.2*(x-5)+.05*(x-5)**5
        signal=np.exp(-((x-5)/.03)**2)*10;y=np.stack([base+signal,2*base+signal])
        regions=[[2,4.7],[5.3,8]]
        z,info=baseline_fit(x,y,regions,5,per_layer=True)
        np.testing.assert_allclose(z,[base,2*base],atol=1e-10)
        self.assertLess(info['condition'],1e4)
        r=synthetic();o=default_options(r);original=process_fids(r,o)
        o.update(baseline_enabled=True,baseline_per_layer=True,baseline_regions=[[3.68,3.9],[5.1,5.3]],baseline_degree=2)
        corrected=process_fids(r,o);self.assertFalse(np.array_equal(corrected.spectrum,original.spectrum))
        o['baseline_enabled']=False;undo=process_fids(r,o);np.testing.assert_array_equal(undo.spectrum,original.spectrum)
        o['baseline_enabled']=True
        with self.assertRaisesRegex(ValueError,'przecina'):process_fids(r,o,[[3.8,3.9]])

    def test_processing_profile_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            o=default_options(synthetic());o['phase1_deg']=7200;p=Path(d)/'profile.json';save_profile(p,o)
            self.assertEqual(load_profile(p),o)


class FitStatisticsExportTests(unittest.TestCase):
    def test_known_diffusion_all_active_and_explicit_exclusion(self):
        b=np.linspace(0,4e9,16);y=50*np.exp(-b*3.86e-10)
        fit=fit_decay(b,y);self.assertTrue(fit['ok']);self.assertAlmostEqual(fit['D_m2_s']/3.86e-10,1.,places=8)
        corrupt=y.copy();corrupt[5]*=1.8;bad=fit_decay(b,corrupt);mask=np.ones(16,bool);mask[5]=False
        good=fit_decay(b,corrupt,mask);self.assertGreater(abs(bad['D_m2_s']/3.86e-10-1),.02)
        self.assertAlmostEqual(good['D_m2_s']/3.86e-10,1.,places=7);self.assertEqual(len(good['prediction']),16)
        self.assertFalse(good['fit_mask'][5])

    def test_negative_integrals_nan_and_small_n(self):
        b=np.arange(6)*1e9;y=np.exp(-b*3e-10);y[-1]=-.1
        fit=fit_decay(b,y);self.assertTrue(fit['ok']);self.assertTrue(any('niedodatnie' in w for w in fit['warnings']))
        self.assertFalse(fit_decay(b,-np.ones(6))['ok']);self.assertFalse(fit_decay(b[:2],y[:2])['ok'])
        y[1]=np.nan;self.assertFalse(fit_decay(b,y)['ok'])

    def test_weighted_uncertainty_and_noise_SE(self):
        rng=np.random.default_rng(3);b=np.linspace(0,5e9,32);sigma=np.full(32,.002)
        y=np.exp(-b*3.86e-10)+rng.normal(0,sigma)
        f=fit_decay(b,y,sigma=sigma);self.assertTrue(f['weighted']);self.assertIsNotNone(f['D_CI95_m2_s'])
        self.assertLess(abs(f['D_m2_s']/3.86e-10-1),.02)
        self.assertFalse(fit_decay(b,y,sigma=np.zeros(32))['ok'])

    def test_statistics_do_not_mix_compounds_and_preserve_N01(self):
        s={'records':[dict(id='A',group='compound1',eligible=True,D=3e-10,SE=1e-12),dict(id='B',group='compound2',eligible=True,D=9e-10,SE=1e-12)]}
        self.assertEqual(summarize(s,{})['N'],0)
        self.assertEqual(summarize(s,{},group='compound1')['mean'],3e-10)
        self.assertEqual(summarize(s,{'A':False},group='compound1')['N'],0)
        self.assertIsNone(descriptive_statistics([3e-10])['SD'])
        z=descriptive_statistics([2e-10,4e-10]);self.assertAlmostEqual(z['SD']/1e-10,np.sqrt(2))

    def test_export_csv_xlsx_snapshots_and_complex_arrays(self):
        from openpyxl import load_workbook
        r=synthetic();o=default_options(r);o.update(phase0_deg=12.,phase1_deg=30.,pivot_ppm=4.6)
        d=process_fids(r,o);entries=[dict(id='S1',name='signal',group='compound',bounds=[4.61,4.72],color='#67e8f9',fit_mask=[True]*16,sigma=[None]*16,weighted=False)]
        rows=analyze(d,entries)
        with tempfile.TemporaryDirectory() as tmp:
            out=export_bundle(Path(tmp)/'export',d,rows)
            self.assertTrue((out/'Raport.md').exists());self.assertTrue((out/'processing.json').exists())
            saved=json.loads((out/'processing.json').read_text(encoding='utf-8'))
            self.assertEqual([saved[k] for k in ('phase0_deg','phase1_deg','pivot_ppm')],[12.,30.,4.6])
            book=load_workbook(out/'DOSY_results.xlsx',read_only=True,data_only=True);self.assertEqual(book['ROI']['A2'].value,'S1');book.close()
            with np.load(out/'complex_data.npz') as data:np.testing.assert_array_equal(data['fid'],r.fid)
            snap,choices=import_statistics(out/'statistics_snapshot.json');self.assertEqual(snap['records'][0]['group'],'compound')
            pts,_=diagnostic_points(snap,snap['records'][0]);self.assertEqual(len(pts),16)

    @real_sucrose
    def test_actual_project_roundtrip_and_real_regression(self):
        for path in (STE,DSTE):
            r=read_jdf(path);o=default_options(r)
            # Explicit validation profile: phase from the reference spectrum only,
            # not from a fit to D. The observed delayed FID starts near sample 20.
            delay=int(np.argmax(abs(r.fid[0,:100])))
            o.update(phase1_deg=360.*delay,lb_Hz=.3,zero_fill_factor=2)
            first=process_fids(r,o)
            o['phase0_deg']=-float(np.angle(signed_integrals(first.ppm,first.phased,[5.26,5.33])[0],deg=True))
            d=process_fids(r,o);entry=dict(id='S1',name='3.60–3.68 ppm',group='sucrose',bounds=[3.60,3.68],fit_mask=[True]*16,sigma=[None]*16,color='#67e8f9')
            rows=analyze(d,[entry]);f=rows[0]['fit'];self.assertTrue(f['ok']);self.assertGreater(f['R2'],.995)
            self.assertGreater(f['D_m2_s'],3.7e-10);self.assertLess(f['D_m2_s'],4.05e-10)
            with tempfile.TemporaryDirectory() as tmp:
                p=Path(tmp)/'project.json';save_project(p,d,[entry]);rr,oo,ee=read_project(p)
                self.assertEqual(oo,o);self.assertEqual(ee[0]['bounds'],entry['bounds']);np.testing.assert_array_equal(rr.fid,r.fid)


if __name__=='__main__':unittest.main()
