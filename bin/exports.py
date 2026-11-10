"""Reproducible projects, CSV/XLSX tables and a compact analysis report."""
from pathlib import Path, PureWindowsPath
from datetime import datetime,timezone
import csv,json,uuid
from sys import float_info
import numpy as np
from core import clean_json,validate_options
from d_statistics import build_snapshot,groups,summarize,export_statistics,SCIENTIFIC_NOTE
from jeol import read_jdf,manual_b
from storage import DATA, EXAMPLES, WORK, require_extension


def project_payload(dataset,entries,*,integral_scale=1.):
    definitions=[{k:v for k,v in e.items() if k not in ('fit','result','error')} for e in entries]
    return dict(schema='dosy_lab_project',version=1,
                source_path=dataset.metadata['source_path'],source_sha256=dataset.metadata['source_sha256'],
                processing=dataset.options,rois=definitions,b=dataset.b,b_source=dataset.metadata['b_source'],
                integral_scale=integral_scale)


def save_project(path,dataset,entries,*,integral_scale=1.):
    packet=project_payload(dataset,entries,integral_scale=integral_scale)
    require_extension(path,'.proj').write_text(json.dumps(clean_json(packet),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


def integral_definitions(rows):
    if not isinstance(rows, list) or not rows:
        raise ValueError('Lista całek jest pusta lub nieprawidłowa.')
    seen = set()
    definitions = []
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError('Nieprawidłowa definicja całki.')
        sid, group, bounds = row.get('id'), row.get('group'), row.get('bounds')
        if not isinstance(sid, str) or not sid.strip() or not isinstance(group, str):
            raise ValueError('Nieprawidłowe ID lub grupa całki.')
        if sid in seen:
            raise ValueError('Powtórzone ID całki.')
        if (not isinstance(bounds, (list, tuple)) or len(bounds) != 2
                or not all(type(value) in (int, float) and abs(value) <= float_info.max for value in bounds)
                or bounds[0] >= bounds[1]):
            raise ValueError('Nieprawidłowe granice całki.')
        seen.add(sid)
        definitions.append(dict(id=sid, group=group, bounds=list(bounds)))
    return definitions


def save_integrals(path, entries):
    packet = dict(schema='dosy_lab_integrals', version=1, integrals=integral_definitions(entries))
    require_extension(path, '.ints').write_text(json.dumps(packet, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def read_integrals(path):
    packet = json.loads(require_extension(path, '.ints').read_text(encoding='utf-8'))
    if not isinstance(packet, dict) or packet.get('schema') != 'dosy_lab_integrals' or packet.get('version') != 1:
        raise ValueError('Nieobsługiwany plik całek DOSY Lab.')
    return integral_definitions(packet.get('integrals'))

def _project_source(project, packet):
    original = packet['source_path']
    name = PureWindowsPath(original).name
    saved = Path(original)
    if not saved.is_absolute():
        saved = project.parent / saved
    directories = (DATA, WORK, EXAMPLES)
    def candidates():
        if not original.startswith(('\\\\', '//')):
            yield saved
        yield project.parent / name
        for directory in directories:
            yield directory / name
        for directory in directories:
            for candidate in directory.rglob('*'):
                if candidate.name.casefold() == name.casefold():
                    yield candidate
        yield saved

    seen = set()
    mismatch = False
    read_error = None
    for candidate in candidates():
        if candidate in seen:
            continue
        seen.add(candidate)
        try:
            raw = read_jdf(candidate)
        except FileNotFoundError:
            continue
        except (OSError, ValueError) as exc:
            read_error = exc
            continue
        if raw.metadata['source_sha256'] == packet['source_sha256']:
            return raw
        mismatch = True
    if mismatch:
        raise ValueError('SHA źródłowego JDF różni się od zapisanego projektu.')
    if read_error is not None:
        raise ValueError(f'Nie udało się odczytać źródłowego JDF {name}.') from read_error
    raise FileNotFoundError(f'Nie znaleziono źródłowego JDF {name}. Umieść go w katalogu data albo przy pliku projektu.')


def read_project(path,*,with_integral_scale=False):
    path=require_extension(path,'.proj').resolve()
    packet=json.loads(path.read_text(encoding='utf-8'))
    if packet.get('schema')!='dosy_lab_project' or packet.get('version')!=1:raise ValueError('Nieobsługiwany projekt DOSY Lab.')
    raw=_project_source(path,packet)
    if packet['b_source'].get('kind')=='manual':raw=manual_b(raw,packet['b'],packet['b_source']['description'])
    elif packet['b'] is not None:
        if raw.b is None or not np.allclose(raw.b,packet['b'],rtol=1e-12,atol=0):raise ValueError('Lista b projektu różni się od ponownie odczytanego JDF.')
    options=validate_options(packet['processing'],len(raw.fid));seen=set()
    for e in packet['rois']:
        if e['id'] in seen:raise ValueError('Powtórzone ID ROI w projekcie.')
        seen.add(e['id'])
        if len(e['bounds'])!=2 or not np.isfinite(e['bounds']).all():raise ValueError('Niepoprawny ROI w projekcie.')
        if len(e.get('fit_mask',[]))!=len(raw.fid):raise ValueError('Maska fitu nie pasuje do warstw.')
    if with_integral_scale:
        scale=float(packet.get('integral_scale',1.))
        if not np.isfinite(scale) or scale<=0:
            raise ValueError('Nieprawidłowa skala całek w projekcie.')
        return raw,options,packet['rois'],scale
    return raw,options,packet['rois']


def _table(path,header,rows):
    with path.open('w',newline='',encoding='utf-8-sig') as h:
        w=csv.writer(h);w.writerow(header);w.writerows(rows)


def export_bundle(directory,dataset,entries,choices=None,*,integral_scale=1.):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=False)
    choices=dict(choices or {})
    snapshot=build_snapshot(dataset,entries,uuid.uuid4().hex,dataset.name)
    roi_header=['ID','name','compound','low_ppm','high_ppm','D_m2_s','D_SE_m2_s','CI95_low','CI95_high','I0','R2','RMSE','N_active','QC','notes']
    roi_rows=[];point_rows=[]
    for e,r in zip(entries,snapshot['records']):
        f=e.get('fit') or {};ci=f.get('D_CI95_m2_s') or [None,None]
        roi_rows.append([e['id'],e.get('name',''),e.get('group',''),*sorted(e['bounds']),r['D'],r['SE'],*ci,
                         f.get('amplitude'),r['R2'],f.get('RMSE'),f.get('n_active'),r['qc'],r['qc_detail']])
        result=e.get('result') or {};pre=result.get('raw_integral');y=result.get('corrected_integral')
        if y is None:continue
        mask=f.get('fit_mask',e.get('fit_mask',np.ones(len(y),bool)))
        pred=f.get('prediction',[None]*len(y));res=f.get('residuals',[None]*len(y))
        sigma=e.get('sigma')
        if sigma is None:sigma=[None]*len(y)
        for i in range(len(y)):
            point_rows.append([e['id'],i+1,float(dataset.gradients[i]),None if dataset.b is None else float(dataset.b[i]),
                               float(pre[i]),float(y[i]),float(pred[i]) if pred[i] is not None else None,
                               float(res[i]) if res[i] is not None else None,bool(mask[i]),sigma[i]])
    point_header=['ID','acquired_layer','g_T_m','b_s_m2','pre_baseline_integral','processed_signed_integral','prediction','residual','fit_included','sigma_I_user']
    _table(directory/'roi_results.csv',roi_header,roi_rows)
    _table(directory/'integrals_and_fits.csv',point_header,point_rows)
    group_rows=[]
    for group in groups(snapshot):
        s=summarize(snapshot,choices,group=group)
        group_rows.append([group,s['N'],s['mean'],s['median'],s['SD'],s['CV_percent'],s['SEM'],s['minimum'],s['maximum']])
    group_header=['compound','N','mean_D_m2_s','median_D_m2_s','SD_sample_m2_s','CV_percent','SEM_descriptive_m2_s','min_D_m2_s','max_D_m2_s']
    _table(directory/'compound_statistics.csv',group_header,group_rows)
    metadata=dict(dataset.metadata);metadata.pop('context',None)
    (directory/'experiment.json').write_text(json.dumps(clean_json(metadata),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    (directory/'processing.json').write_text(json.dumps(clean_json(dataset.options),ensure_ascii=False,indent=2),encoding='utf-8')
    (directory/'fit_results.json').write_text(json.dumps(clean_json(entries),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    (directory/'embedded_sequence.txt').write_text(dataset.metadata.get('context',''),encoding='utf-8')
    save_project(directory/'project.proj',dataset,entries,integral_scale=integral_scale)
    np.savez_compressed(directory/'complex_data.npz',fid=dataset.raw.fid,ppm=dataset.ppm,
                        phased=dataset.phased,spectrum=dataset.spectrum,baseline=dataset.baseline,
                        gradients=dataset.gradients,b=np.array([]) if dataset.b is None else dataset.b)
    export_statistics(directory,snapshot,choices)
    from openpyxl import Workbook
    from openpyxl.styles import Font,PatternFill
    book=Workbook();book.remove(book.active)
    options=[(k,json.dumps(v,ensure_ascii=False)) for k,v in dataset.options.items()]
    exp=[('source',dataset.name),('SHA256',dataset.metadata['source_sha256']),
         ('sequence',dataset.metadata['sequence']),('b_source',json.dumps(dataset.metadata['b_source'],ensure_ascii=False)),
         ('units','g T/m; b s/m2; D m2/s; integral a.u. ppm')]
    for title,header,rows in [('ROI',roi_header,roi_rows),('Punkty',point_header,point_rows),
                              ('Statystyka',group_header,group_rows),('Processing',['parameter','value'],options),
                              ('Eksperyment',['parameter','value'],exp)]:
        sheet=book.create_sheet(title);sheet.append(header)
        for row in rows:
            sheet.append(clean_json(row))
            for cell in sheet[sheet.max_row]:
                if isinstance(cell.value,str):cell.data_type='s'
        for cell in sheet[1]:cell.fill=PatternFill('solid',fgColor='26364E');cell.font=Font(color='FFFFFF',bold=True)
        sheet.freeze_panes='A2';sheet.auto_filter.ref=sheet.dimensions
        for col in sheet.columns:sheet.column_dimensions[col[0].column_letter].width=min(65,max(14,max(len(str(c.value or '')) for c in col)+2))
    book.save(directory/'DOSY_results.xlsx')
    report=['# DOSY Lab — raport analizy','',dataset.name,'',
            'Sekwencja: '+str(dataset.metadata['sequence']),
            'Źródło b: '+dataset.metadata['b_source']['description'],
            'SHA256: '+dataset.metadata['source_sha256'],'',
            'Processing: `'+json.dumps(dataset.options,ensure_ascii=False)+'`','',
            '| ROI | Grupa | D [m²/s] | SE(D) | R² | QC |','|---|---|---:|---:|---:|---|']
    for r in snapshot['records']:report.append(f"| {r['id']} | {r['group']} | {r['D']} | {r['SE']} | {r['R2']} | {r['qc']} |")
    report+=['',SCIENTIFIC_NOTE,'','Całki pre_baseline pochodzą ze sfazowanego widma przed korekcją bazy. '
             'Podpisanych całek nie przycinano. Maski fitów i processing zapisano w projekcie i tabelach. '
             'SE/CI są warunkowe względem modelu; nie obejmują kalibracji b.',
             '', 'Ważne: '+ '; '.join(dataset.metadata.get('limitations',[]))]
    (directory/'Raport.md').write_text('\n'.join(report),encoding='utf-8')
    return directory
