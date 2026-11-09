"""Read-only JEOL classical DOSY reader. No inference from an expected D."""
from dataclasses import dataclass, replace
from pathlib import Path
import hashlib,re,struct
import numpy as np

GAMMA_H=267522187.44


@dataclass
class RawData:
    fid: np.ndarray                 # acquired layer order, complex; never modified
    gradients: np.ndarray           # T/m, NaN when unavailable
    b: np.ndarray | None
    metadata: dict

    def __post_init__(self):
        self.fid.setflags(write=False)
        self.gradients.setflags(write=False)
        if self.b is not None:self.b.setflags(write=False)


def reference_layer(raw):
    """Lowest b, otherwise lowest known encoding gradient; keep acquired order."""
    if raw.b is not None:return int(np.argmin(raw.b))
    valid=np.isfinite(raw.gradients)
    return int(np.flatnonzero(valid)[np.argmin(abs(raw.gradients[valid]))]) if valid.any() else 0


def validate_b(values,count):
    b=np.asarray(values,dtype=float)
    if b.shape!=(count,) or not np.isfinite(b).all() or (b<0).any():
        raise ValueError(f'Podaj dokładnie {count} skończonych, nieujemnych b w s/m².')
    return b.copy()


def manual_b(raw,values,reason):
    if not reason.strip():raise ValueError('Podaj źródło lub uzasadnienie ręcznej listy b.')
    meta=dict(raw.metadata)
    history=list(meta.get('b_history',[]))
    history.append(meta.get('b_source'))
    meta.update(b_source={'kind':'manual','description':reason.strip(),'calibrated':False},b_history=history)
    return replace(raw,b=validate_b(values,len(raw.fid)),metadata=meta)


def number(params,name,base=None):
    p=params[name]
    active=[u for u in p['units'] if u['base_code'] or u['power']]
    if base is not None and (len(active)!=1 or active[0]['base_code']!=base or active[0]['power']!=1):
        raise ValueError('Niewłaściwa jednostka parametru '+name)
    value=float(p['value'])*10.**p['decimal_scaler']
    # Unit prefixes are powers of 1000 (1=milli, 2=micro).
    if base!=26:value*=10.**(-3*sum(u['prefix_code']*u['power'] for u in active))
    if not np.isfinite(value):raise ValueError('Nieskończony parametr '+name)
    return value


def parameters(blob,start,length,endian):
    if start<1360 or length<16 or start+length>len(blob):raise ValueError('Uszkodzona sekcja parametrów JDF.')
    size,lo,hi,total=struct.unpack_from(endian+'4I',blob,start)
    count=(length-16)//64
    if size!=64 or hi-lo+1!=count or total!=64*count or (length-16)%64:
        raise ValueError('Nieobsługiwany układ parametrów JDF.')
    result={}
    for at in range(start+16,start+length,64):
        kind=struct.unpack_from(endian+'i',blob,at+32)[0]
        units=[]
        for code,base in zip(blob[at+6:at+16:2],blob[at+7:at+16:2]):
            pre,power=code>>4,code&15
            units.append(dict(prefix_code=pre-16 if pre>=8 else pre,power=power-16 if power>=8 else power,base_code=base))
        raw=blob[at+16:at+32]
        if kind==0:value=raw.split(b'\0')[0].decode('utf-8',errors='replace').rstrip()
        elif kind in (1,4):value=struct.unpack_from(endian+'i',raw)[0]
        elif kind==2:value=struct.unpack_from(endian+'d',raw)[0]
        elif kind==3:value=list(struct.unpack_from(endian+'2d',raw))
        else:value=None
        name=blob[at+36:at+64].decode('ascii',errors='replace').strip().lower()
        result[name]=dict(value=value,units=units,decimal_scaler=struct.unpack_from(endian+'h',blob,at+4)[0],value_type=kind)
    return result


def _axis_list(blob,axis,count,start,stop):
    offset=struct.unpack_from('>I',blob,1220+4*axis)[0]
    length=struct.unpack_from('>I',blob,1252+4*axis)[0]
    if not length:return None,None
    valid=stop-start+1
    if length not in (8*count,8*valid) or offset<1360 or offset+length>len(blob):
        raise ValueError('Niezgodna lista osi warstw JDF.')
    values=np.frombuffer(blob,dtype='>f8',count=length//8,offset=offset).copy()
    if len(values)==count:values=values[start:stop+1]
    slot=blob[1328+axis]
    if slot not in (1,2):return values,None
    scaler,*codes=struct.unpack_from('>h5H',blob,1336+12*(slot-1))
    active=[c for c in codes if c]
    if len(active)==2 and active[0]&0xfff==0x11f and active[1]==0xf13:
        return values*10.**(scaler-3*(active[0]>>12)),'gradient_T_m'
    if len(active)==2 and active[0]&0xfff==0x11c and active[1]==0xe13:
        return values*10.**(scaler-3*(active[0]>>12)),'b_s_m2'
    return values,None


def _gradient_source(context,n):
    match=re.search(r'\bg\s*=>\s*y_acq\s*\{([^}]+)\}',context,re.I)
    if not match:return None
    values=[]
    for item in match[1].split(','):
        m=re.fullmatch(r'\s*([+\-\d.eE]+)\s*\[\s*(T/m|mT/m|uT/m)\s*\]\s*',item)
        if not m:return None
        values.append(float(m[1])*{'T/m':1,'mT/m':1e-3,'uT/m':1e-6}[m[2]])
    if len(values)!=n:raise ValueError('Liczba gradientów w programie nie zgadza się z warstwami.')
    return np.array(values)


def factory_b(params,context,g):
    identity=re.search(r'\bexperiment\s*=\s*"([^"]+)"',context,re.I)
    name=Path(identity[1].replace('\\','/')).name.lower() if identity else ''
    families={'bpp_ste_led_diffusion.jxp':'BPP-STE-LED','bpp_dste_led_diffusion.jxp':'BPP-DSTE-LED'}
    family=families.get(name)
    info={'sequence':family or name or 'Nierozpoznana','kind':'unavailable',
          'description':'Brak zweryfikowanego modelu b. Wprowadź własną listę b.', 'calibrated':False}
    if family is None:return None,info
    try:
        if params.get('x_domain',{}).get('value')!='Proton':raise ValueError('Model zweryfikowano dla kodowania ¹H.')
        if not np.isfinite(g).all() or (g<0).any():raise ValueError('Brak wiarygodnej listy gradientów kodujących.')
        if not re.search(r'grad_1\s*=\?\s*delta\s*/\s*2',context,re.I):raise ValueError('Niepotwierdzona konwencja δ/2.')
        shape=params['grad_shape']['value']
        constants={'SQUARE':(1,1/3),'SINE':(4/np.pi**2,5/16),
                   'SQUARE_SINE':(.25,1/3-5/(16*np.pi**2)),
                   'TRAPEZOID':(.81,15817/48600),
                   'S_RECTANGLE':(.81,(6344*np.pi**2-207)/(19440*np.pi**2))}
        if shape not in constants:raise ValueError('Nierozpoznany kształt gradientu: '+str(shape))
        delta=number(params,'delta',28);Delta=number(params,'delta_large',28)
        if family=='BPP-DSTE-LED':
            tau=number(params,'tau_1',28)+number(params,'tau_2',28);pair=2
        else:tau=number(params,'tau',28);pair=1
        factor,correction=constants[shape]
        expected=factor*(Delta-pair*correction*delta-.5*tau)
        if delta<=0 or Delta<=pair*delta or expected<=0:raise ValueError('Nieprawidłowe czasy kodowania.')
        if 'st_coef' in params:
            coefficient=number(params,'st_coef',28)
            if coefficient<=0 or not np.isclose(coefficient,expected,rtol=3e-5,atol=1e-10):
                raise ValueError('Zapisany st_coef nie zgadza się z rozpoznanym modelem fabrycznym.')
            source='JDF st_coef, sprawdzony względem modelu '+family
        else:coefficient=expected;source='Model fabryczny '+family+' z zapisanych timingów'
        b=(GAMMA_H*g*delta)**2*coefficient
        info.update(kind='factory',description=source,gradient_shape=shape,
                    delta_s=delta,Delta_s=Delta,coefficient_s=coefficient,
                    expected_coefficient_s=expected,gamma_rad_s_T=GAMMA_H,
                    note='Nominalne b; kalibracja amplitud/odpowiedzi gradientów nie jest weryfikowana przez program.')
        return validate_b(b,len(g)),info
    except (KeyError,TypeError,ValueError) as exc:
        info['description']=str(exc)+' Wprowadź własną listę b.'
        return None,info


def infer_filter_delay(params, first_point, npoints):
    """Infer a linear-phase decimation delay for the supported raw JEOL metadata."""
    value=lambda key:params.get(key,{}).get('value')
    if (first_point!=0 or value('digital_filter')!='TRUE'
            or value('clipped')!='FALSE' or value('digital_filter_status')!='Np'):
        return None
    try:
        stages,*taps=map(int,value('orders').split())
        factors=list(map(int,value('factors').split()))
    except (AttributeError,ValueError,TypeError):
        return None
    if not 1<=stages<=8 or len(taps)!=stages or len(factors)!=stages:
        return None
    if any(not 1<=n<=65536 for n in taps+factors):
        return None
    decimation=1;delay=0.
    for length,factor in zip(reversed(taps),reversed(factors)):
        decimation*=factor
        delay+=(length-1)/(2*decimation)
    if decimation!=value('filter_factor') or decimation<2 or not 0<delay<=250 or 2*delay+8>=npoints:
        return None
    return dict(samples=delay,taps=taps,factors=factors,
                source='Model kaskady liniowofazowej z JDF orders/factors (oszacowanie)')


def read_jdf(path):
    path=Path(path).resolve();blob=path.read_bytes();digest=hashlib.sha256(blob).hexdigest()
    if len(blob)<1360 or blob[:8]!=b'JEOL.NMR' or blob[8] not in (0,1):raise ValueError('To nie jest poprawny JEOL JDF.')
    dimension,code,dtype_code=blob[12],blob[14]&63,blob[14]>>6
    if dimension not in (1,2):raise ValueError('DOSY Lab otwiera klasyczne FID 1D × warstwy; ten plik nie jest takim pomiarem.')
    if dtype_code not in (0,1) or blob[24] not in (3,4) or (dimension==2 and blob[25]!=1):
        raise ValueError('Nieobsługiwane typy osi lub brak zespolonego FID.')
    endian='<' if blob[8] else '>'
    dtype=np.dtype(endian+('f8' if dtype_code==0 else 'f4'))
    pts=struct.unpack_from('>8I',blob,176)[:dimension]
    starts=struct.unpack_from('>8I',blob,208)[:dimension];stops=struct.unpack_from('>8I',blob,240)[:dimension]
    ds=struct.unpack_from('>I',blob,1284)[0];dl=struct.unpack_from('>Q',blob,1288)[0]
    cs=struct.unpack_from('>Q',blob,1296)[0];cl=struct.unpack_from('>I',blob,1304)[0]
    declared=struct.unpack_from('>Q',blob,1320)[0]
    if declared!=len(blob) or ds<1360 or ds+dl>len(blob) or cs+cl>len(blob):raise ValueError('Uszkodzone rozmiary sekcji JDF.')
    if any(a>b or b>=n for n,a,b in zip(pts,starts,stops)):raise ValueError('Nieprawidłowy zakres ważnych punktów JDF.')
    if dl!=2*int(np.prod(pts))*dtype.itemsize:raise ValueError('Niezgodna liczba sekcji zespolonych JDF.')
    raw=np.frombuffer(blob,dtype=dtype,count=dl//dtype.itemsize,offset=ds)
    if dimension==1 and code==1:
        sections=raw.reshape(2,pts[0])[:,None,:]
    elif dimension==2 and code in (2,12):
        edge=4 if code==12 else 32
        nx,ng=pts
        if nx%edge or ng%edge:raise ValueError('Wymiary nie pasują do kafli JDF.')
        sections=raw.reshape(2,ng//edge,nx//edge,edge,edge).transpose(0,1,3,2,4).reshape(2,ng,nx)
        repacked=sections.reshape(2,ng//edge,edge,nx//edge,edge).transpose(0,1,3,2,4).reshape(-1)
        if repacked.tobytes()!=raw.tobytes():raise ValueError('Kontrola rekonstrukcji kafli nie powiodła się.')
    else:raise ValueError(f'Nieobsługiwany format JDF {code}; nie zgadujemy kolejności danych.')
    pstart,plen=struct.unpack_from('>2I',blob,1212);params=parameters(blob,pstart,plen,endian)
    context='\n'.join(s.decode('ascii') for s in re.findall(rb'[\x09\x0a\x0d\x20-\x7e]{6,}',blob[cs:cs+cl]))
    xsl=slice(starts[0],stops[0]+1);gsl=slice(starts[1],stops[1]+1) if dimension==2 else slice(None)
    fid=np.array((sections[0]-1j*sections[1])[gsl,xsl],copy=True)
    if not np.isfinite(fid).all():raise ValueError('FID zawiera NaN/Inf.')
    ng=len(fid);g=np.full(ng,np.nan);b=None;axis_kind=None
    if dimension==2:
        axis,axis_kind=_axis_list(blob,1,pts[1],starts[1],stops[1])
        source=_gradient_source(context,ng)
        if axis_kind=='gradient_T_m':
            g=axis
            if source is not None and not np.allclose(g,source,atol=5.1e-6,rtol=0):raise ValueError('Lista osi gradientów nie zgadza się z programem.')
        elif source is not None:g=source
        if axis_kind=='b_s_m2':b=validate_b(axis,ng)
    inferred,model=factory_b(params,context,g)
    if b is None:b=inferred
    else:model={'kind':'stored_axis','description':'Jawna lista b z osi JDF [s/m²]','sequence':model['sequence'],'calibrated':False}
    sw=number(params,'x_sweep',13);freq=number(params,'x_freq',13);off=number(params,'x_offset',26)
    if min(sw,freq)<=0:raise ValueError('Nieprawidłowa częstotliwość lub szerokość widma.')
    # Report encoding and spoilers separately; never use spoiler amplitudes as g.
    body=context.split('\nbegin\n',1)[-1].split('end pulse;',1)[0]
    events=[]
    for line in body.splitlines():
        if 'FGZ.GATE' in line.upper():
            role='encoding' if re.search(r'FGZ\.AMP\.\s*-?\s*grad_1_amp',line,re.I) else 'other_or_spoiler'
            events.append({'instruction':line.strip(),'role':role})
    meta={'source_name':path.name,'source_path':str(path),'source_sha256':digest,'parameters':params,
          'sequence':model['sequence'],'b_source':model,'gradient_axis_source':axis_kind or 'source_or_missing',
          'gradient_events':events,'x_sweep_Hz':sw,'x_frequency_Hz':freq,'offset_ppm':off,
          'x_sweep_clipped_Hz':number(params,'x_sweep_clipped',13) if 'x_sweep_clipped' in params else sw,
          'digital_filter':params.get('digital_filter',{}).get('value')=='TRUE',
          'already_clipped':params.get('clipped',{}).get('value')=='TRUE',
          'filter_delay_model':infer_filter_delay(params,starts[0],fid.shape[1]),
          'complex_convention':'section_0 - 1j*section_1','stored_shape':list(pts[::-1]),
          'valid_shape':list(fid.shape),'context':context,
          'limitations':['Oś ppm z x_offset i x_freq; odnośnik chemiczny wymaga potwierdzenia.',
                         'Korekcja fazy filtra to model opóźnienia, bez odtwarzania pełnej odpowiedzi Delty.']}
    return RawData(fid,g,b,meta)
