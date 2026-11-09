"""Generate a synthetic JDF for demonstrating DOSY Lab."""
import argparse
from pathlib import Path
import struct

import numpy as np

from jeol import GAMMA_H


def encode_fixture(fid, gradients, parameters, context, *, edge=4, little=True, dtype_code=0):
    """Pack the subset of JDF needed by the demonstration and reader tests."""
    layers, valid_points = fid.shape
    points = ((valid_points + edge - 1) // edge) * edge
    padded = np.pad(fid, ((0, 0), (0, points - valid_points)))
    packed = []
    for section in (padded.real, -padded.imag):
        for y in range(0, layers, edge):
            for x in range(0, points, edge):
                packed.extend(section[y:y+edge, x:x+edge].ravel())
    endian = '<' if little else '>'
    dtype = endian + ('f8' if dtype_code == 0 else 'f4')
    data = np.array(packed, dtype=dtype).tobytes()

    records = []
    for name, value in parameters.items():
        record = bytearray(64)
        if name in ('delta', 'delta_large', 'tau', 'tau_1', 'tau_2', 'st_coef'):
            record[6:8] = bytes((0x11, 28))
        elif name in ('x_sweep', 'x_freq'):
            record[6:8] = bytes((1, 13))
        elif name == 'x_offset':
            record[6:8] = bytes((1, 26))
        if isinstance(value, str):
            record[16:32] = value.encode('ascii').ljust(16, b'\0')
        else:
            struct.pack_into(endian + 'd', record, 16, value)
            struct.pack_into(endian + 'i', record, 32, 2)
        record[36:64] = name.encode('ascii').ljust(28, b' ')
        records.append(record)
    params = struct.pack(endian + '4I', 64, 0, len(records)-1, 64*len(records)) + b''.join(records)
    context = context.encode('ascii')
    header = bytearray(1360)
    header[:8] = b'JEOL.NMR'
    header[8] = int(little)
    header[12] = 2
    header[14] = (dtype_code << 6) | (12 if edge == 4 else 2)
    header[24:26] = bytes((3, 1))
    struct.pack_into('>8I', header, 176, points, layers, 1, 1, 1, 1, 1, 1)
    struct.pack_into('>8I', header, 240, valid_points-1, layers-1, 0, 0, 0, 0, 0, 0)
    struct.pack_into('>2I', header, 1212, 1360, len(params))
    axis_start = 1360 + len(params)
    axis = (gradients*1000).astype('>f8').tobytes()
    struct.pack_into('>I', header, 1224, axis_start)
    struct.pack_into('>I', header, 1256, len(axis))
    header[1329] = 1
    struct.pack_into('>h5H', header, 1336, 0, 0x111f, 0xf13, 0, 0, 0)
    data_start = axis_start + len(axis)
    struct.pack_into('>I', header, 1284, data_start)
    struct.pack_into('>Q', header, 1288, len(data))
    context_start = data_start + len(data)
    struct.pack_into('>Q', header, 1296, context_start)
    struct.pack_into('>I', header, 1304, len(context))
    struct.pack_into('>Q', header, 1320, context_start + len(context))
    return header + params + axis + data + context


def demo_bytes():
    gradients = np.linspace(.05, .4, 16)
    sw = 3000.
    delta = .002
    correction = (6344*np.pi**2-207)/(19440*np.pi**2)
    coefficient = .81 * (.1 - correction*delta - .001)
    b = (GAMMA_H*gradients*delta)**2 * coefficient
    t = np.arange(8191)/sw
    rng = np.random.default_rng(912)
    fid = np.zeros((len(gradients),len(t)),dtype=complex)
    for ppm, amplitude, diffusion in (
        (3.624,.45,3.8e-10), (3.649,.65,3.8e-10), (3.78,.42,3.8e-10),
        (4.08,.55,3.8e-10), (5.29,.75,3.8e-10), (4.70,1.05,8e-10),
    ):
        signal = amplitude * np.exp(2j*np.pi*(ppm-4.5)*600*t - np.pi*2*t)
        fid += np.exp(-b[:,None]*diffusion)*signal[None,:]
    fid += .0005*(rng.normal(size=fid.shape)+1j*rng.normal(size=fid.shape))
    params = dict(x_domain='Proton',x_sweep=sw,x_freq=600e6,x_offset=4.5,
                  delta=delta*1000,delta_large=100.,tau=2.,grad_shape='S_RECTANGLE',
                  st_coef=coefficient*1000,digital_filter='FALSE')
    context = ('experiment = "bpp_ste_led_diffusion.jxp";\ngrad_1 =? delta / 2;\n'
               'g => y_acq {'+', '.join(f'{g:.12g}[T/m]' for g in gradients)+'};\n'
               'begin\nend pulse;')
    return encode_fixture(fid,gradients,params,context)


def main():
    parser = argparse.ArgumentParser(description='Utwórz syntetyczny przykład DOSY w formacie JDF.')
    parser.add_argument('path',nargs='?',type=Path,default=Path('DOSY_demo.jdf'))
    path = parser.parse_args().path
    with path.open('xb') as stream:
        stream.write(demo_bytes())
    print(path.resolve())


if __name__ == '__main__':
    main()
