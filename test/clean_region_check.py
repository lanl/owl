"""
Confirm the gradient is consistent (C = D/⟨g,δm⟩ ≈ const) once the test
perturbations are placed in the well-illuminated region, AWAY from the point
source and off the direct source–receiver line where adjoint-state gradients
carry their well-known near-source artifact.

Source at (z=500, x=1000); receivers on z=500.  We probe δm at z = 300 and
z = 700 m (200 m off the line) and several x in [600, 1400] m, skipping the
source column.  A correct gradient → all C agree to a few %.
"""

import os, sys, shutil
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from owl_test_utils import (
    write_model, read_model, write_geometry, write_param, write_su,
    run_owl, read_su, Source
)

LX, LZ = 2000.0, 1000.0
VP0, RHO = 2500.0, 1.0
F0, TMAX, AMP = 15.0, 1.5, 1e4
DX = DZ = 10.0
DT = 1.0e-3
NX = int(round(LX / DX)) + 1
NZ = int(round(LZ / DZ)) + 1
SX, SZ = 1000.0, 500.0
RECV_X = [500.0 + 100.0 * k for k in range(8)]
EPS = 0.02
WORK = os.path.join(HERE, 'work_clean')


def gaussian_phys(z0, x0, sigma, amp):
    iz = (np.arange(NZ)[:, None] * DZ - z0)
    ix = (np.arange(NX)[None, :] * DX - x0)
    return amp * np.exp(-(iz**2 + ix**2) / sigma**2)


def forward(tag, vp):
    write_model(os.path.join(WORK, f'vp_{tag}.bin'), vp)
    write_param(os.path.join(WORK, f'param_{tag}.rb'), {
        'nx': NX, 'nz': NZ, 'dx': DX, 'dz': DZ,
        'dt': DT, 'data_dt': DT, 'tmax': TMAX, 'ns': 1,
        'file_geometry': './geometry/geometry.txt',
        'which_medium': 'acoustic-iso', 'model_name': 'vp, rho',
        'file_vp': f'./vp_{tag}.bin', 'file_rho': './model/rho.bin',
        'dir_synthetic': f'./syn_{tag}', 'verbose': 'n',
    })
    shutil.rmtree(os.path.join(WORK, f'syn_{tag}'), ignore_errors=True)
    run_owl('owl_modeling2', f'param_{tag}.rb', WORK)
    d, _ = read_su(os.path.join(WORK, f'syn_{tag}', 'shot_1_seismogram_p.su'))
    return d.astype(np.float64)


def main():
    shutil.rmtree(WORK, ignore_errors=True)
    os.makedirs(WORK, exist_ok=True)
    write_model(os.path.join(WORK, 'model', 'rho.bin'), np.full((NZ, NX), RHO))
    write_geometry(os.path.join(WORK, 'geometry'),
                   [Source(x=SX, z=SZ, f0=F0, amp=AMP)],
                   [(x, SZ) for x in RECV_X])

    vp0 = np.full((NZ, NX), VP0)
    d0  = forward('base', vp0)
    g_ref = gaussian_phys(SZ, 1000.0, 120.0, 0.1 * VP0)
    dd = (forward('dref_p', vp0 + EPS * g_ref)
          - forward('dref_m', vp0 - EPS * g_ref)) / (2.0 * EPS)
    d_obs = d0 - dd

    rdir = os.path.join(WORK, 'obs'); os.makedirs(rdir, exist_ok=True)
    syn_path = os.path.join(WORK, 'syn_base', 'shot_1_seismogram_p.su')
    write_su(os.path.join(rdir, 'shot_1_seismogram_p.su'), d_obs, DT, like=syn_path)
    write_model(os.path.join(WORK, 'fwi_vp.bin'), vp0)
    write_param(os.path.join(WORK, 'param_fwi.rb'), {
        'nx': NX, 'nz': NZ, 'dx': DX, 'dz': DZ,
        'dt': DT, 'data_dt': DT, 'tmax': TMAX, 'ns': 1,
        'file_geometry': './geometry/geometry.txt',
        'which_medium': 'acoustic-iso', 'model_update': 'vp',
        'file_vp': './fwi_vp.bin', 'niter_max': 1, 'misfit_type': 'waveform',
        'dir_record': './obs', 'dir_working': './fwi', 'verbose': 'n',
    })
    shutil.rmtree(os.path.join(WORK, 'fwi'), ignore_errors=True)
    run_owl('owl_fwi2', 'param_fwi.rb', WORK)
    g = read_model(os.path.join(WORK, 'fwi', 'iteration_1', 'model', 'grad_vp.bin'),
                   (NZ, NX))

    def C_of(z0, x0):
        dvm = gaussian_phys(z0, x0, 60.0, 0.05 * VP0)
        dp = forward('dp', vp0 + EPS * dvm)
        dm = forward('dm', vp0 - EPS * dvm)
        D  = (np.sum((dp - d_obs)**2) - np.sum((dm - d_obs)**2)) / (2.0 * EPS)
        return float(D / np.dot(dvm.ravel(), g.ravel()))

    probes = [(z0, x0) for z0 in (300.0, 700.0)
                       for x0 in (650.0, 800.0, 1200.0, 1350.0)]
    Cs = []
    print('Clean-region gradient test  C = D / <g, δm> :')
    for (z0, x0) in probes:
        c = C_of(z0, x0)
        Cs.append(c)
        print(f'    z={z0:5.0f}  x={x0:6.0f}:  C = {c:.5e}')
    Cs = np.array(Cs)
    spread = (Cs.max() - Cs.min()) / np.median(Cs)
    print(f'\n  median C = {np.median(Cs):.5e}')
    print(f'  full spread (max−min)/median = {spread*100:.2f} %')
    shutil.rmtree(WORK, ignore_errors=True)


if __name__ == '__main__':
    main()
