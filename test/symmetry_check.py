"""
Decisive gradient check via SYMMETRY (no Born operator, no source-singularity
confound).

The geometry — single source plus all receivers on the line z = 500 m — is
exactly mirror-symmetric under  z → (LZ − z).  Therefore a correct gradient
MUST satisfy  g(z, x) = g(LZ − z, x)  to machine precision, and any pair of
mirror-image model perturbations must yield identical directional derivatives.

  (A) Direct test:   ‖g − flipud(g)‖ / ‖g‖   should be ~1e-6 (float32 round-off).
  (B) Cross test:    for δm at (z0,x) and its mirror (LZ−z0, x),
                     C = D/⟨g,δm⟩ must match.

A real asymmetry in the adjoint implementation shows up here regardless of
grid or source effects; discretisation cannot break a symmetry the scheme
respects.
"""

import os, sys, shutil
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from owl_test_utils import (
    write_model, read_model, write_geometry, write_param, write_su,
    run_owl, read_su, Source
)

LX, LZ   = 2000.0, 1000.0
VP0, RHO = 2500.0, 1.0
F0, TMAX, AMP = 15.0, 1.5, 1e4
DX = DZ = 10.0
DT = 1.0e-3
NX = int(round(LX / DX)) + 1
NZ = int(round(LZ / DZ)) + 1
SX, SZ = 1000.0, 500.0
RECV_X = [500.0 + 100.0 * k for k in range(8)]
EPS = 0.02

WORK = os.path.join(HERE, 'work_sym')


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

    # (A) direct symmetry of the gradient about z = LZ/2  (row index NZ-1-i)
    g_mirror = g[::-1, :]
    asym = np.linalg.norm(g - g_mirror) / np.linalg.norm(g)
    print(f'(A) gradient mirror asymmetry  ‖g−flipud(g)‖/‖g‖ = {asym:.3e}')
    print(f'    (float32 round-off floor ~1e-6; >>1e-3 indicates a real asymmetry)')

    # (B) cross test with mirror-pair perturbations off the receiver line
    def C_of(z0, x0):
        dvm = gaussian_phys(z0, x0, 60.0, 0.05 * VP0)
        dp = forward('dp', vp0 + EPS * dvm)
        dm = forward('dm', vp0 - EPS * dvm)
        D  = (np.sum((dp - d_obs)**2) - np.sum((dm - d_obs)**2)) / (2.0 * EPS)
        return float(D / np.dot(dvm.ravel(), g.ravel()))

    print('(B) mirror-pair directional-derivative ratios:')
    for x0 in (700.0, 1000.0, 1300.0):
        c_up = C_of(400.0, x0)   # 100 m above the line
        c_dn = C_of(600.0, x0)   # 100 m below the line
        rel  = abs(c_up - c_dn) / abs(c_up)
        print(f'    x={x0:6.0f} m:  C(z=400)={c_up:.5e}  C(z=600)={c_dn:.5e}  '
              f'rel-diff={rel*100:.3f} %')

    shutil.rmtree(WORK, ignore_errors=True)


if __name__ == '__main__':
    main()
