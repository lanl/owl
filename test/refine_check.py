"""
Diagnostic: is the 8% across-δm spread in the gradient test a discretisation
effect (correct but "consistent" adjoint) or a genuine inconsistency (bug)?

Method — pure gradient (Taylor) test, NO Born operator required:
  • Build a synthetic "observed" record  d_obs = d0 − δd  (δd arbitrary), so the
    misfit  φ(m') = ‖F(m') − d_obs‖²  has residual δd at the background model m.
  • g = one-iteration FWI gradient of φ  (OWL adjoint-state gradient).
  • For each δm_k, measure the directional derivative of φ purely from the
    NONLINEAR forward by a central difference:
        D_k = [φ(m+ε·δm_k) − φ(m−ε·δm_k)] / (2ε)
    and form  C_k = D_k / ⟨g, δm_k⟩.
  • If g ∝ ∇φ with a single constant, C_k is the same for every δm_k.

We measure the spread  max_k|C_k − C_1| / |C_1|  on two grids (Δx, Δt) and
(Δx/2, Δt/2).  A correct-but-consistent adjoint → spread shrinks under
refinement.  A real inconsistency → spread stays put.

The whole physical setup is held fixed in metres/seconds; only the sampling
changes, so the two runs are directly comparable.
"""

import os, sys, shutil
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from owl_test_utils import (
    write_model, read_model, write_geometry, write_param, write_su,
    run_owl, read_su, Source
)

# ── Fixed physical setup (metres, seconds) ──────────────────────────────────────
LX, LZ   = 2000.0, 1000.0     # domain size
VP0, RHO = 2500.0, 1.0
F0       = 15.0
TMAX     = 1.5
AMP      = 1e4
SX, SZ   = 1000.0, 500.0      # source at centre
RECV_X   = [500.0 + 100.0 * k for k in range(8)]   # receivers at z = SZ
EPS_BORN = 0.02

# δd reference perturbation and the two test δm (physical centres/widths)
DD_REF  = dict(zx=(SZ, 1000.0), sigma=120.0, amp=0.1 * VP0)
TEST_DM = [
    dict(zx=(SZ, 600.0), sigma=60.0, amp=0.05 * VP0),   # far from source (x=1000)
    dict(zx=(SZ, 900.0), sigma=60.0, amp=0.06 * VP0),   # near source
]


def gaussian_phys(nz, nx, dz, dx, z0, x0, sigma, amp):
    iz = (np.arange(nz)[:, None] * dz - z0)
    ix = (np.arange(nx)[None, :] * dx - x0)
    return amp * np.exp(-(iz**2 + ix**2) / sigma**2)


def make_grid(factor):
    dx = dz = 10.0 / factor
    dt = 1.0e-3 / factor
    nx = int(round(LX / dx)) + 1
    nz = int(round(LZ / dz)) + 1
    return dict(dx=dx, dz=dz, dt=dt, nx=nx, nz=nz)


def run_factor(factor):
    G = make_grid(factor)
    dx, dz, dt, nx, nz = G['dx'], G['dz'], G['dt'], G['nx'], G['nz']
    work = os.path.join(HERE, f'work_refine_f{factor}')
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work, exist_ok=True)

    rho = np.full((nz, nx), RHO)
    write_model(os.path.join(work, 'model', 'rho.bin'), rho)
    write_geometry(os.path.join(work, 'geometry'),
                   [Source(x=SX, z=SZ, f0=F0, amp=AMP)],
                   [(x, SZ) for x in RECV_X])

    def forward(tag, vp):
        write_model(os.path.join(work, f'vp_{tag}.bin'), vp)
        write_param(os.path.join(work, f'param_{tag}.rb'), {
            'nx': nx, 'nz': nz, 'dx': dx, 'dz': dz,
            'dt': dt, 'data_dt': dt, 'tmax': TMAX, 'ns': 1,
            'file_geometry': './geometry/geometry.txt',
            'which_medium': 'acoustic-iso', 'model_name': 'vp, rho',
            'file_vp': f'./vp_{tag}.bin', 'file_rho': './model/rho.bin',
            'dir_synthetic': f'./syn_{tag}', 'verbose': 'n',
        })
        shutil.rmtree(os.path.join(work, f'syn_{tag}'), ignore_errors=True)
        run_owl('owl_modeling2', f'param_{tag}.rb', work)
        d, _ = read_su(os.path.join(work, f'syn_{tag}', 'shot_1_seismogram_p.su'))
        return d.astype(np.float64)

    vp0 = np.full((nz, nx), VP0)
    d0  = forward('base', vp0)

    # reference δd via central difference (high SNR, O(ε²))
    g_ref = gaussian_phys(nz, nx, dz, dx, DD_REF['zx'][0], DD_REF['zx'][1],
                          DD_REF['sigma'], DD_REF['amp'])
    dd = (forward('dref_p', vp0 + EPS_BORN * g_ref)
          - forward('dref_m', vp0 - EPS_BORN * g_ref)) / (2.0 * EPS_BORN)

    d_obs = d0 - dd

    # gradient g = ∇φ via 1 FWI iteration
    rdir = os.path.join(work, 'obs'); os.makedirs(rdir, exist_ok=True)
    syn_path = os.path.join(work, 'syn_base', 'shot_1_seismogram_p.su')
    write_su(os.path.join(rdir, 'shot_1_seismogram_p.su'), d_obs, dt, like=syn_path)
    write_model(os.path.join(work, 'fwi_vp.bin'), vp0)
    write_param(os.path.join(work, 'param_fwi.rb'), {
        'nx': nx, 'nz': nz, 'dx': dx, 'dz': dz,
        'dt': dt, 'data_dt': dt, 'tmax': TMAX, 'ns': 1,
        'file_geometry': './geometry/geometry.txt',
        'which_medium': 'acoustic-iso', 'model_update': 'vp',
        'file_vp': './fwi_vp.bin', 'niter_max': 1, 'misfit_type': 'waveform',
        'dir_record': './obs', 'dir_working': './fwi', 'verbose': 'n',
    })
    shutil.rmtree(os.path.join(work, 'fwi'), ignore_errors=True)
    run_owl('owl_fwi2', 'param_fwi.rb', work)
    g = read_model(os.path.join(work, 'fwi', 'iteration_1', 'model', 'grad_vp.bin'),
                   (nz, nx))

    # gradient test: C_k = directional-derivative(φ) / <g, δm_k>
    Cs = []
    for k, spec in enumerate(TEST_DM):
        dvm = gaussian_phys(nz, nx, dz, dx, spec['zx'][0], spec['zx'][1],
                            spec['sigma'], spec['amp'])
        dp = forward(f'dm{k}_p', vp0 + EPS_BORN * dvm)
        dm = forward(f'dm{k}_m', vp0 - EPS_BORN * dvm)
        # φ(m') = ||F(m') - d_obs||^2 ; central diff directional derivative
        phi_p = float(np.sum((dp - d_obs)**2))
        phi_m = float(np.sum((dm - d_obs)**2))
        D_k   = (phi_p - phi_m) / (2.0 * EPS_BORN)
        gd    = float(np.dot(dvm.ravel(), g.ravel()))
        Cs.append(D_k / gd)

    Cs = np.array(Cs)
    spread = np.max(np.abs(Cs - Cs[0])) / np.abs(Cs[0])
    print(f'  factor {factor}: dx={dx:.2f} m  dt={dt*1e3:.3f} ms  nz×nx={nz}×{nx}')
    for k, c in enumerate(Cs):
        print(f'      C_{k+1} = {c:.6e}')
    print(f'      across-δm spread = {spread*100:.3f} %')
    shutil.rmtree(work, ignore_errors=True)
    return spread


if __name__ == '__main__':
    print('Grid-refinement gradient test (no Born operator):')
    s1 = run_factor(1)
    s2 = run_factor(2)
    print()
    print(f'  spread(Δx)    = {s1*100:.3f} %')
    print(f'  spread(Δx/2)  = {s2*100:.3f} %')
    print(f'  ratio s1/s2   = {s1/s2:.2f}   '
          f'({"shrinks → discretisation (gradient OK)" if s2 < 0.7*s1 else "does NOT shrink → investigate"})')
