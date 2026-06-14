"""
Test 5: Jacobian operator linearity test with respect to model.

The Jacobian (Born/linearised forward) operator J(m) is linear in its
perturbation argument:
    (a) Scaling:        J(α·δm) = α·J(δm)
    (b) Superposition:  J(δm₁ + δm₂) = J(δm₁) + J(δm₂)

Both are verified numerically using finite-difference Born modelling:
    J(δm) · δm ≈ [F(m + ε·δm) − F(m)] / ε

The linearisation error is O(ε), so we use a small but not tiny ε to stay
well above round-off (ε = 1e-3).  The residuals (J(δm₁+δm₂)−J(δm₁)−J(δm₂))
and (J(α·δm₁)−α·J(δm₁)) should vanish as O(ε) → much smaller than the
Born data itself.

Pass criterion: relative residual < 5·ε = 5e-3 for both sub-tests.
"""

import os, sys, shutil
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from owl_test_utils import (
    write_model, write_geometry, write_param,
    run_owl, read_su, rel_err, report, Source
)

WORK = os.path.join(HERE, 'work', os.path.splitext(os.path.basename(__file__))[0])
PLOT = os.path.join(HERE, 'plots')
os.makedirs(WORK, exist_ok=True)
os.makedirs(PLOT, exist_ok=True)

NZ, NX = 101, 201
DZ, DX = 10.0, 10.0
VP0  = 2500.0
RHO  = 1.0
F0   = 15.0
DT   = 1.0e-3
TMAX = 1.5
EPS  = 0.02     # large enough for Born data >> float32 noise (~500× floor)

ALPHA    = 2.0      # scaling factor for sub-test 5a
PASS_TOL = 0.05     # 5% – accounts for O(ε) nonlinear truncation error

SX, SZ = NX // 2 * DX, NZ // 2 * DZ
RECVS  = [(DX * (NX // 4 + 10 * k), SZ) for k in range(8)]


def gaussian_perturbation(nz, nx, iz0, ix0, sigma, amp=1.0):
    iz = np.arange(nz)[:, None]
    ix = np.arange(nx)[None, :]
    return amp * np.exp(-((iz - iz0)**2 + (ix - ix0)**2) / sigma**2)


def _forward(tag, vp_arr):
    write_model(os.path.join(WORK, f'vp_{tag}.bin'), vp_arr)
    write_param(os.path.join(WORK, f'param_{tag}.rb'), {
        'nx': NX, 'nz': NZ, 'dx': DX, 'dz': DZ,
        'dt': DT, 'data_dt': DT, 'tmax': TMAX,
        'ns': 1,
        'file_geometry': './geometry/geometry.txt',
        'which_medium': 'acoustic-iso',
        'model_name': 'vp, rho',
        'file_vp':  f'./vp_{tag}.bin',
        'file_rho': './model/rho.bin',
        'dir_synthetic': f'./syn_{tag}',
        'verbose': 'n',
    })
    shutil.rmtree(os.path.join(WORK, f'syn_{tag}'), ignore_errors=True)
    run_owl('owl_modeling2', f'param_{tag}.rb', WORK)
    d, _ = read_su(os.path.join(WORK, f'syn_{tag}', 'shot_1_seismogram_p.su'))
    return d.astype(np.float64)


def run():
    vp0 = np.full((NZ, NX), VP0)
    rho = np.full((NZ, NX), RHO)
    write_model(os.path.join(WORK, 'model', 'rho.bin'), rho)
    write_geometry(
        os.path.join(WORK, 'geometry'),
        [Source(x=SX, z=SZ, f0=F0, amp=1.0)],
        RECVS
    )

    # Two independent Gaussian perturbations
    dvm1 = gaussian_perturbation(NZ, NX, NZ // 3,     NX // 3,     sigma=10.0, amp=0.2 * VP0)
    dvm2 = gaussian_perturbation(NZ, NX, 2 * NZ // 3, 2 * NX // 3, sigma=12.0, amp=0.2 * VP0)

    # Background data
    d0 = _forward('base', vp0)

    # Born data for δm₁, δm₂, δm₁+δm₂, α·δm₁
    d1   = (_forward('dm1',    vp0 + EPS * dvm1)        - d0) / EPS
    d2   = (_forward('dm2',    vp0 + EPS * dvm2)        - d0) / EPS
    d12  = (_forward('dm12',   vp0 + EPS * (dvm1+dvm2)) - d0) / EPS
    d1a  = (_forward('dm1a',   vp0 + EPS * ALPHA * dvm1) - d0) / EPS   # = J(α·δm₁)

    # Sub-test 5a: scaling  J(α·δm₁) = α·J(δm₁)
    err_5a = rel_err(d1a, ALPHA * d1)
    pass_5a = err_5a < PASS_TOL

    # Sub-test 5b: superposition  J(δm₁+δm₂) = J(δm₁)+J(δm₂)
    err_5b = rel_err(d12, d1 + d2)
    pass_5b = err_5b < PASS_TOL

    passed = pass_5a and pass_5b
    report('Test 5a – Jacobian scaling',       pass_5a, f'err = {err_5a:.2e}  (tol {PASS_TOL:.1e})')
    report('Test 5b – Jacobian superposition', pass_5b, f'err = {err_5b:.2e}  (tol {PASS_TOL:.1e})')

    # ── Visualization ──────────────────────────────────────────────────────────
    t_ax = np.arange(d0.shape[0]) * DT
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    def _show(ax, d, title, cmap='RdBu_r'):
        im = ax.imshow(d.T, aspect='auto', cmap=cmap,
                       vmin=-np.abs(d).max(), vmax=np.abs(d).max(),
                       extent=[0, d.shape[0] * DT, d.shape[1], 0])
        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.set_xlabel('Time (s)'); ax.set_ylabel('Receiver')
        ax.set_title(title, fontsize=9)

    _show(axes[0, 0], d1,  'J·δm₁  (Born data 1)')
    _show(axes[0, 1], d2,  'J·δm₂  (Born data 2)')
    _show(axes[0, 2], d12, 'J·(δm₁+δm₂)  (simultaneous)')

    _show(axes[1, 0], d1a - ALPHA * d1,
          f'J(α·δm₁) − α·J(δm₁)\nerr={err_5a:.1e}  '
          f'[{"PASS" if pass_5a else "FAIL"}]')
    _show(axes[1, 1], d12 - (d1 + d2),
          f'J(δm₁+δm₂) − J(δm₁) − J(δm₂)\nerr={err_5b:.1e}  '
          f'[{"PASS" if pass_5b else "FAIL"}]')

    # Model perturbations
    ax = axes[1, 2]
    im = ax.imshow(dvm1 + dvm2, aspect='auto', cmap='RdBu_r',
                   extent=[0, NX * DX / 1e3, NZ * DZ / 1e3, 0])
    plt.colorbar(im, ax=ax, label='δvp (m/s)', shrink=0.8)
    ax.set_title('δm₁ + δm₂\n(two Gaussian anomalies)')

    fig.suptitle('Test 5: Jacobian Linearity Test  '
                 f'[{"PASSED" if passed else "FAILED"}]\n'
                 f'ε={EPS:.0e}, α={ALPHA}', fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(PLOT, 'test5_jacobian_linearity.png'), dpi=150)
    plt.close()

    return passed, err_5a, err_5b


if __name__ == '__main__':
    passed, *_ = run()
    sys.exit(0 if passed else 1)
