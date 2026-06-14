"""
Test 4: Nonlinear operator linearization test with respect to model.

Taylor's theorem guarantees that if F is Fréchet differentiable:
    F(m + ε·δm) = F(m) + ε·J(m)·δm + O(ε²)

Define the Taylor errors:
    err₀(ε) = ||F(m + ε·δm) − F(m)||        (first-order residual)
    err₁(ε) = ||F(m + ε·δm) − F(m) − ε·δd|| (second-order residual)

where δd = J(m)·δm is approximated by finite differences with a tiny ε_ref.

Expected convergence as ε → 0 (halving ε each step):
    err₀ ~ O(ε)     → slope ≈ 1 on a log–log plot
    err₁ ~ O(ε²)    → slope ≈ 2 on a log–log plot

Pass criterion: the slope of err₁ in log–log space is ≥ 1.9 (≈ 2).
"""

import os, sys, shutil
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from owl_test_utils import (
    write_model, read_model, write_geometry, write_param,
    run_owl, read_su, report, Source
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

SX, SZ = NX // 2 * DX, NZ // 2 * DZ
RECVS  = [(DX * (NX // 4 + 10 * k), SZ) for k in range(8)]

# Reference eps for computing the Born (Jacobian) data
EPS_REF  = 5e-4
# Test eps values: 0.5, 0.25, 0.125, 0.0625, 0.03125
EPS_LIST = [0.5 ** k for k in range(5)]


def gaussian_perturbation(nz, nx, iz0, ix0, sigma, amp=1.0):
    iz = np.arange(nz)[:, None]
    ix = np.arange(nx)[None, :]
    return amp * np.exp(-((iz - iz0)**2 + (ix - ix0)**2) / sigma**2)


def _forward(tag, vp_arr):
    geom_dir = 'geometry'
    syn_dir  = f'syn_{tag}'
    write_model(os.path.join(WORK, f'vp_{tag}.bin'), vp_arr)
    write_param(os.path.join(WORK, f'param_{tag}.rb'), {
        'nx': NX, 'nz': NZ, 'dx': DX, 'dz': DZ,
        'dt': DT, 'data_dt': DT, 'tmax': TMAX,
        'ns': 1,
        'file_geometry': f'./{geom_dir}/geometry.txt',
        'which_medium': 'acoustic-iso',
        'model_name': 'vp, rho',
        'file_vp':  f'./vp_{tag}.bin',
        'file_rho': './model/rho.bin',
        'dir_synthetic': f'./{syn_dir}',
        'verbose': 'n',
    })
    shutil.rmtree(os.path.join(WORK, syn_dir), ignore_errors=True)
    run_owl('owl_modeling2', f'param_{tag}.rb', WORK)
    d, _ = read_su(os.path.join(WORK, syn_dir, 'shot_1_seismogram_p.su'))
    return d.astype(np.float64)


def run():
    # Models
    vp0 = np.full((NZ, NX), VP0)
    rho = np.full((NZ, NX), RHO)
    write_model(os.path.join(WORK, 'model', 'rho.bin'), rho)

    # Model perturbation: smooth Gaussian anomaly in the centre
    dvm_norm = gaussian_perturbation(NZ, NX, NZ // 2, NX // 2, sigma=15.0)
    dvm_amp  = 0.1 * VP0   # 10 % of background velocity
    dvm      = dvm_norm * dvm_amp

    write_geometry(
        os.path.join(WORK, 'geometry'),
        [Source(x=SX, z=SZ, f0=F0, amp=1.0)],
        RECVS
    )

    # Background data
    d0 = _forward('base', vp0)

    # Reference Born data at very small eps
    d_ref = _forward('ref', vp0 + EPS_REF * dvm)
    d_born = (d_ref - d0) / EPS_REF   # Jacobian action approximation

    err0_list, err1_list = [], []
    for eps in EPS_LIST:
        d_eps = _forward(f'eps{eps:.5f}'.replace('.', 'p'), vp0 + eps * dvm)
        e0 = np.linalg.norm((d_eps - d0).ravel())
        e1 = np.linalg.norm((d_eps - d0 - eps * d_born).ravel())
        err0_list.append(e0)
        err1_list.append(e1)
        print(f'  ε={eps:.4f}  err₀={e0:.3e}  err₁={e1:.3e}  ratio={e1/e0:.3e}')

    eps_arr  = np.array(EPS_LIST)
    err0_arr = np.array(err0_list)
    err1_arr = np.array(err1_list)

    # Log-log slopes via linear regression (skip first point which can be large)
    slope0 = np.polyfit(np.log(eps_arr), np.log(err0_arr), 1)[0]
    slope1 = np.polyfit(np.log(eps_arr), np.log(err1_arr), 1)[0]

    passed = slope1 >= 1.9
    report('Test 4 – Linearization (Taylor test)', passed,
           f'slope(err₀)={slope0:.2f} (expect ≈1)  '
           f'slope(err₁)={slope1:.2f} (expect ≈2)')

    # ── Visualization ──────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    ax.loglog(eps_arr, err0_arr, 'bs-', label=f'err₀  (slope={slope0:.2f})')
    ax.loglog(eps_arr, err1_arr, 'r^-', label=f'err₁  (slope={slope1:.2f})')
    # Reference lines
    ref = err0_arr[0] * eps_arr / eps_arr[0]
    ax.loglog(eps_arr, ref,            'b--', alpha=0.4, label='O(ε)')
    ax.loglog(eps_arr, ref * eps_arr / eps_arr[0], 'r--', alpha=0.4, label='O(ε²)')
    ax.set_xlabel('ε'); ax.set_ylabel('‖residual‖₂')
    ax.set_title('Taylor test: convergence rates')
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)

    ax = axes[1]
    # Show δm perturbation
    im = ax.imshow(dvm, aspect='auto', cmap='RdBu_r',
                   extent=[0, NX * DX / 1e3, NZ * DZ / 1e3, 0])
    plt.colorbar(im, ax=ax, label='δvp (m/s)')
    ax.set_xlabel('X (km)'); ax.set_ylabel('Z (km)')
    ax.set_title('Model perturbation δm\n(Gaussian, 10% of vp)')

    fig.suptitle('Test 4: Nonlinear Operator Linearization Test  '
                 f'[{"PASSED" if passed else "FAILED"}]', fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(PLOT, 'test4_linearization.png'), dpi=150)
    plt.close()

    return passed, slope0, slope1


if __name__ == '__main__':
    passed, *_ = run()
    sys.exit(0 if passed else 1)
