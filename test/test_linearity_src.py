"""
Test 2: Modeling operator linearity with respect to source.

A linear forward operator F satisfies:
  (a) Amplitude scaling:  F(α·s) = α·F(s)
  (b) Superposition:      F(s₁ + s₂) = F(s₁) + F(s₂)

Sub-test 2a – amplitude scaling:
  Run with source amplitude 1 → d₁
  Run with source amplitude 2 → d₂
  Check: ||d₂ − 2·d₁|| / ||d₁|| < ε_mach  (machine precision)

Sub-test 2b – superposition of two sources:
  Run with source A only  → d_A
  Run with source B only  → d_B
  Run with sources A+B    → d_AB
  Check: ||d_AB − (d_A + d_B)|| / ||d_A|| < ε_mach

Both should hold to floating-point precision (round-off only).

Pass criterion: relative error < 1e-5.
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

WORK  = os.path.join(HERE, 'work', os.path.splitext(os.path.basename(__file__))[0])
PLOT  = os.path.join(HERE, 'plots')
os.makedirs(WORK, exist_ok=True)
os.makedirs(PLOT, exist_ok=True)

NZ, NX = 101, 201
DZ, DX = 10.0, 10.0
VP      = 2500.0
RHO     = 1.0
F0      = 15.0
DT      = 1.0e-3
TMAX    = 1.5

SX_A = NX // 2 * DX            # source A at horizontal center
SX_B = (NX // 2 + 20) * DX     # source B shifted by 200 m
SZ   = NZ // 2 * DZ
RECVS = [(DX * (NX // 4 + 10 * k), SZ) for k in range(8)]

PASS_TOL = 1e-5


def _base_params(tag, ns, geom_dir, synth_dir, extra=None):
    p = {
        'nx': NX, 'nz': NZ, 'dx': DX, 'dz': DZ,
        'dt': DT, 'data_dt': DT, 'tmax': TMAX,
        'ns': ns,
        'file_geometry': f'./{geom_dir}/geometry.txt',
        'which_medium': 'acoustic-iso',
        'model_name': 'vp, rho',
        'file_vp':  './model/vp.bin',
        'file_rho': './model/rho.bin',
        'dir_synthetic': f'./{synth_dir}',
        'verbose': 'n',
    }
    if extra:
        p.update(extra)
    return p


def run():
    # Write constant model once
    vp  = np.full((NZ, NX), VP)
    rho = np.full((NZ, NX), RHO)
    write_model(os.path.join(WORK, 'model', 'vp.bin'),  vp)
    write_model(os.path.join(WORK, 'model', 'rho.bin'), rho)

    results = {}

    # ─ Sub-test 2a: amplitude scaling ─────────────────────────────────────────
    for amp, tag in [(1.0, 'amp1'), (2.0, 'amp2')]:
        write_geometry(
            os.path.join(WORK, f'geom_{tag}'),
            [Source(x=SX_A, z=SZ, f0=F0, amp=amp)],
            RECVS
        )
        write_param(
            os.path.join(WORK, f'param_{tag}.rb'),
            _base_params(tag, 1, f'geom_{tag}', f'syn_{tag}')
        )
        shutil.rmtree(os.path.join(WORK, f'syn_{tag}'), ignore_errors=True)
        run_owl('owl_modeling2', f'param_{tag}.rb', WORK)
        results[tag], _ = read_su(
            os.path.join(WORK, f'syn_{tag}', 'shot_1_seismogram_p.su'))

    d1, d2 = results['amp1'].astype(np.float64), results['amp2'].astype(np.float64)
    err_2a = rel_err(d2, 2.0 * d1)
    pass_2a = err_2a < PASS_TOL
    report('Test 2a – amplitude scaling', pass_2a,
           f'||d₂ − 2·d₁|| / ||2·d₁|| = {err_2a:.2e}')

    # ─ Sub-test 2b: superposition ──────────────────────────────────────────────
    # Source A only (shot 1), Source B only (shot 1 with different geometry)
    for tag, sx in [('A', SX_A), ('B', SX_B)]:
        write_geometry(
            os.path.join(WORK, f'geom_{tag}'),
            [Source(x=sx, z=SZ, f0=F0, amp=1.0)],
            RECVS
        )
        write_param(
            os.path.join(WORK, f'param_{tag}.rb'),
            _base_params(tag, 1, f'geom_{tag}', f'syn_{tag}')
        )
        shutil.rmtree(os.path.join(WORK, f'syn_{tag}'), ignore_errors=True)
        run_owl('owl_modeling2', f'param_{tag}.rb', WORK)
        results[tag], _ = read_su(
            os.path.join(WORK, f'syn_{tag}', 'shot_1_seismogram_p.su'))

    # A+B simultaneously: use 2-source geometry for a single shot
    write_geometry(
        os.path.join(WORK, 'geom_AB'),
        [[Source(x=SX_A, z=SZ, f0=F0, amp=1.0),
          Source(x=SX_B, z=SZ, f0=F0, amp=1.0)]],
        RECVS
    )
    write_param(
        os.path.join(WORK, 'param_AB.rb'),
        _base_params('AB', 1, 'geom_AB', 'syn_AB')
    )
    shutil.rmtree(os.path.join(WORK, 'syn_AB'), ignore_errors=True)
    run_owl('owl_modeling2', 'param_AB.rb', WORK)
    results['AB'], _ = read_su(
        os.path.join(WORK, 'syn_AB', 'shot_1_seismogram_p.su'))

    dA  = results['A'].astype(np.float64)
    dB  = results['B'].astype(np.float64)
    dAB = results['AB'].astype(np.float64)
    err_2b = rel_err(dAB, dA + dB)
    pass_2b = err_2b < PASS_TOL
    report('Test 2b – superposition', pass_2b,
           f'||d_AB − (d_A+d_B)|| / ||d_A+d_B|| = {err_2b:.2e}')

    passed = pass_2a and pass_2b

    # ── Visualization ──────────────────────────────────────────────────────────
    t_ax = np.arange(d1.shape[0]) * DT

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    ax = axes[0]
    ax.plot(t_ax, d1[:, 0], 'b-',  lw=1.2, label='amp=1')
    ax.plot(t_ax, d2[:, 0] / 2.0, 'r--', lw=1.2, label='amp=2 (÷2)')
    ax.set_title(f'2a: Amplitude scaling\nTrace 1  err={err_2a:.1e}  '
                 f'[{"PASS" if pass_2a else "FAIL"}]')
    ax.set_xlabel('Time (s)'); ax.set_ylabel('Pressure')
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.plot(t_ax, dA[:, 0], 'b-', lw=1.2, label='Source A')
    ax.plot(t_ax, dB[:, 0], 'g-', lw=1.2, label='Source B')
    ax.plot(t_ax, dAB[:, 0], 'r--', lw=1.2, label='A+B (simultaneous)')
    ax.plot(t_ax, dA[:, 0] + dB[:, 0], 'k:', lw=1.5, label='A + B (summed)')
    ax.set_title(f'2b: Superposition\nTrace 1  err={err_2b:.1e}  '
                 f'[{"PASS" if pass_2b else "FAIL"}]')
    ax.set_xlabel('Time (s)'); ax.legend(fontsize=7)

    ax = axes[2]
    diff_2a = (d2 - 2.0 * d1).ravel()
    diff_2b = (dAB - (dA + dB)).ravel()
    ax.semilogy(np.abs(diff_2a) + 1e-30, 'b.', ms=2, alpha=0.4, label='2a residual')
    ax.semilogy(np.abs(diff_2b) + 1e-30, 'r.', ms=2, alpha=0.4, label='2b residual')
    ax.set_xlabel('Sample index'); ax.set_ylabel('|residual|')
    ax.set_title('Residual magnitude (all traces)')
    ax.legend(fontsize=8)

    fig.suptitle('Test 2: Modeling Operator Linearity w.r.t. Source  '
                 f'[{"PASSED" if passed else "FAILED"}]', fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(PLOT, 'test2_linearity_src.png'), dpi=150)
    plt.close()

    return passed, err_2a, err_2b


if __name__ == '__main__':
    passed, *_ = run()
    sys.exit(0 if passed else 1)
