"""
Test 3: Modeling operator adjoint test (source–receiver reciprocity).

For the acoustic wave equation in a medium with no absorption, the Green's
function is symmetric:
    G(xr, xs, ω) = G(xs, xr, ω)

i.e., swapping source and receiver yields the same recorded signal.  This
is the physical manifestation of the adjoint (transpose) of the forward
modeling operator being its own forward operator applied from the receiver
side.

Numerically, for each source–receiver pair (A→B) and its reciprocal (B→A):
    ||d(A→B) − d(B→A)|| / ||d(A→B)|| < ε_tol

We test several geometrically distinct pairs to stress-test the staggered-
grid interpolation.

Pass criterion: relative L2 error < 1e-3 (interpolation error from non-
integer grid positions).
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

NZ, NX = 201, 201
DZ, DX = 10.0, 10.0
VP   = 2500.0
RHO  = 1.0
F0   = 15.0
DT   = 1.0e-3
TMAX = 1.5
PASS_TOL = 1e-3   # sub-sample interpolation introduces ~0.1% error

# Source–receiver pairs to test (coordinates in metres, on-grid positions)
PAIRS = [
    ((500.0,  500.0), (1500.0, 1000.0)),
    ((300.0,  300.0), (1700.0,  700.0)),
    ((1000.0, 200.0), (1000.0, 1800.0)),
    ((600.0,  800.0), (1400.0,  400.0)),
]


def _run_shot(tag, sx, sz, rx, rz):
    geom_dir = f'geom_{tag}'
    syn_dir  = f'syn_{tag}'
    write_geometry(
        os.path.join(WORK, geom_dir),
        [Source(x=sx, z=sz, f0=F0, amp=1.0)],
        [(rx, rz)]
    )
    write_param(os.path.join(WORK, f'param_{tag}.rb'), {
        'nx': NX, 'nz': NZ, 'dx': DX, 'dz': DZ,
        'dt': DT, 'data_dt': DT, 'tmax': TMAX,
        'ns': 1,
        'file_geometry': f'./{geom_dir}/geometry.txt',
        'which_medium': 'acoustic-iso',
        'model_name': 'vp, rho',
        'file_vp':  './model/vp.bin',
        'file_rho': './model/rho.bin',
        'dir_synthetic': f'./{syn_dir}',
        'verbose': 'n',
    })
    shutil.rmtree(os.path.join(WORK, syn_dir), ignore_errors=True)
    run_owl('owl_modeling2', f'param_{tag}.rb', WORK)
    d, _ = read_su(os.path.join(WORK, syn_dir, 'shot_1_seismogram_p.su'))
    return d[:, 0].astype(np.float64)


def run():
    vp  = np.full((NZ, NX), VP)
    rho = np.full((NZ, NX), RHO)
    write_model(os.path.join(WORK, 'model', 'vp.bin'),  vp)
    write_model(os.path.join(WORK, 'model', 'rho.bin'), rho)

    errs = []
    traces_fwd, traces_rec = [], []

    for k, ((sx, sz), (rx, rz)) in enumerate(PAIRS):
        d_fwd = _run_shot(f'fwd{k}', sx, sz, rx, rz)
        d_rec = _run_shot(f'rec{k}', rx, rz, sx, sz)
        e = rel_err(d_fwd, d_rec)
        errs.append(e)
        traces_fwd.append(d_fwd)
        traces_rec.append(d_rec)
        print(f'  Pair {k+1}: src=({sx:.0f},{sz:.0f}) rec=({rx:.0f},{rz:.0f})  '
              f'err = {e:.2e}')

    max_err = max(errs)
    passed  = max_err < PASS_TOL
    report('Test 3 – Source reciprocity (modeling adjoint)', passed,
           f'max rel-L2 err = {max_err:.2e}')

    # ── Visualization ──────────────────────────────────────────────────────────
    t_ax = np.arange(traces_fwd[0].size) * DT
    fig, axes = plt.subplots(2, len(PAIRS), figsize=(14, 6))

    for k in range(len(PAIRS)):
        (sx, sz), (rx, rz) = PAIRS[k]
        ax_top = axes[0, k]
        ax_top.plot(t_ax, traces_fwd[k], 'b-',  lw=1.2, label=f'A→B')
        ax_top.plot(t_ax, traces_rec[k], 'r--', lw=1.2, label=f'B→A')
        ax_top.set_title(f'Pair {k+1}\nsrc=({sx:.0f},{sz:.0f})\nrec=({rx:.0f},{rz:.0f})',
                         fontsize=8)
        ax_top.set_xlabel('Time (s)')
        if k == 0:
            ax_top.set_ylabel('Pressure')
            ax_top.legend(fontsize=8)

        ax_bot = axes[1, k]
        ax_bot.plot(t_ax, traces_fwd[k] - traces_rec[k], 'k-', lw=1.0)
        ax_bot.set_title(f'Residual  err={errs[k]:.1e}\n'
                         f'[{"PASS" if errs[k]<PASS_TOL else "FAIL"}]', fontsize=8)
        ax_bot.set_xlabel('Time (s)')
        if k == 0:
            ax_bot.set_ylabel('d(A→B) − d(B→A)')

    fig.suptitle('Test 3: Source–Receiver Reciprocity (Modeling Operator Adjoint)  '
                 f'[{"PASSED" if passed else "FAILED"}]', fontsize=10)
    plt.tight_layout()
    fig.savefig(os.path.join(PLOT, 'test3_adjoint_src.png'), dpi=150)
    plt.close()

    return passed, errs


if __name__ == '__main__':
    passed, _ = run()
    sys.exit(0 if passed else 1)
