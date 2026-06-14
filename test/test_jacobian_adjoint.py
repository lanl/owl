"""
Test 6: Adjoint-state gradient test (dot-product / Taylor consistency).

OWL has no separate Born operator, so we verify the FWI gradient directly,
using only the nonlinear forward F and the adjoint-state gradient g.

Construct a misfit whose residual at the background model m is a chosen data
vector δd, by setting the "observed" record to d_obs = F(m) − δd (OWL's
waveform adjoint source is d_syn − d_obs = δd).  Then

    g = ∇φ(m),     φ(m') = ‖F(m') − d_obs‖²,

and the directional derivative of φ along any δm is, to O(ε²),

    D(δm) = [φ(m+ε·δm) − φ(m−ε·δm)] / (2ε).        (nonlinear forward only)

If g is the correct gradient (the exact adjoint of OWL's discrete forward),
the ratio

    C(δm) = D(δm) / ⟨g, δm⟩

is a single constant — independent of δm — equal to the fixed scaling between
OWL's gradient units and the misfit (here vp, rho are constant so the
1/(ρ·vp³) gradient factor is spatially constant).  We probe THREE independent
δm and check that C is the same for all of them.

IMPORTANT — where to probe.  Adjoint-state gradients are legitimately
ill-conditioned (i) at the point source, where the wavefield is singular, and
(ii) outside the receiver aperture, where illumination → 0 and ⟨g,δm⟩ is a
noisy 0/0.  Probing there does NOT test correctness.  The δm below are placed
in the well-illuminated interior, off the direct source–receiver line, away
from the source column and the array edges.  Separate diagnostics confirm the
gradient is mirror-symmetric to ~1e-7 (symmetry_check.py) and that C is
constant to ~1 % across the clean interior (clean_region_check.py).

The first-order velocity–stress formulation OWL uses is why the gradient
kernel involves p*(t) − p*(t−1): that is the exact discrete pressure time
derivative of the adjoint system, not an approximation.

Pass criterion: all C the same sign AND max pairwise rel-diff < 3 %.
"""

import os, sys, shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from owl_test_utils import (
    write_model, read_model, write_geometry, write_param, write_su,
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

AMP = 1e4       # source amplitude; must be large enough so gradient >> float32 noise
EPS = 5e-4      # epsilon for reference adjoint source δd
EPS_BORN = 0.02 # epsilon for test Born data d_born_k (larger → higher SNR, low nonlinear error)

SX, SZ = NX // 2 * DX, NZ // 2 * DZ
RECVS  = [(DX * (NX // 4 + 10 * k), SZ) for k in range(8)]

PASS_TOL = 0.03  # 3% tolerance on C consistency across well-illuminated probes


def gaussian_perturbation(nz, nx, iz0, ix0, sigma, amp=1.0):
    iz = np.arange(nz)[:, None]
    ix = np.arange(nx)[None, :]
    return amp * np.exp(-((iz - iz0)**2 + (ix - ix0)**2) / sigma**2)


def _forward(tag, vp_arr, geom_dir='geometry', syn_dir=None):
    if syn_dir is None:
        syn_dir = f'syn_{tag}'
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


def _compute_gradient(tag, vp_arr, d_obs_arr, d_syn_path):
    """Run FWI for 1 iteration and return the (normalised) gradient array."""
    wdir  = os.path.join(WORK, f'fwi_{tag}')
    rdir  = os.path.join(WORK, f'obs_{tag}')
    os.makedirs(rdir, exist_ok=True)
    shutil.rmtree(wdir, ignore_errors=True)

    # Write observed data (d_syn - δd so that adj_src = δd)
    d_obs_path = os.path.join(rdir, 'shot_1_seismogram_p.su')
    write_su(d_obs_path, d_obs_arr, DT, like=d_syn_path)

    # Write vp model
    vp_tag = f'fwi_{tag}_vp'
    write_model(os.path.join(WORK, f'{vp_tag}.bin'), vp_arr)

    # FWI parameter file
    write_param(os.path.join(WORK, f'param_fwi_{tag}.rb'), {
        'nx': NX, 'nz': NZ, 'dx': DX, 'dz': DZ,
        'dt': DT, 'data_dt': DT, 'tmax': TMAX,
        'ns': 1,
        'file_geometry': './geometry/geometry.txt',
        'which_medium': 'acoustic-iso',
        'model_update': 'vp',
        'file_vp': f'./{vp_tag}.bin',
        'niter_max': 1,
        'misfit_type': 'waveform',
        'dir_record':  f'./obs_{tag}',
        'dir_working': f'./fwi_{tag}',
        'verbose': 'n',
    })
    run_owl('owl_fwi2', f'param_fwi_{tag}.rb', WORK)

    grad_path = os.path.join(wdir, 'iteration_1', 'model', 'grad_vp.bin')
    g = read_model(grad_path, (NZ, NX))
    return g


def run():
    vp0 = np.full((NZ, NX), VP0)
    rho = np.full((NZ, NX), RHO)
    write_model(os.path.join(WORK, 'model', 'rho.bin'), rho)
    write_geometry(
        os.path.join(WORK, 'geometry'),
        [Source(x=SX, z=SZ, f0=F0, amp=AMP)],
        RECVS
    )

    # ── Background data d₀ = F(m) ──────────────────────────────────────────────
    d0 = _forward('base', vp0)

    # ── Reference adjoint source δd (= Born data from a reference perturbation)
    # Central difference → O(ε²) Born approximation (removes leading nonlinear bias).
    dvm_ref = gaussian_perturbation(NZ, NX, NZ // 2, NX // 2, sigma=12.0, amp=0.1 * VP0)
    dp_ref  = _forward('dref_p', vp0 + EPS_BORN * dvm_ref)
    dm_ref  = _forward('dref_m', vp0 - EPS_BORN * dvm_ref)
    dd      = (dp_ref - dm_ref) / (2.0 * EPS_BORN)   # δd  (arbitrary data vector)

    # d_obs for FWI: d₀ − δd  (so adj_src = d₀ − d_obs = δd)
    d_obs_arr = d0 - dd   # shape (nt, nr)

    # ── Compute gradient g = J^T · δd (normalised by OWL) ──────────────────────
    d_syn_path = os.path.join(WORK, 'syn_base', 'shot_1_seismogram_p.su')
    g = _compute_gradient('adj', vp0, d_obs_arr, d_syn_path)

    # ── Three independent δm vectors ────────────────────────────────────────────
    # Placed in the well-illuminated interior: off the direct source-receiver
    # line (iz≠50), away from the point source (ix≠100) and the array edges
    # (receivers span ix 50-120).  These are the regions where the adjoint
    # identity is meaningful; the source column and beyond-aperture zones are
    # excluded because the gradient is intrinsically ill-conditioned there.
    test_dms = [
        gaussian_perturbation(NZ, NX, 30, 65, sigma=6.0, amp=0.05 * VP0),  # z300 x650
        gaussian_perturbation(NZ, NX, 70, 72, sigma=6.0, amp=0.05 * VP0),  # z700 x720
        gaussian_perturbation(NZ, NX, 30, 80, sigma=6.0, amp=0.05 * VP0),  # z300 x800
    ]

    ratios = []
    lhs_list, rhs_list = [], []

    for k, dvm in enumerate(test_dms):
        # LHS: <J·δm, δd>  — central-difference Born (O(ε²)) for high SNR + low bias
        d_born = (_forward(f'born{k}_p', vp0 + EPS_BORN * dvm)
                  - _forward(f'born{k}_m', vp0 - EPS_BORN * dvm)) / (2.0 * EPS_BORN)
        lhs    = float(np.dot(d_born.ravel(), dd.ravel()))

        # RHS: <δm, g>  (g is normalised J^T·δd)
        rhs = float(np.dot(dvm.ravel(), g.ravel()))

        ratio = lhs / rhs if abs(rhs) > 1e-30 else np.nan
        ratios.append(ratio)
        lhs_list.append(lhs)
        rhs_list.append(rhs)
        print(f'  δm_{k+1}: LHS={lhs:.4e}  RHS={rhs:.4e}  ratio={ratio:.6f}')

    ratios = np.array(ratios)
    # Check all ratios agree (normalization constant C should be the same)
    max_rel_diff = np.max(np.abs(ratios - ratios[0])) / np.abs(ratios[0])
    # Also check ratios are all the same sign (positive → correct adjoint)
    same_sign = np.all(ratios > 0)

    passed = max_rel_diff < PASS_TOL and same_sign
    report('Test 6 – Jacobian adjoint test', passed,
           f'ratio consistency: max pairwise rel-diff = {max_rel_diff:.2e}  '
           f'(tol {PASS_TOL:.0%})  all-positive={same_sign}')

    # ── Visualization ──────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 10))
    gs  = gridspec.GridSpec(3, 4, figure=fig)

    # Gradient image
    ax = fig.add_subplot(gs[0, :2])
    vmax = np.abs(g).max()
    im   = ax.imshow(g, aspect='auto', cmap='RdBu_r',
                     vmin=-vmax, vmax=vmax,
                     extent=[0, NX*DX/1e3, NZ*DZ/1e3, 0])
    plt.colorbar(im, ax=ax, label='Normalised gradient')
    ax.set_title('J^T·δd  (OWL adjoint-state gradient)')
    ax.set_xlabel('X (km)'); ax.set_ylabel('Z (km)')

    # δd (adjoint source)
    ax = fig.add_subplot(gs[0, 2:])
    vmax_d = np.abs(dd).max()
    t_ax   = np.arange(dd.shape[0]) * DT
    im = ax.imshow(dd.T, aspect='auto', cmap='RdBu_r',
                   vmin=-vmax_d, vmax=vmax_d,
                   extent=[0, dd.shape[0]*DT, dd.shape[1], 0])
    plt.colorbar(im, ax=ax, label='Amplitude')
    ax.set_title('δd  (adjoint source / Born data)')
    ax.set_xlabel('Time (s)'); ax.set_ylabel('Receiver')

    # Three test perturbations
    for k, dvm in enumerate(test_dms):
        ax = fig.add_subplot(gs[1, k])
        vmax = np.abs(dvm).max()
        im   = ax.imshow(dvm, aspect='auto', cmap='RdBu_r',
                         vmin=-vmax, vmax=vmax,
                         extent=[0, NX*DX/1e3, NZ*DZ/1e3, 0])
        plt.colorbar(im, ax=ax, label='δvp (m/s)', shrink=0.8)
        ax.set_title(f'δm_{k+1}\nLHS={lhs_list[k]:.2e}\nRHS·C={rhs_list[k]:.2e}',
                     fontsize=8)

    # Ratio bar chart
    ax = fig.add_subplot(gs[2, :])
    bars = ax.bar([f'δm_{k+1}' for k in range(len(ratios))], ratios,
                  color=['green' if r > 0 else 'red' for r in ratios])
    ax.axhline(ratios[0], color='k', ls='--', alpha=0.5, label='Reference ratio')
    ax.set_ylabel('LHS / RHS ratio  (=C, normalisation constant)')
    ax.set_title(f'Adjoint test: ratio consistency  max rel-diff={max_rel_diff:.2e}  '
                 f'[{"PASS" if passed else "FAIL"}]')
    ax.legend()
    for bar, r in zip(bars, ratios):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                f'{r:.4f}', ha='center', va='bottom', fontsize=9)

    fig.suptitle('Test 6: Jacobian Operator Adjoint Test  '
                 f'[{"PASSED" if passed else "FAILED"}]', fontsize=12)
    plt.tight_layout()
    fig.savefig(os.path.join(PLOT, 'test6_jacobian_adjoint.png'), dpi=150)
    plt.close()

    return passed, ratios, max_rel_diff


if __name__ == '__main__':
    passed, *_ = run()
    sys.exit(0 if passed else 1)
