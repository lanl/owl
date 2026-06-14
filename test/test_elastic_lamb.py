"""
Test 7: 2D elastic modeling vs. the Lamb problem on a TILTED free surface.

A homogeneous isotropic elastic half-space with a *planar* free surface tilted
by angle theta is exactly Lamb's problem rotated by theta.  With a point force
applied along the surface normal and receivers along the surface, the response
— expressed in surface-normal / surface-tangential particle velocity — must be
identical to the flat (theta=0) Lamb solution for EVERY theta.  Any
theta-dependence exposes an error in OWL's topographic free-surface scheme.

Reference: native discrete-wavenumber half-space solution (lamb_reference.py).
owl side: owl_modeling2, which_medium=elastic-tti / anisotropy_type=iso, an
immersed free surface from a tilted topography file with near-surface dz
refinement (free_surface_dz_refine), source/receiver depths measured from the
surface so they sit exactly on it.

Geometry: surface deepens downdip, z_surf(x) = z0 + (x-xc)*tan(theta).  Force is
along the inward normal (OWL polar = -theta).  Receivers at fixed along-surface
offsets r, so x_r = xc + r*cos(theta) and the SAME reference serves all angles.
Rotation of owl velocities into the surface frame:
    v_normal = -sin(theta)*vx + cos(theta)*vz
    v_tang   =  cos(theta)*vx + sin(theta)*vz

Pass criterion: for every tilt angle and both components, the normalised
waveform misfit vs. the (single) reference stays small and does NOT grow with
theta.
"""

import os, sys, shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from owl_test_utils import (write_model, read_model, write_param, write_geometry,
                            run_owl, read_su, Source)
from lamb_reference import lamb_halfspace_velocity, rayleigh_speed

PLOT = os.path.join(HERE, 'plots')
os.makedirs(PLOT, exist_ok=True)

# ── medium / acquisition ────────────────────────────────────────────────────
VP, VS, RHO = 3000.0, 1732.0, 2000.0     # Poisson solid (vp/vs = sqrt3)
F0   = 12.5
DX = DZ = 10.0
DT   = 2.0e-4
DATA_DT = 1.0e-3
TMAX = 0.6
REFINE = 4.0
AMP  = 1.0e6

NX, NZ = 121, 121                         # 1200 m x 1200 m
XC   = 250.0                              # source position (updip side of spread)
EXT  = 200.0                              # topo extension beyond the domain (m)
OFFSETS = np.array([150., 250., 350., 450., 550.])   # along-surface, downdip (m)

ANGLES = [0.0, 5.0, 10.0, 20.0, 30.0]
POLAR_SIGN = -1.0     # descending-downdip surface: slopex=-dz/dx on the receiver
                      # side, so the inward surface normal is (-sin th, cos th)
SRC_DEPTH = 0.0       # source depth below the free surface (m)
REC_DEPTH = 0.0       # receiver depth below the free surface (m)
MECH = 'force'        # 'force' (surface-normal) or 'explosion' (isotropic)

SNAPS = '0.0, 0.05, 0.6'   # OWL snapshot times: 0, 0.05, ..., 0.6 s (index 1..13)
SNAP_L = 6                 # snapshot index to draw; t = (SNAP_L-1)*0.05 s


def write_topo(path, theta_deg):
    """Continuous planar free surface, anchored at elevation z = 0 at x = 0 and
    descending downdip (toward the receivers) with slope tan(theta).  It extends
    LINEARLY beyond the model on both sides -- a true tilted plane, no flat
    clamp.  Flat case (theta=0) -> z = 0 everywhere."""
    th = np.radians(theta_deg)
    xmax = (NX - 1) * DX
    x = np.arange(-EXT, xmax + EXT + DX, DX)
    z = -x * np.tan(th)                          # continuous plane, z = 0 at x = 0
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        for xi, zi in zip(x, z):
            f.write(f"{xi:.4f} {zi:.6f}\n")


def run_owl_angle(theta_deg):
    """Run one tilted-surface elastic simulation; return (vx, vz, nt)."""
    tag  = f"lamb_t{int(round(theta_deg))}"
    work = os.path.join(HERE, 'work', tag)
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work, exist_ok=True)
    th = np.radians(theta_deg)

    # constant half-space models
    for nm, val in [('vp', VP), ('vs', VS), ('rho', RHO)]:
        write_model(os.path.join(work, 'model', f'{nm}.bin'),
                    np.full((NZ, NX), val))

    write_topo(os.path.join(work, 'model', 'topo.txt'), theta_deg)

    # source on/below surface.  MECH selects 'force ...' (normal force) or
    # 'explosion' (isotropic, for the localisation test).
    mech = (f"force {POLAR_SIGN * theta_deg:.6f} 0.0" if MECH == 'force'
            else 'explosion')
    src = Source(x=XC, z=SRC_DEPTH, f0=F0, amp=AMP, mechanism=mech)
    recvs = [(XC + r * np.cos(th), REC_DEPTH) for r in OFFSETS]
    write_geometry(os.path.join(work, 'geometry'), [src], recvs)

    write_param(os.path.join(work, 'param.rb'), {
        'nx': NX, 'nz': NZ, 'dx': DX, 'dz': DZ,
        'dt': DT, 'data_dt': DATA_DT, 'tmax': TMAX, 'ns': 1,
        'file_geometry': './geometry/geometry.txt',
        'which_medium': 'elastic-tti', 'anisotropy_type': 'iso',
        'model_name': 'vp, vs, rho',
        'file_vp': './model/vp.bin', 'file_vs': './model/vs.bin',
        'file_rho': './model/rho.bin',
        'yn_free_surface': 'y', 'file_topo': './model/topo.txt',
        'free_surface_dz_refine': REFINE,
        'measure_source_depth_from_surface': 'y',
        'measure_receiver_depth_from_surface': 'y',
        'dir_synthetic': './syn', 'verbose': 'n',
        'snaps': SNAPS, 'dir_snapshot': './snap',
    })
    run_owl('owl_modeling2', 'param.rb', work)
    vx, _ = read_su(os.path.join(work, 'syn', 'shot_1_seismogram_x.su'))
    vz, _ = read_su(os.path.join(work, 'syn', 'shot_1_seismogram_z.su'))
    return vx.astype(np.float64), vz.astype(np.float64)


def read_snapshot(theta_deg, comp, l):
    """Read an OWL regular-mesh wavefield snapshot (NZ x NX); z=0 row is at the
    highest topography elevation, the air above the surface is zero."""
    tag = f"lamb_t{int(round(theta_deg))}"
    p = os.path.join(HERE, 'work', tag, 'snap',
                     f'shot_1_forward_wavefield_{comp}_{l}.bin')
    return read_model(p, (NZ, NX))


def best_scale_misfit(a, ref):
    """Global least-squares scale (allowing sign) then normalised L2 misfit."""
    s = np.dot(a.ravel(), ref.ravel()) / np.dot(a.ravel(), a.ravel())
    mis = np.linalg.norm(s * a - ref) / np.linalg.norm(ref)
    return s, mis


def shape_misfit(a, b):
    """Per-trace shape misfit: amplitude-normalise each trace (sign allowed),
    then mean normalised L2.  Measures waveform SHAPE, not absolute level."""
    m = []
    for i in range(a.shape[1]):
        ai, bi = a[:, i], b[:, i]
        na = np.linalg.norm(ai)
        if na == 0:
            continue
        s = np.dot(ai, bi) / np.dot(ai, ai)
        m.append(np.linalg.norm(s * ai - bi) / np.linalg.norm(bi))
    return float(np.mean(m))


def run():
    nt = int(round(TMAX / DATA_DT)) + 1
    # single reference (flat Lamb, surface-aligned frame), rotated per-angle to
    # OWL's (x,z) lab frame.  Surface normal (force) = (sin th, cos th);
    # along-surface +r tangent = (cos th, -sin th)  (elevation-up convention).
    ref_n, ref_t = lamb_halfspace_velocity(OFFSETS, nt, DATA_DT, VP, VS, RHO, F0)
    vR = rayleigh_speed(VP, VS)

    def fit(a):
        b = np.zeros((nt, a.shape[1]))
        b[:min(a.shape[0], nt)] = a[:nt]
        return b

    owl = {}
    for th in ANGLES:
        vx, vz = run_owl_angle(th)
        owl[th] = (fit(vx), fit(vz))

    # one global calibration scale K from the flat case (no per-trace, no per-
    # angle normalisation): the solutions are then compared as-is.
    vx0, vz0 = owl[0.0]
    a = np.concatenate([vx0.ravel(), vz0.ravel()])
    b = np.concatenate([ref_t.ravel(), ref_n.ravel()])   # th=0: vx=tang, vz=norm
    K = np.dot(a, b) / np.dot(a, a)

    results = {}
    for th in ANGLES:
        s, c = np.sin(np.radians(th)), np.cos(np.radians(th))
        # descending-downdip surface: normal=(-sin,cos), tangent(+r)=(cos,sin)
        ref_vx = c * ref_t - s * ref_n
        ref_vz = s * ref_t + c * ref_n
        vx, vz = owl[th]        
        owlx, owlz = K * vx, K * vz                # as-is, single global scale
        # both reported as RELATIVE ERRORS (0 = perfect):
        #   shape err = per-trace amplitude-normalised waveform L2 error
        #   amp   err = |OWL/ref energy ratio - 1|
        sh = 0.5 * (shape_misfit(vx, ref_vx) + shape_misfit(vz, ref_vz))
        eo = np.sqrt((owlx ** 2 + owlz ** 2).sum())
        er = np.sqrt((ref_vx ** 2 + ref_vz ** 2).sum())
        amperr = abs(eo / er - 1.0)
        results[th] = dict(vx=owlx, vz=owlz, rvx=ref_vx, rvz=ref_vz,
                           shape=sh, amperr=amperr)
        print(f'  tilt angle={th:4.0f} deg :  shape err={sh:.3f}   amp err={amperr:.3f}')

    worst = max(max(r['shape'], r['amperr']) for r in results.values())
    passed = worst < 0.20
    print(f'\n  worst-case error (shape or amplitude) = {worst:.3f}')
    print(f'  {"PASSED" if passed else "FAILED"} (tol 0.20)')

    snaps = {th: read_snapshot(th, 'z', SNAP_L) for th in ANGLES}
    _plot(results, nt, vR, passed, snaps)
    return passed, results


def _plot(results, nt, vR, passed, snaps):
    t = np.arange(nt) * DATA_DT
    allref = np.concatenate([np.abs(r['rvz']) for r in results.values()])
    sc = 80.0 / allref.max()
    ncol = len(ANGLES)
    xext, zext = (NX - 1) * DX, (NZ - 1) * DZ

    fig = plt.figure(figsize=(4.2 * ncol, 12))
    gs = gridspec.GridSpec(3, ncol, height_ratios=[1.25, 1, 1],
                           hspace=0.38, wspace=0.28)

    # ── top row: OWL wavefield snapshot per tilt (shows each tilted surface) ──
    for j, thd in enumerate(ANGLES):
        axs = fig.add_subplot(gs[0, j])
        th = np.radians(thd)
        snap = snaps[thd]
        clip = 0.25 * (np.percentile(np.abs(snap), 99.5) or 1.0)
        # Surface depth below the regular-grid peak at the source column (the
        # masked air above the surface is exactly zero), then recover the peak
        # elevation above the z=0 datum so we can plot depth-below-z=0.
        scol = np.abs(snap[:, int(round(XC / DX))])
        sd = np.argmax(scol > 1e-3 * (scol.max() or 1.0)) * DZ
        peak = sd - XC * np.tan(th)              # peak elevation above z=0 datum
        axs.imshow(snap, aspect='auto', cmap='bwr', vmin=-clip, vmax=clip,
                   extent=[0, xext, zext - peak, -peak])
        # free surface in depth-below-z=0: z(x)=-x*tan(th) -> depth = x*tan(th)
        xs = np.linspace(0, xext, 400)
        axs.plot(xs, xs * np.tan(th), 'k-', lw=1.5,
                 label='free surface' if j == 0 else None)
        axs.plot(XC, XC * np.tan(th), 'y*', ms=13, mec='k',
                 label='source' if j == 0 else None)
        rxs = XC + OFFSETS * np.cos(th)
        axs.plot(rxs, rxs * np.tan(th), 'gv', ms=6, mec='k',
                 label='receivers' if j == 0 else None)
        axs.set_xlim(0, xext); axs.set_ylim(1200.0, -100.0)
        axs.set_title(f'tilt={thd:.0f}°  v_z snapshot @ t={(SNAP_L-1)*0.05:.2f}s',
                      fontsize=8)
        if j == 0:
            axs.set_ylabel('depth below z=0 (m)')
            axs.legend(fontsize=6, loc='lower right', framealpha=0.9)
        axs.set_xlabel('x (m)', fontsize=8)

    # ── waveform overlays: OWL vs analytic Lamb ──────────────────────────────
    for j, thd in enumerate(ANGLES):
        for row, key, rkey, lab in [(1, 'vz', 'rvz', 'vertical v_z'),
                                    (2, 'vx', 'rvx', 'horizontal v_x')]:
            ax = fig.add_subplot(gs[row, j])
            V, ref = results[thd][key], results[thd][rkey]
            for i, r in enumerate(OFFSETS):
                ax.plot(t, r + sc * ref[:, i], 'b-', lw=1.3,
                        label='reference (Lamb)' if i == 0 else None)
                ax.plot(t, r + sc * V[:, i], 'r--', lw=1.0,
                        label='OWL (1 global scale)' if i == 0 else None)
            ax.plot(OFFSETS / vR, OFFSETS, 'g:', lw=1)
            ax.set_title(f'tilt={thd:.0f}°  {lab}\nshape err={results[thd]["shape"]:.2f}  '
                         f'amp err={results[thd]["amperr"]:.2f}', fontsize=8)
            if j == 0:
                ax.set_ylabel('offset (m)')
            if row == 2:
                ax.set_xlabel('time (s)')
            if row == 1 and j == 0:
                ax.legend(fontsize=6, loc='upper right')
            ax.set_xlim(0, TMAX); ax.invert_yaxis()

    fig.suptitle('Test 7: 2D elastic Lamb on a tilted free surface   '
                 f'[shape {"PASS" if passed else "FAIL"}]   '
                 'OWL (one global scale, as-is) vs analytic Lamb reference   '
                 f'vp={VP} vs={VS} rho={RHO} f0={F0}Hz', fontsize=11)
    fig.savefig(os.path.join(PLOT, 'test7_elastic_lamb_tilt.png'),
                dpi=130, bbox_inches='tight')
    plt.close()
    print('figure ->', os.path.join(PLOT, 'test7_elastic_lamb_tilt.png'))


if __name__ == '__main__':
    passed, _ = run()
    sys.exit(0 if passed else 1)
