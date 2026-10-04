"""
Test 9: free-surface modeling with sources and receivers at and near the free surface.

Elastic (2D): owl_modeling2 with elastic-vhtiort (SSG) and elastic-tti (FSG), a flat free surface
with the default near-surface mesh refinement, and an explosion, a vertical force and a general
moment tensor at depths 0 m (on the surface), 1.25 m (between the first two refined rows) and 20 m.
The particle velocities at the surface and at 150 m depth are compared with the analytical
half-space solution (halfspace_reference.py) in absolute amplitude, after one time shift fitted per
solver and mechanism on the 20 m source.

Acoustic (2D): with a free surface (p = 0), the pressure is exactly the whole-space field of the
source minus that of its mirror image. owl_modeling2 with acoustic-iso, sources 5 m and 15 m below
the surface (off the grid) and receivers 15 m and 50 m below it (the stencils of the former reach
above the surface) is compared with two whole-space (all-PML) runs.

Pass criteria: elastic amplitude ratios within 12% at the surface receivers and within 6% at
150 m, waveform residuals below 10%, and acoustic relative errors below 1e-4. The free-surface
treatment is first-order accurate near the surface; a source exactly on the surface is the worst
case (with the default free_surface_dz_refine = 4, up to about 11% at the surface receivers and 5%
at depth for the FSG solver), and the error halves with each doubling of free_surface_dz_refine.
"""

import os
import sys
import shutil
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from owl_test_utils import write_model, read_su, write_param, run_owl, report
from halfspace_reference import halfspace_velocity

WORK = os.path.join(HERE, 'work', 'free_surface')
PLOT = os.path.join(HERE, 'plots')

# Elastic half-space (Poisson solid) and acquisition
VP, VS, RHO = 3000.0, 1732.0, 2000.0
F0, DT, TMAX = 12.5, 4.0e-4, 0.5
NX, NZ, D, PML = 121, 51, 10.0, 20
XS = 300.0
OFFSETS = np.arange(200.0, 601.0, 100.0)
REC_Z = [0.0, 150.0]
DEPTHS = [0.0, 1.25, 20.0]
# geometry-file mechanism, and the same source for the reference (owl's force amplitude is a force
# density on one grid cell of the regular mesh, i.e. a line force of dx*dz)
MECHS = {'explosion': ('explosion', dict(mxx=1.0, mzz=1.0)),
         'force': ('force 0 0', dict(fz=1.0, cell=D*D)),
         'moment tensor': ('mt 1 0 -0.5 0 0.7 0', dict(mxx=1.0, mzz=-0.5, mxz=0.7))}
SOLVERS = ['elastic-vhtiort', 'elastic-tti']
TOL_SURFACE, TOL_DEEP, TOL_SHAPE = 0.12, 0.06, 0.10
NT = int(round(TMAX/DT)) + 1

# Acoustic whole space and image test
AC_VP, AC_RHO, AC_F0, AC_DT, AC_TMAX = 2000.0, 1000.0, 10.0, 1.0e-3, 0.8
AC_NX, AC_XS, AC_SHIFT = 101, 500.0, 400.0
AC_DEPTHS = [5.0, 15.0]
AC_OFFSETS = np.arange(0.0, 401.0, 100.0)
AC_REC_Z = [15.0, 50.0]
TOL_ACOUSTIC = 1.0e-4


def elastic_run(solver, mech):
    """One owl run with a shot per source depth; returns {depth: (nt, 2*nrec) [vx, vz]}."""
    d = os.path.join(WORK, solver, mech.replace(' ', '_'))
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, 'geometry'))
    for name, v in (('vp', VP), ('vs', VS), ('rho', RHO)):
        write_model(os.path.join(d, 'model', f'{name}.bin'), np.full((NZ, NX), v))
    with open(os.path.join(d, 'geometry', 'geometry.txt'), 'w') as f:
        for i in range(len(DEPTHS)):
            f.write(f'shot_{i + 1}_geometry.txt\n')
    for i, h in enumerate(DEPTHS):
        with open(os.path.join(d, 'geometry', f'shot_{i + 1}_geometry.txt'), 'w') as f:
            f.write(f'{i + 1}\n1\n{XS} 0.0 {h}\n{MECHS[mech][0]}\nricker {F0} 1.0 0.0\n0 0\n')
            f.write(f'{len(OFFSETS)*len(REC_Z)}\n')
            for rz in REC_Z:
                for r in OFFSETS:
                    f.write(f'{XS + r} 0.0 {rz} 1.0\n')
    write_param(os.path.join(d, 'param.rb'), {
        'nx': NX, 'nz': NZ, 'dx': D, 'dz': D, 'dt': DT, 'tmax': TMAX, 'ns': len(DEPTHS),
        'file_geometry': './geometry/geometry.txt', 'which_medium': solver,
        'model_name': 'vp, vs, rho', 'file_vp': './model/vp.bin', 'file_vs': './model/vs.bin',
        'file_rho': './model/rho.bin', 'npml': PML, 'yn_free_surface': 'y', 'dir_synthetic': 'data',
        'verbose': 'n'})
    run_owl('owl_modeling2', 'param.rb', d)
    out = {}
    for i, h in enumerate(DEPTHS):
        out[h] = np.concatenate([read_su(os.path.join(d, 'data', f'shot_{i + 1}_seismogram_{c}.su'))[0]
                                 for c in 'xz'], axis=1)
    return out


def shift(a, tau, dt):
    """Delay the columns of a by tau seconds (fractional, via FFT)."""
    n = a.shape[0]
    f = np.fft.rfftfreq(2*n, dt)
    return np.fft.irfft(np.fft.rfft(a, n=2*n, axis=0)*np.exp(-2j*np.pi*f*tau)[:, None], n=2*n, axis=0)[:n]


def fit(o, r):
    """Least-squares amplitude ratio of o to r, and the residual after scaling."""
    c = np.sum(o*r)/np.sum(r*r)
    return c, np.linalg.norm(o - c*r)/np.linalg.norm(c*r)


def acoustic_run(tag, nz, fs, src_z, rec_dz):
    d = os.path.join(WORK, 'acoustic', tag)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, 'geometry'))
    write_model(os.path.join(d, 'model', 'vp.bin'), np.full((nz, AC_NX), AC_VP))
    write_model(os.path.join(d, 'model', 'rho.bin'), np.full((nz, AC_NX), AC_RHO))
    with open(os.path.join(d, 'geometry', 'geometry.txt'), 'w') as f:
        f.write('shot_1_geometry.txt\n')
    with open(os.path.join(d, 'geometry', 'shot_1_geometry.txt'), 'w') as f:
        f.write(f'1\n1\n{AC_XS} 0.0 {src_z}\nexplosion\nricker {AC_F0} 1.0 0.0\n0 0\n'
                f'{len(AC_OFFSETS)*len(AC_REC_Z)}\n')
        for rz in AC_REC_Z:
            for r in AC_OFFSETS:
                f.write(f'{AC_XS + r} 0.0 {rz + rec_dz} 1.0\n')
    write_param(os.path.join(d, 'param.rb'), {
        'nx': AC_NX, 'nz': nz, 'dx': D, 'dz': D, 'dt': AC_DT, 'tmax': AC_TMAX, 'ns': 1,
        'file_geometry': './geometry/geometry.txt', 'which_medium': 'acoustic-iso',
        'model_name': 'vp, rho', 'file_vp': './model/vp.bin', 'file_rho': './model/rho.bin',
        'npml': 30, 'yn_free_surface': 'y' if fs else 'n', 'dir_synthetic': 'data', 'verbose': 'n'})
    run_owl('owl_modeling2', 'param.rb', d)
    return read_su(os.path.join(d, 'data', 'shot_1_seismogram_p.su'))[0]


def run():
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(PLOT, exist_ok=True)
    n, nz = len(OFFSETS), len(REC_Z)
    cols = lambda ic, iz: slice((ic*nz + iz)*n, (ic*nz + iz + 1)*n)
    xr = np.concatenate([OFFSETS for _ in REC_Z])
    zr = np.concatenate([np.full(n, z) for z in REC_Z])

    # Elastic: OWL/analytical amplitude ratio per solver, mechanism, depth, receiver depth, component
    ratios = {}
    worst = dict(surface=0.0, deep=0.0, shape=0.0)
    for mech, (_, kw) in MECHS.items():
        ref = {}
        for h in DEPTHS:
            vx, vz = halfspace_velocity(xr, zr, h, NT, DT, VP, VS, RHO, F0, **kw)
            ref[h] = np.concatenate([vx, vz], axis=1)
        for solver in SOLVERS:
            owl = elastic_run(solver, mech)
            taus = np.linspace(-2*DT, 2*DT, 41)
            tau = taus[int(np.argmin([fit(owl[20.0], shift(ref[20.0], s, DT))[1] for s in taus]))]
            print(f'  {solver}, {mech}: OWL/analytical amplitude ratio [vx, vz] at the surface | at 150 m')
            for h in DEPTHS:
                rs = shift(ref[h], tau, DT)
                line = f'    source at {h:5.2f} m:'
                for iz, rz in enumerate(REC_Z):
                    for ic in range(2):
                        c, e = fit(owl[h][:, cols(ic, iz)], rs[:, cols(ic, iz)])
                        ratios[(solver, mech, h, iz, ic)] = c
                        key = 'surface' if rz == 0.0 else 'deep'
                        worst[key] = max(worst[key], abs(c - 1.0))
                        worst['shape'] = max(worst['shape'], e)
                        line += f'  {c:6.3f}'
                    line += '  |' if iz == 0 else ''
                print(line)

    # Acoustic: free surface vs direct minus image (exact)
    ac_err = {}
    for h in AC_DEPTHS:
        fs = acoustic_run(f'fs_{h:g}', 51, True, h, 0.0)
        direct = acoustic_run(f'direct_{h:g}', 101, False, AC_SHIFT + h, AC_SHIFT)
        image = acoustic_run(f'image_{h:g}', 101, False, AC_SHIFT - h, AC_SHIFT)
        ac_err[h] = np.linalg.norm(fs - (direct - image))/np.linalg.norm(direct - image)
        print(f'  acoustic-iso, source at {h:g} m: relative error vs direct - image = {ac_err[h]:.2e}')

    passed = (worst['surface'] <= TOL_SURFACE and worst['deep'] <= TOL_DEEP
              and worst['shape'] <= TOL_SHAPE and max(ac_err.values()) <= TOL_ACOUSTIC)

    # Figure: amplitude ratios at the two receiver depths
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for iz, (ax, tol) in enumerate(zip(axs, (TOL_SURFACE, TOL_DEEP))):
        ax.axhspan(1 - tol, 1 + tol, color='0.9')
        for solver, ls in zip(SOLVERS, ('-', '--')):
            for mech, color in zip(MECHS, ('C0', 'C1', 'C2')):
                for ic, marker in enumerate(('o', 's')):
                    ax.plot(range(len(DEPTHS)), [ratios[(solver, mech, h, iz, ic)] for h in DEPTHS], ls,
                            color=color, marker=marker, label=f'{solver}, {mech}, v{"xz"[ic]}')
        ax.set_xticks(range(len(DEPTHS)))
        ax.set_xticklabels([f'{h:g}' for h in DEPTHS])
        ax.set_xlabel('Source depth (m)')
        ax.set_title(f'Receivers at {REC_Z[iz]:g} m')
    axs[0].set_ylabel('OWL / analytical')
    axs[1].legend(fontsize=7, loc='lower right')
    fig.suptitle(f'Test 9: free surface [{"PASSED" if passed else "FAILED"}]')
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT, 'test9_free_surface.png'), dpi=150)
    plt.close(fig)

    report('Free surface: sources and receivers at and near the surface', passed,
           f'max amplitude error {worst["surface"]:.3f} at the surface (tol {TOL_SURFACE}), '
           f'{worst["deep"]:.3f} at 150 m (tol {TOL_DEEP}), max residual {worst["shape"]:.3f} '
           f'(tol {TOL_SHAPE}), acoustic image error {max(ac_err.values()):.1e} (tol {TOL_ACOUSTIC})')
    return passed, worst, ac_err


if __name__ == '__main__':
    ok = run()[0]
    sys.exit(0 if ok else 1)
