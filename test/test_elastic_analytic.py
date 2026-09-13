"""
Test 8: Far-field analytic wholespace response, 2D ELASTIC medium.

Elastic counterpart of Test 1.  Compares OWL's 2D isotropic elastic FD solution
(which_medium = elastic-iso) against the exact 2D elastodynamic Green's tensor
for a point (line) force in a homogeneous, unbounded medium -- no free surface,
PML on all four sides, so any discrepancy is attributable to the solver itself.

2D elastodynamic Green's tensor (plane strain, e^{+iwt} convention).  With the
scalar Helmholtz Green's functions

  G^P = -(i/4) H0^(2)(w r / vp),   G^S = -(i/4) H0^(2)(w r / vs),

the displacement at x_i due to a unit point force along j is

  G_ij(r, w) = 1/(rho w^2) * [ d_i d_j (G^S - G^P) + (w/vs)^2 delta_ij G^S ]

with  d_i d_j f(r) = f''(r) gamma_i gamma_j + (f'(r)/r)(delta_ij - gamma_i gamma_j)
and gamma_i = x_i / r.  In the far field (k r >> 1) this reduces to the familiar
P/S separation

  G_ij -> 1/rho [ gamma_i gamma_j G^P / vp^2 + (delta_ij - gamma_i gamma_j) G^S / vs^2 ],

i.e. P motion polarised along gamma with a cos radiation pattern and S motion
polarised perpendicular to gamma with a sin pattern.  The receivers below sit at
r = 600 m ~ 5 P-wavelengths / 9 S-wavelengths from the source, and the exact
tensor (not the far-field limit) is used for the comparison.

OWL records particle VELOCITY, so the analytic displacement is differentiated in
time (multiplied by i*w).  OWL's force source injects  dv = dt * amp * stf / rho
at one grid cell, which is the discrete form of a point force of
amp * stf * DX * DZ -- the same discrete-to-continuum factor as in Test 1.  The
force direction is (sin(polar), cos(polar)) in (x, z).

Geometry: a vertical point force (polar = 0) with receivers on a half circle of
radius r around the source, at azimuths phi measured from the +x axis.  This
sweeps the full P and S radiation patterns in one shot: at phi = 90 deg (directly
below the source) the motion is pure P on vz with a null S, while at phi = 0 deg
(broadside) it is pure S on vz with a null P.

Pass criterion: per-receiver relative L2 error of the VECTOR (vx, vz) trace
< 2%.
"""

import os, sys, shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hankel2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from owl_test_utils import (
    write_model, write_geometry, write_param, run_owl, read_su, ricker, report, Source
)

WORK = os.path.join(HERE, 'work', os.path.splitext(os.path.basename(__file__))[0])
PLOT = os.path.join(HERE, 'plots')
os.makedirs(WORK, exist_ok=True)
os.makedirs(PLOT, exist_ok=True)

# ── Model / acquisition ───────────────────────────────────────────────────────
NZ, NX = 401, 401
DZ, DX = 5.0, 5.0           # 5 m: ~23 points per S wavelength at f0
VP, VS = 3000.0, 1732.0     # Poisson solid (vp/vs = sqrt 3)
RHO    = 2000.0
F0     = 15.0
DT     = 4.0e-4
TMAX   = 0.55               # < 0.67 s, the earliest possible PML round trip
PML    = 30
AMP    = 1.0e6

SX, SZ = NX // 2 * DX, NZ // 2 * DZ   # source at the grid centre
POLAR  = 0.0                          # vertical force: (sin 0, cos 0) = +z
RADIUS = 600.0                        # receiver circle radius (m)
PHIS   = np.array([0.0, 30.0, 60.0, 90.0, 120.0, 150.0])   # azimuth from +x

# Staggered-grid source-injection / recording offset.  The force increment is
# added to v at the start of a step and v is recorded after it, so the numerical
# trace lags the continuum solution by exactly half a time step.  This is a
# property of the velocity-stress scheme, not a fitted parameter: an empirical
# scan of the shift puts the error minimum at 0.4 dt, i.e. within 0.1 dt of the
# predicted value, the small residue being FD dispersion.
TIME_SHIFT = 0.5 * DT

# Tighter than Test 1's 5%: all receivers sit at one radius, so there is no
# accumulating-dispersion spread across the spread, and the observed error is
# ~0.6%.  2% still leaves a 3x margin.
PASS_TOL = 0.02


def recv_positions():
    ph = np.radians(PHIS)
    return [(SX + RADIUS * np.cos(p), SZ + RADIUS * np.sin(p)) for p in ph]


def build_inputs():
    for nm, val in [('vp', VP), ('vs', VS), ('rho', RHO)]:
        write_model(os.path.join(WORK, 'model', f'{nm}.bin'), np.full((NZ, NX), val))
    write_geometry(
        os.path.join(WORK, 'geometry'),
        [Source(x=SX, z=SZ, f0=F0, amp=AMP, mechanism=f'force {POLAR:.6f} 0.0')],
        recv_positions())
    write_param(os.path.join(WORK, 'param.rb'), {
        'nx': NX, 'nz': NZ, 'dx': DX, 'dz': DZ,
        'npml': PML,
        'dt': DT, 'data_dt': DT, 'tmax': TMAX,
        'ns': 1,
        'file_geometry': './geometry/geometry.txt',
        'which_medium': 'elastic-iso',
        'model_name': 'vp, vs, rho',
        'file_vp':  './model/vp.bin',
        'file_vs':  './model/vs.bin',
        'file_rho': './model/rho.bin',
        'yn_free_surface': 'n',
        'dir_synthetic': './syn',
        'verbose': 'n',
    })


def owl_source_spectrum(npad, dt, f0, amp):
    """Spectrum of the equivalent point force  F(t) = amp * DX * DZ * stf(t).

    OWL's Ricker STF has its peak at t0 = 1/f0, a length of 2/f0, and is
    normalised to unit maximum.
    """
    nstf = int(round(2.0 / f0 / dt))
    stf = np.zeros(npad)
    for i in range(min(nstf + 1, npad)):
        stf[i] = ricker(i * dt, f0, 1.0 / f0)
    stf /= np.max(np.abs(stf))
    return np.fft.rfft(stf) * dt * amp * DX * DZ


def analytic_velocity(dxr, dzr, nt, dt, f0=F0, vp=VP, vs=VS, rho=RHO,
                      polar=POLAR, amp=AMP, time_shift=TIME_SHIFT):
    """Exact 2D elastic particle velocity (vx, vz) at offset (dxr, dzr).

    Uses the same Fourier conventions as Test 1: numpy's rfft is the forward
    transform, so the reconstruction carries e^{+iwt}, outgoing waves are
    H0^(2), and a time derivative is a factor +i*w.
    """
    npad = 8 * nt
    F = owl_source_spectrum(npad, dt, f0, amp)
    w = 2.0 * np.pi * np.fft.rfftfreq(npad, dt)

    r = float(np.hypot(dxr, dzr))
    g = np.array([dxr / r, dzr / r])              # gamma
    n = np.array([np.sin(polar), np.cos(polar)])  # force direction
    gn = float(g @ n)

    ka, kb = w / vp, w / vs
    c = -0.25j
    with np.errstate(divide='ignore', invalid='ignore'):
        HP0, HP1 = hankel2(0, ka * r), hankel2(1, ka * r)
        HS0, HS1 = hankel2(0, kb * r), hankel2(1, kb * r)
        # d_i d_j f = A gamma_i gamma_j + B (delta_ij - gamma_i gamma_j)
        AP, BP = -c * ka ** 2 * (HP0 - HP1 / (ka * r)), -c * ka * HP1 / r
        AS, BS = -c * kb ** 2 * (HS0 - HS1 / (kb * r)), -c * kb * HS1 / r
        GS = c * HS0
        # G_ij n_j = [ d_i d_j (G^S - G^P) + kb^2 delta_ij G^S ] n_j / (rho w^2)
        radial = (AS - AP) * gn                        # along gamma
        trans = (BS - BP)                              # along (n - gn*gamma)
        pref = 1.0 / (rho * w ** 2)
        U = [pref * (radial * g[i] + trans * (n[i] - gn * g[i])
                     + kb ** 2 * GS * n[i]) for i in range(2)]

    phase = np.exp(-1j * w * time_shift)
    out = []
    for i in range(2):
        V = (1j * w) * np.nan_to_num(U[i], nan=0.0, posinf=0.0, neginf=0.0) * F * phase
        V[0] = 0.0
        out.append(np.fft.irfft(V, n=npad)[:nt] / dt)
    return out[0], out[1]


def run():
    build_inputs()
    shutil.rmtree(os.path.join(WORK, 'syn'), ignore_errors=True)
    run_owl('owl_modeling2', 'param.rb', WORK)

    vx, _ = read_su(os.path.join(WORK, 'syn', 'shot_1_seismogram_x.su'))
    vz, _ = read_su(os.path.join(WORK, 'syn', 'shot_1_seismogram_z.su'))
    vx = vx.astype(np.float64)
    vz = vz.astype(np.float64)
    nt = vx.shape[0]

    recs = recv_positions()
    errs, ana = [], []
    for i, (rx, rz) in enumerate(recs):
        ax_, az_ = analytic_velocity(rx - SX, rz - SZ, nt, DT)
        ana.append((ax_, az_))
        num = np.concatenate([vx[:, i], vz[:, i]])
        ref = np.concatenate([ax_, az_])
        errs.append(np.linalg.norm(num - ref) / np.linalg.norm(ref))

    max_err = max(errs)
    passed = max_err < PASS_TOL
    report('Test 8 – Elastic analytic wholespace response', passed,
           f'max vector rel-L2 err = {max_err:.2%}  (tol {PASS_TOL:.0%})')
    for phi, e in zip(PHIS, errs):
        print(f'    phi = {phi:5.1f} deg :  rel err = {e:.2%}')

    _plot(vx, vz, ana, errs, nt, passed)
    return passed, errs


def _pick(trace_x, trace_z, t, t_arrival):
    """Peak vector amplitude of the phase arriving at t_arrival.

    OWL's Ricker STF peaks at t0 = 1/f0 after onset and has a total length of
    2/f0, so the wavelet occupies [t_arrival, t_arrival + 2/f0].  The P and S
    windows built this way are disjoint for the radius used here.
    """
    m = (t >= t_arrival) & (t <= t_arrival + 2.0 / F0)
    return np.max(np.hypot(trace_x[m], trace_z[m])) if m.any() else 0.0


def _plot(vx, vz, ana, errs, nt, passed):
    t = np.arange(nt) * DT
    tp, ts = RADIUS / VP, RADIUS / VS
    n = len(PHIS)

    fig, axes = plt.subplots(2, n, figsize=(3.0 * n, 6.0), sharex=True)
    for j, phi in enumerate(PHIS):
        # one amplitude scale per azimuth, shared by both components, so that a
        # radiation-pattern null plots as a flat line rather than as autoscaled
        # round-off noise
        col_max = max(np.abs(vx[:, j]).max(), np.abs(vz[:, j]).max(),
                      np.abs(ana[j][0]).max(), np.abs(ana[j][1]).max())
        for row, (num, ref, lab) in enumerate(
                [(vx[:, j], ana[j][0], 'v$_x$'), (vz[:, j], ana[j][1], 'v$_z$')]):
            ax = axes[row, j]
            ax.plot(t, num, 'b-', lw=1.3, label='OWL FD' if (row + j) == 0 else None)
            ax.plot(t, ref, 'r--', lw=1.2, label='Analytic' if (row + j) == 0 else None)
            ax.set_ylim(-1.15 * col_max, 1.15 * col_max)
            for tt, cc, nm in [(tp, 'g', 'P'), (ts, 'm', 'S')]:
                ax.axvline(tt, color=cc, ls=':', lw=0.9)
                if row == 0 and j == 0:
                    ax.text(tt, 0.97 * 1.15 * col_max, nm, color=cc, fontsize=8,
                            ha='center', va='top', fontweight='bold')
            # flag the exact nulls of the radiation pattern
            leak = np.abs(num).max() / col_max
            if np.abs(ref).max() / col_max < 1e-3:
                ax.text(0.5, 0.5, f'radiation-pattern null\nOWL leakage {leak:.0e}',
                        transform=ax.transAxes, ha='center', va='center',
                        fontsize=7, color='0.35')
            ax.set_xlim(0, TMAX)
            ax.tick_params(labelsize=7)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax.yaxis.get_offset_text().set_fontsize(6)
            if row == 0:
                ax.set_title(f'$\\phi$ = {phi:.0f}°  (err {errs[j]:.1%})', fontsize=9)
            else:
                ax.set_xlabel('Time (s)', fontsize=8)
            if j == 0:
                ax.set_ylabel(f'{lab}  (m/s)', fontsize=9)
            if (row + j) == 0:
                ax.legend(fontsize=7, loc='upper left')

    fig.suptitle('Test 8: OWL vs. 2D Analytic Elastic Wholespace Response  '
                 f'[{"PASSED" if passed else "FAILED"}]\n'
                 f'vp={VP:.0f} m/s, vs={VS:.0f} m/s, rho={RHO:.0f} kg/m³, '
                 f'f0={F0:.0f} Hz, dx=dz={DX:.0f} m, vertical force, r={RADIUS:.0f} m',
                 fontsize=10)
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(PLOT, 'test8_elastic_analytic.png'), dpi=150)
    plt.close()

    # ── error by azimuth + P/S radiation patterns ─────────────────────────────
    fig2 = plt.figure(figsize=(11, 4.2))
    ax1 = fig2.add_subplot(1, 2, 1)
    ax1.bar(PHIS, [e * 100 for e in errs], width=9.0,
            color=['green' if e < PASS_TOL else 'red' for e in errs])
    ax1.axhline(PASS_TOL * 100, color='red', ls='--', label=f'{PASS_TOL:.0%} threshold')
    ax1.set_xlabel('Receiver azimuth $\\phi$ (deg)')
    ax1.set_ylabel('Vector relative L2 error (%)')
    ax1.set_title('Test 8: waveform error by azimuth')
    ax1.set_xticks(PHIS)
    ax1.legend(fontsize=8)

    tp, ts = RADIUS / VP, RADIUS / VS
    phi_f = np.linspace(0, 180, 181)
    ax2 = fig2.add_subplot(1, 2, 2, projection='polar')
    num_p = np.array([_pick(vx[:, j], vz[:, j], t, tp) for j in range(n)])
    num_s = np.array([_pick(vx[:, j], vz[:, j], t, ts) for j in range(n)])
    ana_p = np.array([_pick(ana[j][0], ana[j][1], t, tp) for j in range(n)])
    ana_s = np.array([_pick(ana[j][0], ana[j][1], t, ts) for j in range(n)])
    # P and S carry different 1/vp^2, 1/vs^2 and geometrical-spreading factors,
    # so each is normalised by its own maximum; what is compared here is the
    # angular SHAPE of the two patterns.  For a force along +z the P amplitude
    # goes as |gamma.n| = |sin phi| and the S amplitude as |n - gamma(gamma.n)|
    # = |cos phi|, and both maxima are sampled by PHIS.
    ax2.plot(np.radians(phi_f), np.abs(np.sin(np.radians(phi_f))), 'g-', lw=1,
             alpha=0.4, label='P  $|\\sin\\phi|$')
    ax2.plot(np.radians(phi_f), np.abs(np.cos(np.radians(phi_f))), 'm-', lw=1,
             alpha=0.4, label='S  $|\\cos\\phi|$')
    ax2.plot(np.radians(PHIS), ana_p / ana_p.max(), 'g^', ms=9, mfc='none', label='P analytic')
    ax2.plot(np.radians(PHIS), num_p / ana_p.max(), 'gx', ms=7, label='P OWL')
    ax2.plot(np.radians(PHIS), ana_s / ana_s.max(), 'mo', ms=9, mfc='none', label='S analytic')
    ax2.plot(np.radians(PHIS), num_s / ana_s.max(), 'm+', ms=9, label='S OWL')
    # Two distinct quantities: how well OWL reproduces the EXACT pattern (the
    # actual test), and how far the exact pattern itself sits from the pure
    # far-field lobes (a physical O(1/kr) near-field residue, not an FD error).
    owl_dev = max(np.abs(num_p - ana_p).max() / ana_p.max(),
                  np.abs(num_s - ana_s).max() / ana_s.max())
    ff_dev = max(np.abs(ana_p / ana_p.max() - np.abs(np.sin(np.radians(PHIS)))).max(),
                 np.abs(ana_s / ana_s.max() - np.abs(np.cos(np.radians(PHIS)))).max())
    print(f'    radiation pattern: OWL vs exact = {owl_dev:.2%}; '
          f'exact vs far-field |sin|/|cos| = {ff_dev:.2%} '
          f'(near-field residue, k_p*r = {2 * np.pi * F0 * RADIUS / VP:.1f})')
    ax2.set_thetamin(0)
    ax2.set_thetamax(180)
    ax2.set_title('P / S radiation pattern\n(peak vector amplitude, normalised)',
                  fontsize=10)
    ax2.legend(fontsize=7, loc='center left', bbox_to_anchor=(1.02, 0.5))
    plt.tight_layout()
    fig2.savefig(os.path.join(PLOT, 'test8_elastic_errors.png'), dpi=150,
                 bbox_inches='tight')
    plt.close()
    print('figures ->', os.path.join(PLOT, 'test8_elastic_analytic.png'))


if __name__ == '__main__':
    ok, _ = run()
    sys.exit(0 if ok else 1)
