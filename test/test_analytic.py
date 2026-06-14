"""
Test 1: Far-field analytic wholespace response.

Compares OWL's 2D acoustic FD solution against the analytic 2D acoustic
Green's function in a homogeneous full-space medium (no free surface, no
reflections from PML).

2D frequency-domain Green's function for constant density media:
  p̂(r, ω) = -(i/4) * H₀⁽²⁾(ω·r/c) * Ŝ(ω)

where H₀⁽²⁾ is the Hankel function of the second kind, r is source–receiver
distance, and Ŝ(ω) is the Fourier transform of the source time function.

Pass criterion: relative L2 error per trace < 5% (FD spatial dispersion at
15 Hz, 10 m grid: ~4 points per wavelength at 50 Hz is tight, but the
dominant energy is at lower frequencies where dispersion is smaller).
"""

import os, sys, shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hankel2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from owl_test_utils import (
    write_model, read_model, read_su, write_geometry, write_param,
    run_owl, ricker, rel_err, report, Source
)

WORK = os.path.join(HERE, 'work', os.path.splitext(os.path.basename(__file__))[0])
PLOT = os.path.join(HERE, 'plots')
os.makedirs(WORK, exist_ok=True)
os.makedirs(PLOT, exist_ok=True)

# ── Model parameters ──────────────────────────────────────────────────────────
NZ, NX = 201, 201
DZ, DX = 10.0, 10.0
VP   = 2000.0
RHO  = 1.0          # unit density: simplifies Green's function formula
F0   = 15.0
DT   = 1.0e-3
TMAX = 2.0

SX, SZ = NX // 2 * DX, NZ // 2 * DZ   # source at center

# Receivers at horizontal offsets from 100 to 800 m (same depth as source)
OFFSETS = np.arange(100, 900, 100, dtype=float)   # 8 receivers

def build_model():
    vp  = np.full((NZ, NX), VP)
    rho = np.full((NZ, NX), RHO)
    write_model(os.path.join(WORK, 'model', 'vp.bin'),  vp)
    write_model(os.path.join(WORK, 'model', 'rho.bin'), rho)

def build_geometry():
    recvs = [(SX + r, SZ) for r in OFFSETS]
    write_geometry(
        os.path.join(WORK, 'geometry'),
        [Source(x=SX, z=SZ, f0=F0, amp=1.0)],
        recvs
    )

def build_params():
    write_param(os.path.join(WORK, 'param.rb'), {
        'nx': NX, 'nz': NZ, 'dx': DX, 'dz': DZ,
        'dt': DT, 'data_dt': DT, 'tmax': TMAX,
        'ns': 1,
        'file_geometry': './geometry/geometry.txt',
        'which_medium': 'acoustic-iso',
        'model_name': 'vp, rho',
        'file_vp':  './model/vp.bin',
        'file_rho': './model/rho.bin',
        'verbose': 'n',
    })

def analytic_trace(offset, nt, dt, f0, vp, amp=1.0, time_shift=1.5*DT):
    """2D analytic pressure for a Ricker explosion source at distance r.

    Accounts for a 1.5-sample time shift arising from the staggered-grid
    FD scheme's source injection ordering.
    """
    npad = 8 * nt
    t    = np.arange(npad) * dt
    # OWL Ricker STF: peak at 1/f0, normalized to max 1, length 2/f0
    nstf = int(round(2.0 / f0 / dt))
    stf  = np.zeros(npad)
    for i in range(min(nstf + 1, npad)):
        stf[i] = ricker(i * dt, f0, 1.0 / f0)
    stf /= np.max(np.abs(stf))   # OWL normalizes STF to max 1

    S = np.fft.rfft(stf) * dt
    f = np.fft.rfftfreq(npad, dt)
    w = 2.0 * np.pi * f

    r = float(offset)
    with np.errstate(divide='ignore', invalid='ignore'):
        # G is the 2D acoustic pressure Green's function including ∂S/∂t source factor
        # and the DX*DZ discrete-to-continuum normalisation.
        # Derivation: ∂²p/∂t² - c²∇²p = ∂S/∂t·δ(r)  →
        #   p̂ = (−iω·Ŝ)·(i/4)·H₀⁽²⁾(kr)/c²  [e⁻ⁱωᵗ numpy FFT; H₀⁽²⁾ matches OWL]
        G = (-1j / (4.0 * vp ** 2)) * hankel2(0, w * r / vp)
    G[0] = 0.0

    # Phase shift for FD staggered-grid source injection ordering (empirical 1.5 dt)
    phase = np.exp(-1j * w * time_shift)

    # Full pressure spectrum: p̂ = DX·DZ · amp · (−iω) · Ŝ · G · phase
    P = DX * DZ * amp * (1j * w) * S * G * phase
    p = np.fft.irfft(P, n=npad)[:nt] / dt
    return p

def run():
    build_model()
    build_geometry()
    build_params()
    shutil.rmtree(os.path.join(WORK, 'data_synthetic'), ignore_errors=True)
    run_owl('owl_modeling2', 'param.rb', WORK)

    d, _ = read_su(os.path.join(WORK, 'data_synthetic', 'shot_1_seismogram_p.su'))
    nt   = d.shape[0]

    errs = []
    traces_num = []
    traces_ana = []
    for itr, r in enumerate(OFFSETS):
        p_num = d[:, itr].astype(np.float64)
        p_ana = analytic_trace(r, nt, DT, F0, VP)
        e     = rel_err(p_num, p_ana)
        errs.append(e)
        traces_num.append(p_num)
        traces_ana.append(p_ana)

    max_err = max(errs)
    passed  = max_err < 0.05     # < 5 % relative L2 error per trace
    report('Test 1 – Analytic wholespace response',
           passed, f'max rel-L2 err = {max_err:.2%}')

    # ── Visualization ──────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 4, figsize=(14, 6))
    t_ax = np.arange(nt) * DT
    for idx, (ax, r, pn, pa, e) in enumerate(
            zip(axes.ravel(), OFFSETS, traces_num, traces_ana, errs)):
        ax.plot(t_ax, pn, 'b-',  lw=1.2, label='OWL FD')
        ax.plot(t_ax, pa, 'r--', lw=1.2, label='Analytic')
        ax.set_xlim(0, TMAX)
        ax.set_title(f'r = {r:.0f} m  (err {e:.1%})', fontsize=9)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Pressure')
        ax.tick_params(labelsize=8)
        if idx == 0:
            ax.legend(fontsize=8)

    fig.suptitle('Test 1: OWL vs. 2D Analytic Wholespace Response\n'
                 f'vp={VP} m/s, rho={RHO} kg/m³, f0={F0} Hz, '
                 f'dx=dz={DX} m   [{"PASSED" if passed else "FAILED"}]',
                 fontsize=10)
    plt.tight_layout()
    fig.savefig(os.path.join(PLOT, 'test1_analytic.png'), dpi=150)
    plt.close()

    # Error-vs-offset bar chart
    fig2, ax2 = plt.subplots(figsize=(7, 4))
    ax2.bar(OFFSETS, [e * 100 for e in errs], color=['green' if e < 0.05 else 'red' for e in errs])
    ax2.axhline(5.0, color='red', ls='--', label='5 % threshold')
    ax2.set_xlabel('Offset (m)')
    ax2.set_ylabel('Relative L2 error (%)')
    ax2.set_title('Test 1: Waveform error by offset')
    ax2.legend()
    fig2.savefig(os.path.join(PLOT, 'test1_errors.png'), dpi=150)
    plt.close()

    return passed, errs

if __name__ == '__main__':
    passed, errs = run()
    sys.exit(0 if passed else 1)
