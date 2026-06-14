"""
Native reference solution for the 2D Lamb problem (line source on the free
surface of a homogeneous, isotropic elastic half-space), computed by
frequency–wavenumber (discrete-wavenumber) integration.

Geometry (surface-aligned frame): half-space occupies z >= 0 with a free
surface at z = 0.  A point (line) force is applied at the origin on the
surface; receivers sit on the surface at along-surface offsets r.  We return
PARTICLE VELOCITY (to match owl_modeling2, which records velocity), for a
force directed along the surface NORMAL (the +z axis of this frame).

Plane-strain P-SV.  For a surface vertical (normal) point force the surface
displacement Green's functions in the (k, omega) domain are

    u_z(k) = (omega^2/beta^2) * nu_p / (mu * R)            (even in k)
    u_x(k) =  i k (Q - 2 nu_p nu_s) / (mu * R)             (odd  in k)

with  nu_p = sqrt(k^2 - (omega/alpha)^2),  nu_s = sqrt(k^2 - (omega/beta)^2),
Q = 2 k^2 - (omega/beta)^2,  R = Q^2 - 4 k^2 nu_p nu_s  (Rayleigh function),
mu = rho beta^2.  The Rayleigh pole (R = 0) is moved off the real-k axis by
giving omega a small negative imaginary part (Bouchon's method): omega ->
omega - i*omega_I, integrate over real k, then multiply the time series by
exp(+omega_I t).

The absolute scaling (1/2pi, source constants, the body-force-vs-traction
equivalence) is irrelevant for waveform validation — callers fit a single
global scalar — but is kept for cleanliness.
"""

import numpy as np


def ricker_owl(nt, dt, f0):
    """OWL's Ricker source time function: peak at t0 = 1/f0, length 2/f0,
    normalised to unit maximum (matches the wavelet owl injects)."""
    t0 = 1.0 / f0
    nstf = int(round(2.0 / f0 / dt))
    s = np.zeros(nt)
    for i in range(min(nstf + 1, nt)):
        a = (np.pi * f0 * (i * dt - t0)) ** 2
        s[i] = (1.0 - 2.0 * a) * np.exp(-a)
    return s / np.max(np.abs(s))


def lamb_halfspace_velocity(offsets, nt, dt, vp, vs, rho, f0,
                            nk=6000, kmax_fac=6.0, damp=4.0):
    """2D Lamb half-space surface particle-velocity seismograms.

    Parameters
    ----------
    offsets : 1d array of along-surface source-receiver distances (m), > 0.
    nt, dt  : time samples and step (s) — match the owl recording.
    vp,vs,rho : half-space properties.
    f0      : Ricker centre frequency (same wavelet as owl).
    nk      : number of wavenumbers in the integration.
    kmax_fac: k integrated over [-kmax, kmax], kmax = kmax_fac * omega_max / vs.
    damp    : imaginary-frequency factor; omega_I = damp / (nt*dt).

    Returns
    -------
    vz, vx : arrays (nt, n_offset) — surface-normal (z) and surface-tangential
             (x) particle velocity.  vz is the symmetric component, vx the
             antisymmetric one (sign follows +r direction).
    """
    offsets = np.atleast_1d(np.asarray(offsets, dtype=np.float64))
    alpha, beta = float(vp), float(vs)
    mu = rho * beta ** 2

    T = nt * dt
    omega_I = damp / T

    # frequency grid (real FFT)
    freqs = np.fft.rfftfreq(nt, dt)
    omega_r = 2.0 * np.pi * freqs
    omega_max = omega_r[-1]

    # source spectrum (same Ricker as owl)
    stf = ricker_owl(nt, dt, f0)
    STF = np.fft.rfft(stf)

    # wavenumber grid
    kmax = kmax_fac * omega_max / beta
    k = np.linspace(-kmax, kmax, nk)
    dk = k[1] - k[0]
    expo = np.exp(1j * np.outer(k, offsets))        # (nk, n_off)

    Vz = np.zeros((len(freqs), len(offsets)), dtype=np.complex128)
    Vx = np.zeros_like(Vz)

    for iw in range(1, len(freqs)):                 # skip DC (Ricker DC ~ 0)
        w = omega_r[iw] - 1j * omega_I
        nup = np.sqrt(k ** 2 - (w / alpha) ** 2)    # principal sqrt -> Re>=0
        nus = np.sqrt(k ** 2 - (w / beta) ** 2)
        Q = 2.0 * k ** 2 - (w / beta) ** 2
        R = Q ** 2 - 4.0 * k ** 2 * nup * nus

        uz_hat = (w ** 2 / beta ** 2) * nup / (mu * R)        # even in k
        ux_hat = 1j * k * (Q - 2.0 * nup * nus) / (mu * R)    # odd  in k

        # surface displacement Green's functions G(r) = (1/2pi) ∫ u_hat e^{ikr} dk
        Gz = (uz_hat @ expo) * dk / (2.0 * np.pi)
        Gx = (ux_hat @ expo) * dk / (2.0 * np.pi)

        # velocity = iω · STF(ω) · G(r)
        Vz[iw] = 1j * w * STF[iw] * Gz
        Vx[iw] = 1j * w * STF[iw] * Gx

    # back to time, undo the imaginary-frequency damping
    vz = np.fft.irfft(Vz, n=nt, axis=0)
    vx = np.fft.irfft(Vx, n=nt, axis=0)
    grow = np.exp(omega_I * np.arange(nt) * dt)[:, None]
    return (vz * grow).real, (vx * grow).real


def rayleigh_speed(vp, vs):
    """Rayleigh wave speed (real root of the Rayleigh equation)."""
    # solve (2 - x)^2 = 4 sqrt(1-x) sqrt(1 - x vs^2/vp^2), x = (vR/vs)^2
    g = (vs / vp) ** 2
    xs = np.linspace(1e-6, 0.9999, 200000)
    f = (2.0 - xs) ** 2 - 4.0 * np.sqrt(1.0 - xs) * np.sqrt(1.0 - g * xs)
    i = np.where(np.diff(np.sign(f)))[0][-1]
    x = xs[i]
    return vs * np.sqrt(x)


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    vp, vs, rho = 3000.0, 1732.0, 2000.0     # Poisson solid (vp/vs = sqrt3)
    f0 = 12.5
    dt = 4.0e-4
    nt = 2500
    T = nt * dt
    offsets = np.arange(200.0, 1001.0, 100.0)

    vz, vx = lamb_halfspace_velocity(offsets, nt, dt, vp, vs, rho, f0)
    t = np.arange(nt) * dt

    vR = rayleigh_speed(vp, vs)
    print(f'vp={vp} vs={vs}  Rayleigh speed vR = {vR:.1f} m/s  (vR/vs={vR/vs:.4f})')
    print(f'expected Rayleigh arrivals (s): '
          + ', '.join(f'{r/vR:.3f}' for r in offsets))
    print(f'expected P arrivals       (s): '
          + ', '.join(f'{r/vp:.3f}' for r in offsets))

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    for comp, ax, V, title in [('z', axes[0], vz, 'vertical (normal) v_z'),
                               ('x', axes[1], vx, 'horizontal (tang.) v_x')]:
        sc = 70.0 / np.max(np.abs(V))
        for j, r in enumerate(offsets):
            ax.plot(t, r + sc * V[:, j], 'k', lw=0.8)
        ax.plot(offsets / vp, offsets, 'b--', lw=1, label='P  (r/vp)')
        ax.plot(offsets / vs, offsets, 'g--', lw=1, label='S  (r/vs)')
        ax.plot(offsets / vR, offsets, 'r--', lw=1, label='Rayleigh (r/vR)')
        ax.set_title(f'Lamb half-space — {title}')
        ax.set_xlabel('time (s)'); ax.set_xlim(0, T)
        ax.invert_yaxis()
        if comp == 'z':
            ax.set_ylabel('offset (m)'); ax.legend(fontsize=8)
    plt.tight_layout()
    out = __file__.replace('lamb_reference.py', 'plots/lamb_reference_selftest.png')
    fig.savefig(out, dpi=130)
    print('self-test figure ->', out)
