"""
Reference solution for a line source at any depth in a homogeneous, isotropic elastic
half-space with a flat free surface (or in a whole space), computed by
frequency-wavenumber (discrete-wavenumber) integration.

This generalizes lamb_reference.py (a vertical force on the surface, surface receivers) to

    * a line force (fx, fz) or a line moment tensor (mxx, mzz, mxz), e.g. an explosion
      mxx = mzz = 1, at depth h >= 0,
    * receivers at any depth z >= 0,
    * absolute amplitudes: the source and the receivers follow owl_modeling2's units, so
      that OWL and the reference can be compared without fitting a scale.

Geometry: plane strain (P-SV), x horizontal, z downward, half-space z >= 0 with a
traction-free surface at z = 0. Fields are expanded in exp(i k x - i omega t). In a
homogeneous layer the displacement-stress vector b = (u_x, u_z, s_xz, s_zz) is a sum of
four plane waves exp(i q z), q = +nu (down-going) or -nu (up-going) for P and S, with
nu = sqrt(omega^2/c^2 - k^2), Im(nu) > 0. The source makes b jump across z = h
(Aki & Richards; Kennett):

    force       [s_xz] = -fx,  [s_zz] = -fz
    moment      [u_x] = mxz/mu,  [u_z] = mzz/(lam + 2 mu),
                [s_xz] = i k (mxx - lam mzz/(lam + 2 mu)),  [s_zz] = 0

The jump fixes the four whole-space amplitudes; the free surface adds down-going P and S
waves that cancel the traction at z = 0. omega is moved into the upper half plane
(omega + i omega_I) so that the integrand is smooth on the real k axis (Bouchon's method);
the time series is multiplied by exp(omega_I t) at the end.

Units (owl_modeling2 conventions; amp = 1 in the geometry file):
    force       owl adds dt*stf(t)/rho to the particle velocity at the source node, i.e. a
                body force density stf(t) on one grid cell; the equivalent line force is
                F(t) = stf(t)*dx*dz (pass cell = dx*dz).
    explosion   owl adds -dt*stf'(t)/(dx*dz) to sxx and szz at the source node, i.e. a
                line moment M(t) = stf(t) (pass cell = 1).
"""

import numpy as np

from lamb_reference import ricker_owl


def _plane_waves(k, w, vp, vs, rho):
    """P and S eigenvectors of b for exp(i q z), q = +nu (down) and -nu (up)."""
    mu = rho*vs**2
    nup = np.sqrt((w/vp)**2 - k**2)
    nus = np.sqrt((w/vs)**2 - k**2)
    # principal square root; keep Im(nu) >= 0 (decay with depth for down-going waves)
    nup = np.where(nup.imag < 0, -nup, nup)
    nus = np.where(nus.imag < 0, -nus, nus)

    def ep(q):
        return np.stack([1j*k, 1j*q, -2.0*mu*k*q, mu*(2.0*k**2 - (w/vs)**2)], axis=-1)

    def es(q):
        return np.stack([-1j*q, 1j*k, mu*(q**2 - k**2), -2.0*mu*k*q], axis=-1)

    return nup, nus, ep(nup), es(nus), ep(-nup), es(-nus)


def halfspace_velocity(xr, zr, h, nt, dt, vp, vs, rho, f0, fx=0.0, fz=0.0, mxx=0.0, mzz=0.0, mxz=0.0,
                       cell=1.0, free_surface=True, stf=None, fmax=None, smooth=1.0, damp=4.0, pad=2,
                       period=None):
    """
    Particle velocity at receivers (xr, zr) from a line source at (0, h).

    xr, zr      receiver coordinates (m), arrays of equal length; zr >= 0, zr != h
    h           source depth (m), >= 0
    fx, fz      line force components (owl: polar angle of the force, sin and cos)
    mxx ... mxz moment tensor components
    cell        factor converting owl's injected density to the line source (see module doc)
    stf         source time function (nt); default: owl's Ricker of centre frequency f0
    fmax        highest frequency integrated (default: 6 f0)
    smooth      (m) the wavenumber integral is tapered by exp(-(k smooth)^2/4), i.e. the solution
                is smoothed along x by a Gaussian of this width; this regularizes sources and
                receivers on the surface, where the integrand does not decay with k
    damp        imaginary frequency omega_I = damp/T, T the (padded) time window
    pad         the solution is computed over pad*nt samples and truncated, which keeps the
                wrap-around of the discrete Fourier transform out of the returned window

    Returns vx, vz with shape (nt, nrec).
    """
    xr = np.atleast_1d(np.asarray(xr, dtype=np.float64))
    zr = np.atleast_1d(np.asarray(zr, dtype=np.float64))
    lam = rho*(vp**2 - 2.0*vs**2)
    mu = rho*vs**2

    nt0 = nt
    nt = pad*nt0
    T = nt*dt
    omega_I = damp/T
    freqs = np.fft.rfftfreq(nt, dt)
    fmax = 6.0*f0 if fmax is None else fmax
    nw = int(np.searchsorted(freqs, fmax)) + 1

    if stf is None:
        stf = ricker_owl(nt0, dt, f0)
    stf = np.concatenate([np.asarray(stf, dtype=np.float64), np.zeros(nt - nt0)])
    # Spectrum of the damped source time function in the exp(-i w t) convention of the plane
    # waves: S(w) = sum s(t) e^{i w t} dt with w = 2 pi f + i omega_I
    t = np.arange(nt)*dt
    stf_damped = stf*np.exp(-omega_I*t)

    # Wavenumber grid: the spacing 2 pi/period sets the spatial period of the implied source
    # array; the taper exp(-(k smooth)^2/4) is below 1e-4 at kmax
    if period is None:
        period = 2.0*(np.max(np.abs(xr)) + vp*T)
    dk = 2.0*np.pi/period
    kmax = 6.0/smooth
    nk = 2*int(np.ceil(kmax/dk)) + 1
    k = (np.arange(nk) - nk//2)*dk
    taper = np.exp(-(k*smooth)**2/4.0)
    expo = np.exp(1j*np.outer(k, xr))*(taper*dk/(2.0*np.pi))[:, None]   # (nk, nrec)

    Vx = np.zeros((len(freqs), len(xr)), dtype=np.complex128)
    Vz = np.zeros_like(Vx)

    for iw in range(1, nw):
        w = 2.0*np.pi*freqs[iw] + 1j*omega_I
        S = np.sum(stf_damped*np.exp(1j*2.0*np.pi*freqs[iw]*t))*dt*cell

        nup, nus, epd, esd, epu, esu = _plane_waves(k, w, vp, vs, rho)

        # Source jump of b = (u_x, u_z, s_xz, s_zz) across z = h
        s = np.zeros((nk, 4), dtype=np.complex128)
        s[:, 0] = mxz/mu
        s[:, 1] = mzz/(lam + 2.0*mu)
        s[:, 2] = -fx + 1j*k*(mxx - lam/(lam + 2.0*mu)*mzz)
        s[:, 3] = -fz

        # Whole-space amplitudes: below = a0 P_down + a1 S_down, above = a2 P_up + a3 S_up
        A = np.stack([epd, esd, -epu, -esu], axis=-1)     # (nk, 4, 4), columns = waves
        a = np.linalg.solve(A, s[..., None])[..., 0]

        # Free surface: down-going waves r0 P + r1 S that cancel the traction at z = 0
        if free_surface:
            up0 = a[:, 2:3]*epu*np.exp(1j*nup*h)[:, None] + a[:, 3:4]*esu*np.exp(1j*nus*h)[:, None]
            B = np.stack([epd[:, 2:4], esd[:, 2:4]], axis=-1)  # (nk, 2, 2)
            r = np.linalg.solve(B, -up0[:, 2:4, None])[..., 0]
        else:
            r = np.zeros((nk, 2), dtype=np.complex128)

        for j in range(len(xr)):
            z = zr[j]
            if z > h:
                b = a[:, 0:1]*epd*np.exp(1j*nup*(z - h))[:, None] + a[:, 1:2]*esd*np.exp(1j*nus*(z - h))[:, None]
            else:
                b = a[:, 2:3]*epu*np.exp(-1j*nup*(z - h))[:, None] + a[:, 3:4]*esu*np.exp(-1j*nus*(z - h))[:, None]
            b = b + r[:, 0:1]*epd*np.exp(1j*nup*z)[:, None] + r[:, 1:2]*esd*np.exp(1j*nus*z)[:, None]
            # displacement at x: (1/2 pi) sum_k b(k) e^{ikx} dk; velocity = -i w u
            ux = np.sum(b[:, 0]*expo[:, j])
            uz = np.sum(b[:, 1]*expo[:, j])
            Vx[iw, j] = -1j*w*S*ux
            Vz[iw, j] = -1j*w*S*uz

    # Back to time: the spectra are in the exp(-i w t) convention of the plane waves, so the
    # time series of the damped signal is the inverse real FFT of their complex conjugates
    vx = np.fft.irfft(np.conj(Vx), n=nt, axis=0)/dt
    vz = np.fft.irfft(np.conj(Vz), n=nt, axis=0)/dt
    grow = np.exp(omega_I*t)[:, None]
    return (vx*grow)[:nt0], (vz*grow)[:nt0]


def _check_plane_waves():
    """The eigenvectors satisfy d b/dz = A b, i.e. A e = i q e."""
    vp, vs, rho = 3000.0, 1732.0, 2000.0
    lam, mu = rho*(vp**2 - 2*vs**2), rho*vs**2
    k = np.array([0.003, 0.02, 0.5])
    w = 2*np.pi*12.0 + 0.3j
    nup, nus, epd, esd, epu, esu = _plane_waves(k, w, vp, vs, rho)
    for n in range(len(k)):
        kk = k[n]
        A = np.array([[0, -1j*kk, 1/mu, 0],
                      [-1j*kk*lam/(lam + 2*mu), 0, 0, 1/(lam + 2*mu)],
                      [-rho*w**2 + 4*kk**2*mu*(lam + mu)/(lam + 2*mu), 0, 0, -1j*kk*lam/(lam + 2*mu)],
                      [0, -rho*w**2, -1j*kk, 0]])
        for e, q in ((epd[n], nup[n]), (esd[n], nus[n]), (epu[n], -nup[n]), (esu[n], -nus[n])):
            assert np.allclose(A @ e, 1j*q*e, rtol=1e-10, atol=1e-10*np.abs(e).max())
    return True


if __name__ == '__main__':
    print('plane-wave check:', _check_plane_waves())
