"""
Shared utilities for OWL correctness unit tests.

These helpers wrap the OWL executables (owl_modeling2, owl_modeling3, owl_fwi2)
and handle the binary model format, the ASCII geometry format, and the SU
seismogram format used by OWL.

All tests use the 2D/3D isotropic acoustic solver unless noted otherwise.
"""

import os
import subprocess
import numpy as np

# ------------------------------------------------------------------------------
# Binary model I/O (OWL models are float32, Fortran order, z fastest)
# ------------------------------------------------------------------------------


def write_model(path, w):
    """Write a model array to OWL binary format.

    2D: w has shape (nz, nx); 3D: w has shape (nz, ny, nx).
    The file is float32 with z the fastest dimension, which corresponds to
    Fortran-order storage of w(nz, nx) -- i.e., C-order of w[x, z] reversed.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.asarray(w, dtype=np.float32).T.astype(np.float32).tofile(path)


def read_model(path, shape):
    """Read an OWL binary model/gradient file. shape=(nz, nx) or (nz, ny, nx)."""
    w = np.fromfile(path, dtype=np.float32)
    expected = int(np.prod(shape))
    if w.size != expected:
        raise ValueError(f"{path}: expected {expected} floats, got {w.size}")
    return w.reshape(tuple(reversed(shape))).T


# ------------------------------------------------------------------------------
# SU seismogram I/O
# ------------------------------------------------------------------------------


def read_su(path):
    """Read a SU file written by OWL (native little-endian, 240-byte headers).

    Returns (data, dt) where data has shape (nt, ntrace).
    """
    raw = np.fromfile(path, dtype=np.uint8)
    if raw.size < 240:
        raise ValueError(f"{path}: too small to be a SU file")
    ns = int(np.frombuffer(raw[114:116].tobytes(), dtype=np.int16)[0])
    dt_us = int(np.frombuffer(raw[116:118].tobytes(), dtype=np.uint16)[0])
    trace_bytes = 240 + 4 * ns
    if raw.size % trace_bytes != 0:
        raise ValueError(f"{path}: size {raw.size} not divisible by trace size {trace_bytes}")
    ntr = raw.size // trace_bytes
    raw = raw.reshape(ntr, trace_bytes)
    data = raw[:, 240:].copy().view(np.float32).reshape(ntr, ns).T
    return data, dt_us * 1.0e-6


def write_su(path, data, dt, like=None):
    """Write (nt, ntrace) array to a SU file.

    If `like` is given, headers are copied from that file (so OWL sees the
    same geometry headers); otherwise minimal headers (ns, dt) are created.
    """
    data = np.asarray(data, dtype=np.float32)
    nt, ntr = data.shape
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if like is not None:
        raw = np.fromfile(like, dtype=np.uint8)
        ns = int(np.frombuffer(raw[114:116].tobytes(), dtype=np.int16)[0])
        if ns != nt:
            raise ValueError(f"write_su: ns mismatch with template ({ns} vs {nt})")
        trace_bytes = 240 + 4 * ns
        raw = raw.reshape(-1, trace_bytes).copy()
        if raw.shape[0] != ntr:
            raise ValueError("write_su: trace count mismatch with template")
        raw[:, 240:] = data.T.copy().view(np.uint8).reshape(ntr, 4 * nt)
        raw.tofile(path)
    else:
        with open(path, "wb") as f:
            for i in range(ntr):
                hdr = np.zeros(240, dtype=np.uint8)
                hdr[114:116] = np.frombuffer(np.int16(nt).tobytes(), dtype=np.uint8)
                hdr[116:118] = np.frombuffer(np.uint16(round(dt * 1e6)).tobytes(), dtype=np.uint8)
                f.write(hdr.tobytes())
                f.write(data[:, i].tobytes())


# ------------------------------------------------------------------------------
# Geometry files
# ------------------------------------------------------------------------------


class Source:
    def __init__(self, x, z, y=0.0, mechanism="explosion", wavelet="ricker",
                 f0=15.0, amp=1.0, t0=0.0, stf_file=None):
        self.x, self.y, self.z = x, y, z
        self.mechanism = mechanism
        self.wavelet = wavelet
        self.f0, self.amp, self.t0 = f0, amp, t0
        self.stf_file = stf_file


def write_geometry(geom_dir, shots, receivers):
    """Write OWL geometry files.

    shots: list of (list of Source) -- one inner list per shot.
    receivers: list of (x, z) or (x, y, z) tuples (same receivers every shot),
               or a list of such lists (one per shot).
    """
    os.makedirs(geom_dir, exist_ok=True)
    if not isinstance(shots[0], (list, tuple)):
        shots = [[s] for s in shots]
    per_shot_rec = isinstance(receivers[0], (list,)) and isinstance(receivers[0][0], (tuple, list))
    with open(os.path.join(geom_dir, "geometry.txt"), "w") as f:
        for i in range(len(shots)):
            f.write(f"shot_{i + 1}_geometry.txt\n")
    for i, srcs in enumerate(shots):
        recs = receivers[i] if per_shot_rec else receivers
        with open(os.path.join(geom_dir, f"shot_{i + 1}_geometry.txt"), "w") as f:
            f.write(f"{i + 1}\n")
            f.write(f"{len(srcs)}\n")
            for s in srcs:
                f.write(f"{s.x} {s.y} {s.z}\n")
                f.write(f"{s.mechanism}\n")
                f.write(f"{s.wavelet} {s.f0} {s.amp} {s.t0}\n")
                if s.wavelet == "custom":
                    assert s.stf_file is not None
                    f.write(f"{s.stf_file}\n")
                f.write("0 0\n")
            f.write(f"{len(recs)}\n")
            for r in recs:
                if len(r) == 2:
                    f.write(f"{r[0]} 0.0 {r[1]} 1.0\n")
                else:
                    f.write(f"{r[0]} {r[1]} {r[2]} 1.0\n")


def write_stf(path, t, amp):
    """Write a custom source time function as a two-column ASCII file."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for ti, ai in zip(t, amp):
            f.write(f"{ti:.10e} {ai:.10e}\n")


# ------------------------------------------------------------------------------
# Parameter files and runners
# ------------------------------------------------------------------------------


def write_param(path, params):
    if os.path.dirname(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for k, v in params.items():
            f.write(f"{k} = {v}\n")


def run_owl(exe, param_file, cwd, np_ranks=1, quiet=True):
    """Run an OWL executable through mpirun and raise on failure."""
    cmd = ["mpirun", "-np", str(np_ranks), exe, param_file]
    res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(
            f"{' '.join(cmd)} failed (exit {res.returncode}) in {cwd}\n"
            f"--- stdout ---\n{res.stdout[-4000:]}\n--- stderr ---\n{res.stderr[-4000:]}"
        )
    if not quiet:
        print(res.stdout)
    return res


def ricker(t, f0, t0):
    """Ricker wavelet with peak at t0 (matches OWL's analytic 'ricker')."""
    a = (np.pi * f0 * (t - t0)) ** 2
    return (1 - 2 * a) * np.exp(-a)


def rel_err(a, b):
    """Relative L2 error ||a - b|| / ||b||."""
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def report(name, passed, detail=""):
    status = "PASSED" if passed else "FAILED"
    print(f"[{status}] {name}" + (f" -- {detail}" if detail else ""))
    return passed
