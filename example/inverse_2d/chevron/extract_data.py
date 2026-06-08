#!/usr/bin/env python3

import os
import argparse
import struct
from collections import defaultdict

import numpy as np
import segyio


# -----------------------------
# VP MODEL EXPORT
# -----------------------------

def read_vp_model_segy(vp_segy):
    """
    Read a SEG-Y VP model.

    Assumption:
        Each SEG-Y trace is one vertical column.
        Number of traces = nx
        Samples per trace = nz

    Returns:
        vp with shape (nz, nx)
    """

    with segyio.open(vp_segy, "r", ignore_geometry=True) as f:
        f.mmap()

        nx = f.tracecount
        nz = len(f.samples)

        # segyio returns traces as shape (nx, nz)
        traces = segyio.tools.collect(f.trace[:])

        # Convert to vp[z, x]
        vp = traces.T.astype(np.float32)

        # Try to infer dz from SEG-Y sample axis
        samples = np.asarray(f.samples, dtype=np.float64)
        if len(samples) > 1:
            dz = float(np.median(np.diff(samples)))
        else:
            dz = 1.0

        # Try to infer dx from source x coordinates
        dx = infer_dx_from_headers(f)

    return vp, nx, nz, dx, dz


def infer_dx_from_headers(segyfile):
    """
    Try to infer dx from trace headers.
    Returns 1.0 if it cannot be inferred.
    """

    keys = [
        segyio.TraceField.SourceX,
        segyio.TraceField.GroupX,
        segyio.TraceField.CDP_X,
    ]

    for key in keys:
        try:
            x = np.asarray(segyfile.attributes(key)[:], dtype=np.float64)
        except Exception:
            continue

        ux = np.unique(x)

        if len(ux) > 1:
            dx_values = np.diff(np.sort(ux))
            dx_values = dx_values[dx_values != 0]

            if len(dx_values) > 0:
                return float(np.median(np.abs(dx_values)))

    print("Warning: could not infer dx from SEG-Y headers. Using dx = 1.0")
    return 1.0


def export_raw_binary(vp, output_raw):
    """
    Export VP model as little-endian float32 raw binary.

    The vp array has shape (nz, nx).

    Column-based means:
        write vp[:, 0], then vp[:, 1], then vp[:, 2], ...

    This is equivalent to Fortran-order output for an array shaped (nz, nx).
    """

    vp_le_f4 = np.asarray(vp, dtype="<f4", order="F")

    with open(output_raw, "wb") as f:
        f.write(vp_le_f4.tobytes(order="F"))


def write_model_info(output_info, nx, nz, dx, dz, raw_file):
    """
    Write a small text metadata file next to the raw binary file.
    """

    with open(output_info, "w") as f:
        f.write(f"raw_file = {raw_file}\n")
        f.write("dtype = float32\n")
        f.write("endian = little\n")
        f.write("order = column-based\n")
        f.write("array_shape = (nz, nx)\n")
        f.write(f"nx = {nx}\n")
        f.write(f"nz = {nz}\n")
        f.write(f"dx = {dx}\n")
        f.write(f"dz = {dz}\n")


# -----------------------------
# COMMON-SHOT GATHER EXPORT
# -----------------------------

def get_shot_key(header, shot_key_name):
    """
    Return the shot identifier from a SEG-Y trace header.
    """

    if shot_key_name == "FieldRecord":
        return header[segyio.TraceField.FieldRecord]

    if shot_key_name == "EnergySourcePoint":
        return header[segyio.TraceField.EnergySourcePoint]

    if shot_key_name == "SourceX":
        return header[segyio.TraceField.SourceX]

    raise ValueError(f"Unsupported shot key: {shot_key_name}")


def pack_su_trace_header(
    tracl,
    tracr,
    fldr,
    tracf,
    sx,
    gx,
    ns,
    dt_us,
    offset=0,
    scalco=1,
):
    """
    Create a 240-byte SU trace header.

    SU uses the SEG-Y trace header layout, but without the 3200-byte textual
    header or 400-byte binary header.

    This writes a minimal little-endian SU header.
    """

    h = bytearray(240)

    # Bytes are 1-based in SEG-Y/SU docs.
    # struct.pack_into uses 0-based offsets.

    # int32 fields
    struct.pack_into("<i", h, 0, tracl)      # bytes 1-4: trace sequence number within line
    struct.pack_into("<i", h, 4, tracr)      # bytes 5-8: trace sequence number within reel
    struct.pack_into("<i", h, 8, fldr)       # bytes 9-12: field record number
    struct.pack_into("<i", h, 12, tracf)     # bytes 13-16: trace number within field record
    struct.pack_into("<i", h, 36, offset)    # bytes 37-40: offset
    struct.pack_into("<i", h, 72, sx)        # bytes 73-76: source x
    struct.pack_into("<i", h, 80, gx)        # bytes 81-84: group x

    # int16 fields
    struct.pack_into("<h", h, 68, scalco)    # bytes 69-70: coordinate scalar
    struct.pack_into("<h", h, 114, ns)       # bytes 115-116: number of samples
    struct.pack_into("<h", h, 116, dt_us)    # bytes 117-118: sample interval in microseconds

    return h


def infer_dt_us(segyfile):
    """
    Infer sample interval in microseconds.
    SEG-Y binary header stores Interval in microseconds.
    """

    try:
        dt = int(segyfile.bin[segyio.BinField.Interval])
        if dt > 0:
            return dt
    except Exception:
        pass

    samples = np.asarray(segyfile.samples, dtype=np.float64)

    if len(samples) > 1:
        ds = float(np.median(np.diff(samples)))

        # If sample axis is in seconds, convert to microseconds.
        # If already in milliseconds or microseconds, user may override with --dt-us.
        if ds < 1.0:
            return int(round(ds * 1_000_000.0))
        else:
            return int(round(ds))

    return 1000


def extract_common_shot_gathers_to_su(
    input_segy,
    output_dir,
    shot_key_name="FieldRecord",
    dt_us_override=None,
):
    """
    Extract common-shot gathers and write each shot to an SU file.

    Output files:
        output_dir/shot_<shot_id>_seismogram_p.su
    """

    os.makedirs(output_dir, exist_ok=True)

    shot_to_indices = defaultdict(list)

    with segyio.open(input_segy, "r", ignore_geometry=True) as f:
        f.mmap()

        ntr = f.tracecount
        ns = len(f.samples)

        if ns > 32767:
            raise ValueError(
                "SU stores ns as a signed 16-bit integer in the trace header. "
                f"This file has ns = {ns}, which is too large for standard SU."
            )

        if dt_us_override is not None:
            dt_us = int(dt_us_override)
        else:
            dt_us = infer_dt_us(f)

        if dt_us > 32767:
            print(
                f"Warning: dt_us = {dt_us} is larger than 32767. "
                "This may not fit in the standard SU short integer header field."
            )

        print(f"Input seismic SEG-Y: {input_segy}")
        print(f"Number of traces: {ntr}")
        print(f"Samples per trace: {ns}")
        print(f"dt_us: {dt_us}")
        print(f"Grouping traces by: {shot_key_name}")

        # First pass: group trace indices by shot id
        for itr in range(ntr):
            header = f.header[itr]
            shot_id = get_shot_key(header, shot_key_name)
            shot_to_indices[shot_id].append(itr)

        print(f"Number of shots found: {len(shot_to_indices)}")

        # Second pass: write one SU file per shot
        global_trace_counter = 1

        for shot_id in sorted(shot_to_indices.keys()):
            indices = shot_to_indices[shot_id]
            output_su = os.path.join(
                output_dir,
                f"shot_{shot_id}_seismogram_p.su",
            )

            with open(output_su, "wb") as fout:
                for local_trace_number, itr in enumerate(indices, start=1):
                    header = f.header[itr]

                    try:
                        sx = int(header[segyio.TraceField.SourceX])
                    except Exception:
                        sx = 0

                    try:
                        gx = int(header[segyio.TraceField.GroupX])
                    except Exception:
                        gx = 0

                    try:
                        offset = int(header[segyio.TraceField.offset])
                    except Exception:
                        offset = gx - sx

                    try:
                        scalco = int(header[segyio.TraceField.SourceGroupScalar])
                    except Exception:
                        scalco = 1

                    su_header = pack_su_trace_header(
                        tracl=global_trace_counter,
                        tracr=global_trace_counter,
                        fldr=int(shot_id),
                        tracf=local_trace_number,
                        sx=sx,
                        gx=gx,
                        ns=ns,
                        dt_us=dt_us,
                        offset=offset,
                        scalco=scalco,
                    )

                    trace = np.asarray(f.trace[itr], dtype="<f4")

                    fout.write(su_header)
                    fout.write(trace.tobytes(order="C"))

                    global_trace_counter += 1

            print(f"Wrote {output_su} with {len(indices)} traces")


# -----------------------------
# MAIN
# -----------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Export a SEG-Y VP model to little-endian float32 raw binary "
            "and extract common-shot gathers from a SEG-Y file to SU files."
        )
    )

    parser.add_argument(
        "--vp-segy",
        required=True,
        help="Input SEG-Y file containing the VP model.",
    )

    parser.add_argument(
        "--vp-raw",
        required=True,
        help="Output raw binary VP model file.",
    )

    parser.add_argument(
        "--data-segy",
        required=True,
        help="Input SEG-Y file containing all seismic traces.",
    )

    parser.add_argument(
        "--csg-dir",
        default="./csg",
        help="Output directory for common-shot SU files. Default: ./csg",
    )

    parser.add_argument(
        "--shot-key",
        default="FieldRecord",
        choices=["FieldRecord", "EnergySourcePoint", "SourceX"],
        help=(
            "SEG-Y header used to group traces into shots. "
            "Default: FieldRecord"
        ),
    )

    parser.add_argument(
        "--dt-us",
        type=int,
        default=None,
        help=(
            "Override sample interval in microseconds for SU headers. "
            "If omitted, it is inferred from the SEG-Y binary header."
        ),
    )

    args = parser.parse_args()

    # 1. Export VP model
    vp, nx, nz, dx, dz = read_vp_model_segy(args.vp_segy)

    export_raw_binary(vp, args.vp_raw)

    info_file = args.vp_raw + ".txt"
    write_model_info(
        output_info=info_file,
        nx=nx,
        nz=nz,
        dx=dx,
        dz=dz,
        raw_file=args.vp_raw,
    )

    print("Finished VP model export")
    print(f"VP raw file: {args.vp_raw}")
    print(f"VP metadata: {info_file}")
    print(f"VP shape: nz, nx = {nz}, {nx}")
    print(f"dx = {dx}")
    print(f"dz = {dz}")

    # 2. Extract common-shot gathers
    extract_common_shot_gathers_to_su(
        input_segy=args.data_segy,
        output_dir=args.csg_dir,
        shot_key_name=args.shot_key,
        dt_us_override=args.dt_us,
    )

    print("Finished common-shot gather extraction")


if __name__ == "__main__":
    main()
