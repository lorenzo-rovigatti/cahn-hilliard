#!/usr/bin/env python3
"""
Compute circularly averaged structure factor S(q) for each snapshot
of a 2D Cahn–Hilliard trajectory and extract peak height and position.

Input snapshot format:
  # step = ..., t = <float>, size = LxL
  followed by L rows of L floating values

Outputs (in output directory):
  sq_<time>.dat   : q   S(q)
  sq_max.dat      : time   S_max(time)
  q_max.dat       : time   q_max(time)

Assumptions:
  dx = 1
  periodic boundaries
"""

import argparse
import os
import re
import numpy as np


HEADER_RE = re.compile(
    r"t\s*=\s*([0-9eE+\-\.]+).*?size\s*=\s*(\d+)\s*x\s*(\d+)"
)


def sanitize_time(t: float, ndp: int = 6) -> str:
    """Convert time to a filename-safe token."""
    s = f"{t:.{ndp}f}"
    return s.replace(".", "p")


def snapshot_reader(path):
    """Yield (t, L, field) for each snapshot."""
    with open(path, "r") as f:
        while True:
            line = f.readline()
            if not line:
                return
            line = line.strip()
            if not line or not line.startswith("#"):
                continue

            m = HEADER_RE.search(line)
            if not m:
                raise ValueError(f"Cannot parse header:\n{line}")

            t = float(m.group(1))
            Lx = int(m.group(2))
            Ly = int(m.group(3))
            if Lx != Ly:
                raise ValueError("Only square grids supported.")
            L = Lx

            data = []
            for _ in range(L):
                row = f.readline()
                if not row:
                    raise EOFError("Unexpected EOF while reading snapshot.")
                vals = row.split()
                if len(vals) != L:
                    raise ValueError("Wrong number of columns in snapshot.")
                data.append([float(x) for x in vals])

            yield t, L, np.asarray(data, dtype=np.float64)


def radial_average_Sq(field, dx=1.0, subtract_mean=True):
    """
    Compute circularly averaged structure factor S(q).
    """
    L = field.shape[0]
    N = L * L

    phi = field - field.mean() if subtract_mean else field

    F = np.fft.fft2(phi)
    Sk = np.abs(F) ** 2 / N

    k = 2.0 * np.pi * np.fft.fftfreq(L, d=dx)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    K = np.sqrt(KX**2 + KY**2)

    dq = 2.0 * np.pi / (L * dx)
    shell = np.rint(K / dq).astype(np.int64)

    shell_flat = shell.ravel()
    Sk_flat = Sk.ravel()

    max_shell = shell_flat.max()
    sums = np.bincount(shell_flat, weights=Sk_flat, minlength=max_shell + 1)
    counts = np.bincount(shell_flat, minlength=max_shell + 1)

    Sq = np.zeros_like(sums)
    mask = counts > 0
    Sq[mask] = sums[mask] / counts[mask]

    q_vals = np.arange(max_shell + 1) * dq
    return q_vals, Sq


def find_peak(q, Sq, exclude_q0=True):
    """
    Return (q_max, S_max) from S(q).
    """
    if exclude_q0:
        mask = q > 0
        q = q[mask]
        Sq = Sq[mask]

    idx = np.argmax(Sq)
    return q[idx], Sq[idx]


def write_sq(outdir, t, q, Sq, ndp=6):
    os.makedirs(outdir, exist_ok=True)
    t_tok = sanitize_time(t, ndp)
    path = os.path.join(outdir, f"sq_{t_tok}.dat")
    np.savetxt(path, np.column_stack([q, Sq]), header="# q   S(q)")
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("traj", help="Trajectory file")
    ap.add_argument("-o", "--outdir", default="sq_out")
    ap.add_argument("--dx", type=float, default=1.0)
    ap.add_argument("--no-subtract-mean", action="store_true")
    ap.add_argument("--time-ndp", type=int, default=6)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    sqmax_path = os.path.join(args.outdir, "sq_max.dat")
    qmax_path = os.path.join(args.outdir, "q_max.dat")

    with open(sqmax_path, "w") as f_smax, open(qmax_path, "w") as f_qmax:
        f_smax.write("# time   S_max(time)\n")
        f_qmax.write("# time   q_max(time)\n")

        n = 0
        for t, L, field in snapshot_reader(args.traj):
            q, Sq = radial_average_Sq(
                field,
                dx=args.dx,
                subtract_mean=not args.no_subtract_mean,
            )

            write_sq(args.outdir, t, q, Sq, ndp=args.time_ndp)

            q_max, S_max = find_peak(q, Sq, exclude_q0=True)

            f_smax.write(f"{t:.8e} {S_max:.8e}\n")
            f_qmax.write(f"{t:.8e} {q_max:.8e}\n")

            n += 1
            print(
                f"[{n}] t={t:.4f}  "
                f"S_max={S_max:.3e}  q_max={q_max:.3e}"
            )

    print("\nDone.")
    print(f"  → {sqmax_path}")
    print(f"  → {qmax_path}")


if __name__ == "__main__":
    main()

