"""Put two phases into one *common* (near-commensurate) supercell.

    python3 tools/build_common_cell.py --a A7.vasp --b BII.vasp --out . \
        --ma "-2 0 1, 0 2 0, -1 2 0" --mb "-2 -2 1, -1 -1 -1, -1 -2 -1"

Why this exists
---------------
vc-NEB interpolates the *cell* linearly, so the two endpoints must be
described by lattices of the same shape (same number of atoms, same
orientation, near-identical lattice vectors).  If you feed, say, the R-3m
hexagonal conventional cell of A7 (3 atoms) next to the C2/m conventional cell
of Bi-II (4 atoms), the interpolation is meaningless.  This script

  1. takes each input, rebuilds its conventional cell with spglib;
  2. expands it with an integer supercell matrix M (det M = #copies);
  3. Niggli-reduces the resulting lattice, so the basis is not near-parallel;
  4. rigidly rotates the second cell onto the first (SVD, Kabsch);
  5. finds the atom-to-atom mapping that minimises the total displacement
     (linear_sum_assignment over MIC distances, with a search over the origin
     shift);
  6. writes two POSCARs with the same atom count and the same order, plus a
     report of the residual strain and the per-atom displacement.

Pick M_A and M_B so that det(M_A) * Z_A == det(M_B) * Z_B (Z = atoms in the
conventional cell) and so that M_B^-1 M_A is close to the identity.  Searching
for that pair is a small integer problem; brute force over dets up to ~6 is
usually enough.

Needs spglib and scipy (not needed by the main package).
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from ase import Atoms
from ase.build import make_supercell
from ase.io import read, write
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation


def parse_matrix(text: str) -> np.ndarray:
    """'-2 0 1, 0 2 0, -1 2 0' -> 3x3 int array."""
    rows = [r for r in text.replace(";", ",").split(",") if r.strip()]
    M = np.array([[int(x) for x in r.split()] for r in rows], dtype=int)
    if M.shape != (3, 3):
        raise SystemExit(f"supercell matrix must be 3x3, got {M.shape}")
    return M


def conventional(atoms: Atoms, symprec: float = 1e-3) -> tuple[Atoms, str]:
    import spglib

    cell = (atoms.cell[:], atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    std = spglib.standardize_cell(cell, to_primitive=False, no_idealize=True,
                                  symprec=symprec)
    if std is None:
        raise SystemExit("spglib could not symmetrise the input")
    lat, pos, nums = std
    out = Atoms(numbers=nums, scaled_positions=pos, cell=lat, pbc=True)
    ds = spglib.get_symmetry_dataset((lat, pos, nums), symprec=symprec)
    return out, ds.international


def supercell(atoms: Atoms, M: np.ndarray) -> Atoms:
    a = make_supercell(atoms, M)
    import spglib

    H = np.array(a.cell)
    Hr = np.array(spglib.niggli_reduce(H))
    T = Hr @ np.linalg.inv(H)
    if np.abs(T - np.round(T)).max() > 1e-6:
        raise SystemExit("niggli reduction did not give an integer matrix")
    frac = a.get_scaled_positions() @ H @ np.linalg.inv(Hr)
    frac %= 1.0
    a = Atoms(numbers=a.get_atomic_numbers(), scaled_positions=frac,
              cell=Hr, pbc=True)
    keep: list[int] = []
    for i in range(len(a)):
        d = a.get_distances(i, keep, mic=True) if keep else np.array([])
        if not len(d) or d.min() > 1e-3:
            keep.append(i)
    if len(keep) != len(a):
        raise SystemExit(f"lost atoms in the reduced cell: {len(keep)}/{len(a)}")
    order = np.lexsort((a.positions[:, 2], a.positions[:, 1], a.positions[:, 0]))
    return Atoms(numbers=a.get_atomic_numbers()[order],
                 positions=a.get_positions()[order],
                 cell=a.get_cell(), pbc=True)


def nn_distance(atoms: Atoms) -> float:
    d = atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, np.inf)
    return float(d.min())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--a", required=True, help="initial phase (any format ASE reads)")
    ap.add_argument("--b", required=True, help="final phase")
    ap.add_argument("--out", default=".", help="output directory")
    ap.add_argument("--ma", required=True, help="supercell matrix of --a, e.g. '2 0 0, 0 2 0, 0 0 1'")
    ap.add_argument("--mb", required=True, help="supercell matrix of --b")
    ap.add_argument("--name-a", default="initial_common")
    ap.add_argument("--name-b", default="final_common")
    args = ap.parse_args(argv)

    M_A, M_B = parse_matrix(args.ma), parse_matrix(args.mb)
    a_raw, b_raw = read(args.a), read(args.b)
    conv_a, sg_a = conventional(a_raw)
    conv_b, sg_b = conventional(b_raw)

    print(f"A: input {len(a_raw)} atoms -> conventional {len(conv_a)} ({sg_a})"
          f"  V={conv_a.get_volume():.3f}  "
          f"V/atom={conv_a.get_volume() / len(conv_a):.2f}")
    print(f"B: input {len(b_raw)} atoms -> conventional {len(conv_b)} ({sg_b})"
          f"  V={conv_b.get_volume():.3f}  "
          f"V/atom={conv_b.get_volume() / len(conv_b):.2f}")

    sa = supercell(conv_a, M_A)
    sb = supercell(conv_b, M_B)
    print(f"\nA supercell: {len(sa)} atoms  V={sa.get_volume():.3f}  "
          f"det(M)={int(round(np.linalg.det(M_A)))}")
    print(f"B supercell: {len(sb)} atoms  V={sb.get_volume():.3f}  "
          f"det(M)={int(round(np.linalg.det(M_B)))}")
    if len(sa) != len(sb):
        raise SystemExit("the two supercells have a different number of atoms; "
                         "change M_A / M_B")

    # ---- rigid rotation of B onto A -------------------------------------
    Ha, Hb = np.array(sa.cell), np.array(sb.cell)
    rot = Rotation.align_vectors(Ha.reshape(-1, 3), Hb.reshape(-1, 3))[0]
    sb.cell = (rot.as_matrix() @ Hb.T).T
    sb.positions = (rot.as_matrix() @ sb.positions.T).T
    Ha, Hb = np.array(sa.cell), np.array(sb.cell)

    F = Hb @ np.linalg.inv(Ha)
    strain = np.linalg.norm(np.linalg.svd(F)[1] - 1) / np.sqrt(3) * 100
    print(f"\ncell A: a,b,c={np.round(np.linalg.norm(Ha, axis=1), 3)}  "
          f"angles={np.round(angles(Ha), 2)}")
    print(f"cell B: a,b,c={np.round(np.linalg.norm(Hb, axis=1), 3)}  "
          f"angles={np.round(angles(Hb), 2)}")
    print(f"residual strain B->A: {strain:.2f} %   "
          f"(isotropic part {100 * (np.cbrt(1 / np.linalg.det(F)) - 1):+.2f} %)")
    if strain > 15:
        print("  WARNING: more than 15 % residual strain -- look for a better "
              "pair of supercell matrices")

    # ---- atom matching ---------------------------------------------------
    metric = 0.5 * (Ha + Hb)
    fa = sa.get_scaled_positions()
    fb = sb.get_scaled_positions()
    n = len(fa)

    def dist_matrix(shift):
        d = fa[:, None, :] - (fb + shift)[None, :, :]
        d -= np.round(d)
        return np.linalg.norm(d @ metric, axis=2)

    best = None
    for j in range(n):
        shift = fa[0] - fb[j]
        D = dist_matrix(shift)
        r, c = linear_sum_assignment(D)
        cost = D[r, c].sum()
        if best is None or cost < best[0]:
            best = (cost, shift, c)
    cost, shift, perm = best
    D = dist_matrix(shift)
    disp = np.diag(D[:, perm])
    print(f"\nmatching: total |dr| = {cost:.3f} A  (mean {cost / n:.3f} A)")
    print("  per-atom displacement (A):", np.round(disp, 2))
    print(f"  worst pair: {disp.max():.2f} A")

    sb2 = Atoms(numbers=sb.get_atomic_numbers()[perm],
                positions=(sb.get_scaled_positions()[perm] + shift) @ Hb,
                cell=Hb, pbc=True)
    sa = Atoms(numbers=sa.get_atomic_numbers(), scaled_positions=fa,
               cell=Ha, pbc=True)
    sb2.wrap()

    import spglib
    for tag, at in (("A", sa), ("B", sb2)):
        ds = spglib.get_symmetry_dataset(
            (at.cell[:], at.get_scaled_positions(), at.get_atomic_numbers()),
            symprec=1e-3)
        print(f"  {tag}: {len(at)} atoms  V={at.get_volume():.2f}  "
              f"NN={nn_distance(at):.3f} A  "
              f"symmetry={ds.international if ds is not None else '?'}")

    os.makedirs(args.out, exist_ok=True)
    for name, at in ((args.name_a, sa), (args.name_b, sb2)):
        path = os.path.join(args.out, name + ".vasp")
        write(path, at, format="vasp", direct=True)
        print(f"  wrote {path}")
    return 0


def angles(H: np.ndarray) -> np.ndarray:
    out = []
    for i, j in ((1, 2), (0, 2), (0, 1)):
        c = np.dot(H[i], H[j]) / (np.linalg.norm(H[i]) * np.linalg.norm(H[j]))
        out.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return np.array(out)


if __name__ == "__main__":
    sys.exit(main())
