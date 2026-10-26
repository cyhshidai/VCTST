"""Reading / building the two endpoint structures."""

from __future__ import annotations

import os

import numpy as np

from .util import VcnebError


def read_structure(path: str):
    from ase.io import read

    if not os.path.isfile(path):
        raise VcnebError(f"structure file not found: {path}")
    try:
        atoms = read(path)
    except Exception as exc:
        raise VcnebError(
            f"cannot read '{path}' with ASE ({exc}).\n"
            f"  supported: POSCAR/CONTCAR, .cif, .xyz, .traj, .pwi, .vasp, ..."
        ) from exc
    if not isinstance(atoms, list):
        pass
    elif atoms:
        atoms = atoms[0]
    if len(atoms) == 0:
        raise VcnebError(f"'{path}' contains no atoms")
    return atoms


def build_structure(spec: dict):
    """Tiny builder so people can try things without any input file."""
    from ase.build import bulk

    if not isinstance(spec, dict):
        raise VcnebError("structure.build_* must be a mapping, e.g. "
                         "{formula: Cu, crystal: fcc, a: 3.6}")
    spec = dict(spec)
    formula = spec.pop("formula", None)
    if not formula:
        raise VcnebError("structure.build_* needs 'formula'")
    crystal = spec.pop("crystal", None) or spec.pop("crystalstructure", "fcc")
    a = spec.pop("a", None)
    supercell = spec.pop("supercell", None)
    vacuum = spec.pop("vacuum", None)
    try:
        atoms = bulk(formula, crystalstructure=crystal, a=a, **spec)
    except Exception as exc:
        raise VcnebError(
            f"cannot build {formula} in the {crystal} structure ({exc})") from exc
    if supercell:
        atoms = atoms.repeat(supercell)
    if vacuum:
        atoms.center(vacuum=float(vacuum), axis=2)
    return atoms


def _match_endpoints(initial, final, log=None):
    """Reorder `final` so atom i corresponds to atom i of `initial`.

    For homonuclear systems the file order of the two structures is arbitrary;
    without this step a linear interpolation happily makes two atoms collide
    and the whole NEB is nonsense.  Within one species we solve a small
    assignment problem on minimum-image distances.
    """
    from ase import Atoms
    from ase.geometry import find_mic

    sym0 = list(initial.get_chemical_symbols())
    sym1 = list(final.get_chemical_symbols())
    if sorted(sym0) != sorted(sym1):
        raise VcnebError(
            "initial and final structures contain different elements:\n"
            f"  initial: {sorted(set(sym0))}\n"
            f"  final  : {sorted(set(sym1))}")

    pos0 = np.array(initial.get_positions())
    pos1 = np.array(final.get_positions())
    order = [None] * len(sym0)
    for species in sorted(set(sym0)):
        idx0 = [i for i, s in enumerate(sym0) if s == species]
        idx1 = [i for i, s in enumerate(sym1) if s == species]
        if len(idx0) == 1:
            order[idx0[0]] = idx1[0]
            continue
        diff = pos1[idx1][None, :, :] - pos0[idx0][:, None, :]
        # find_mic() only takes 2D arrays in some ASE versions
        flat = diff.reshape(-1, 3)
        mic, _ = find_mic(flat, initial.cell, initial.pbc)
        mic = np.asarray(mic).reshape(diff.shape)
        cost = (mic ** 2).sum(axis=2)
        try:
            from scipy.optimize import linear_sum_assignment

            rows, cols = linear_sum_assignment(cost)
            mapping = {r: c for r, c in zip(rows, cols)}
        except ImportError:  # greedy fallback
            mapping, taken = {}, set()
            for r in np.argsort(cost.min(axis=1)):
                for c in np.argsort(cost[r]):
                    if c not in taken:
                        mapping[r] = c
                        taken.add(c)
                        break
        for a, j0 in enumerate(idx0):
            order[j0] = idx1[mapping[a]]

    if order == list(range(len(sym0))):
        return final
    total = float(np.sum(np.linalg.norm(pos1[order] - pos0, axis=1)))
    if log:
        log(f"  final structure re-ordered to match the initial atom order "
            f"(total displacement {total:.3f} A)")
    out = Atoms(symbols=sym0, positions=pos1[order],
                cell=final.cell, pbc=final.pbc)
    if final.has("initial_magmoms"):
        out.set_initial_magnetic_moments(
            np.asarray(final.get_initial_magnetic_moments())[order])
    return out


def get_endpoints(cfg, log=None):
    """Return (initial, final) Atoms, checked for compatibility."""
    s = cfg.structure
    if s.get("initial"):
        initial = read_structure(s.initial)
    else:
        initial = build_structure(s.get("build_initial"))
    if s.get("final"):
        final = read_structure(s.final)
    else:
        final = build_structure(s.get("build_final"))

    if len(initial) != len(final):
        raise VcnebError(
            f"the two endpoints have a different number of atoms "
            f"({len(initial)} vs {len(final)}); a NEB path needs identical "
            f"atom counts")
    final = _match_endpoints(initial, final, log=log)

    for name, atoms in (("initial", initial), ("final", final)):
        vol = float(np.linalg.det(atoms.cell.array))
        if vol <= 0:
            raise VcnebError(
                f"{name} structure has a non-positive cell volume ({vol:.4f}); "
                f"check the cell / that you did not hand it a left-handed cell")
    return initial, final
