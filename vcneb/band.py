"""Variable-cell nudged elastic band.

The band lives in a joint coordinate space ``x = [cell(3x3) | positions(Nx3)]``
measured in Angstrom, so cell and atoms are optimised at the same time
(that is what makes it a *vc*-NEB).

Physics summary
---------------
* generalized force per image:  ``f = [-dE/dh | -dE/dr]``, i.e. the atomic
  forces plus the stress-derived force on the cell,
* tangent ``tau`` from the neighbouring images (central difference),
* ``f_perp = f - (f.tau) tau`` -- the part that actually lowers the barrier,
* spring force along ``tau`` keeps the images spread out,
* climbing image removes the spring and flips the parallel part, so the top
  image climbs down to the saddle point.
"""

from __future__ import annotations

import copy
import itertools
import json
import math
import os
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
from ase.calculators.calculator import Calculator

from .util import VcnebError, ensure_dir, fmt_duration, jdump, jload


# --------------------------------------------------------------------- helpers
#: 1 GPa in eV/Angstrom^3
GPA_TO_EV_A3 = 1.0 / 160.2176634


class EnthalpyCalculator(Calculator):
    """Add ``P*V`` to the energy and ``P`` to the diagonal of the stress.

    With that offset the plain "minimise E, relax the cell against the stress"
    machinery optimises the enthalpy ``H = E + P*V`` instead:

    * ``dH/dh = -V*(sigma + P*I)*h^-T`` -- exactly what
      :func:`cell_force_from_stress` computes when the stress carries the
      extra ``P``,
    * ``H = E + P*V`` is what the logger/plot/barrier report.

    The engine input itself must stay at zero pressure (`press` in QE,
    `PRESSURE` in CP2K, ...) -- see README, otherwise the pressure is counted
    twice.
    """

    implemented_properties = ["energy", "free_energy", "forces", "stress"]

    def __init__(self, calc, pressure_ev_a3: float):
        super().__init__()
        self.calc = calc
        self.pressure = float(pressure_ev_a3)  # eV/Angstrom^3

    def calculate(self, atoms=None, properties=("energy",),
                  system_changes=None):
        inner = list(getattr(self.calc, "implemented_properties", properties))
        wanted = [p for p in
                  ("energy", "free_energy", "forces", "stress") if p in inner]
        self.calc.calculate(atoms, wanted or list(properties), system_changes)
        res = dict(self.calc.results)
        vol = abs(float(np.linalg.det(atoms.get_cell())))
        shift = self.pressure * vol
        for key in ("energy", "free_energy"):
            if key in res:
                res[key] = float(res[key]) + shift
        if "stress" in res:
            s = np.array(res["stress"], dtype=float)
            s[:3] += self.pressure  # ASE Voigt: xx, yy, zz first
            res["stress"] = s
        self.results = res


def voigt_to_matrix(stress6) -> np.ndarray:
    xx, yy, zz, yz, xz, xy = np.asarray(stress6, dtype=float)
    return np.array([[xx, xy, xz],
                     [xy, yy, yz],
                     [xz, yz, zz]])


def cell_force_from_stress(stress: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Generalized force conjugate to the cell matrix (eV / Angstrom).

    ``-V*stress*cell^-T`` is the derivative of the energy with respect to the
    cell matrix at *fixed fractional* coordinates, i.e. the cell relaxes
    against the stress while the atoms follow the cell.  That is the standard
    variable-cell convention (the atoms are relaxed by their own Cartesian
    forces at the same time), and it is what tests/test_cellforce.py verifies
    against finite differences.

    Do not "fix" this by adding the term ``-frac.T @ forces`` that a fully
    consistent Cartesian chart would need: with free atoms that term is of the
    order of the atom forces, it makes the cell follow the atoms instead of the
    stress, and on the EMT fcc->hcp test the band then ends up with sheared
    cells (~10 deg) and a 4x higher barrier.  See
    tests/test_cellforce_basis.py for the numbers.
    """
    vol = abs(float(np.linalg.det(cell)))
    return -vol * np.asarray(stress, dtype=float) @ np.linalg.inv(cell).T


def max_abs(array) -> float:
    arr = np.asarray(array, dtype=float)
    return float(np.max(np.abs(arr))) if arr.size else 0.0


def _cell_weight(cfg) -> float:
    """Weight of one Angstrom of cell change in the joint path metric."""
    explicit = cfg.neb.get("cell_weight")
    if explicit:
        return max(1e-6, float(explicit))
    cell_max_step = float(cfg.neb.get("cell_max_step") or 0.02)
    max_step = float(cfg.neb.get("max_step") or 0.2)
    if cell_max_step <= 0:
        return 1.0
    return max(1.0, max_step / cell_max_step)


def _get_idpp():
    """ase moved idpp_interpolate to ase.mep in 3.23; support both."""
    try:
        from ase.mep import idpp_interpolate

        return idpp_interpolate
    except ImportError:
        pass
    try:
        from ase.mep.neb import idpp_interpolate

        return idpp_interpolate
    except ImportError:
        pass
    from ase.neb import idpp_interpolate as legacy

    return legacy


def _get_idpp_class():
    """The ASE calculator that evaluates the image dependent pair potential."""
    for module in ("ase.mep.neb", "ase.neb"):
        try:
            mod = __import__(module, fromlist=["IDPP"])
            return getattr(mod, "IDPP")
        except (ImportError, AttributeError):
            continue
    raise VcnebError(
        "this ASE build has no IDPP calculator (ase.mep.neb.IDPP). "
        "Set neb.interpolation: linear instead.")


@dataclass
class ImageState:
    cell: np.ndarray
    pos: np.ndarray
    energy: float | None = None
    forces: np.ndarray | None = None
    stress: np.ndarray | None = None
    #: population analysis (Mulliken) of this image, when the code prints one:
    #: {"atoms": [{index, element, kind, alpha, beta, charge, spin}, ...],
    #:  "total": {alpha, beta, charge, spin}}
    mulliken: dict | None = None

    def copy(self) -> "ImageState":
        return ImageState(self.cell.copy(), self.pos.copy(), self.energy,
                          None if self.forces is None else self.forces.copy(),
                          None if self.stress is None else self.stress.copy(),
                          copy.deepcopy(self.mulliken))


def _atoms_from_state(template, state: ImageState):
    atoms = template.copy()
    atoms.set_cell(state.cell, scale_atoms=False)
    atoms.set_positions(state.pos)
    return atoms


# ----------------------------------------------------------------- point/opt
class JointPoint:
    """Make the flat ``[cell | positions]`` vector look like an Atoms object.

    That way a well tested optimizer from :mod:`ase.optimize` (FIRE, BFGS, ...)
    can be run on the NEB coordinates: one optimizer per image, exactly like
    ASE's own NEB does it.
    """

    def __init__(self, x):
        self.x = np.asarray(x, dtype=float).ravel().copy()
        self.energy = None
        self.forces = np.zeros_like(self.x)

    def __len__(self):
        return len(self.x) // 3

    def get_positions(self):
        return self.x.reshape(-1, 3)

    def set_positions(self, positions):
        self.x = np.asarray(positions, dtype=float).ravel().copy()

    def get_forces(self):
        return np.asarray(self.forces, dtype=float).reshape(-1, 3)

    def get_potential_energy(self):
        return self.energy

    def iterimages(self):
        return iter([self])

    def __ase_optimizable__(self):
        return self


# ---------------------------------------------------------------------- Band
class Band:
    def __init__(self, initial, final, cfg, engine, log, dry_run=False):
        self.cfg = cfg
        self.engine = engine
        self.log = log
        self.dry_run = bool(dry_run)
        self.n_free = int(cfg.neb.images)
        self.variable_cell = bool(cfg.neb.variable_cell)
        self.template = initial.copy()
        self.n_atoms = len(initial)
        self.states: list[ImageState] = [
            ImageState(np.array(initial.get_cell()),
                       np.array(initial.get_positions())),
            ImageState(np.array(final.get_cell()),
                       np.array(final.get_positions())),
        ]
        self.step = 0
        self.history: list[dict] = []
        self.barrier = float("nan")
        self.step_times: list[float] = []
        self.pruned: list[int] = []
        # The joint coordinate is y = [w*cell | positions].  With w > 1 the cell
        # is "stiffer" than the atoms: one Angstrom of cell change counts as w
        # Angstrom of path length.  Derived from max_step/cell_max_step so that
        # a single step cap limits the cell to ~cell_max_step and the atoms to
        # ~max_step.
        self.cell_weight = _cell_weight(cfg)

    # ------------------------------------------------------------- geometry
    def n_images(self) -> int:
        return len(self.states)

    def atoms(self, idx: int):
        return _atoms_from_state(self.template, self.states[idx])

    def flat(self, idx: int) -> np.ndarray:
        """Joint coordinate y = [w*cell | positions] (see cell_weight)."""
        st = self.states[idx]
        n_cell = self.n_cell_components()
        if not n_cell:  # fixed cell: the path is atomic coordinates only
            return st.pos.ravel().copy()
        return np.concatenate([st.cell.ravel() * self.cell_weight,
                               st.pos.ravel()])

    def n_cell_components(self) -> int:
        return 9 if self.variable_cell else 0

    def unpack(self, idx: int, y: np.ndarray) -> None:
        n_cell = self.n_cell_components()
        if n_cell:
            self.states[idx].cell = (y[:n_cell] / self.cell_weight
                                     ).reshape(3, 3).copy()
        self.states[idx].pos = y[n_cell:].reshape(self.n_atoms, 3).copy()

    # ---------------------------------------------------------- evaluations
    def image_dir(self, idx: int) -> str:
        if idx == 0:
            return os.path.join(self.cfg.images_dir, "end_initial")
        if idx == self.n_images() - 1:
            return os.path.join(self.cfg.images_dir, "end_final")
        return os.path.join(self.cfg.images_dir, f"{idx:02d}")

    def _calculator(self, atoms, directory: str):
        """engine.calculator() + a one-time warning if the frozen k-mesh is
        too coarse for this particular cell (the cell moved further than the
        endpoints did)."""
        coarse = self.engine.kpoints_too_coarse(atoms)
        if coarse and not getattr(self, "_kpts_warned", False):
            self._kpts_warned = True
            self.log.warn(
                "frozen k-mesh "
                f"{list(self.engine.kpoints())} is too coarse for this cell "
                f"(axes {coarse}); the cell moved outside the range spanned by "
                "the two endpoints. Re-run with a denser calc.kpts, or switch "
                "to calc.kpts: {density: X} so the mesh is sized for the "
                "whole path")
        calc = self.engine.calculator(atoms, directory)
        pressure = float(self.cfg.neb.get("pressure") or 0.0) * GPA_TO_EV_A3
        if pressure:
            calc = EnthalpyCalculator(calc, pressure)
        return calc

    def evaluate(self, idx: int, retries: int = 2) -> None:
        from ase.calculators.calculator import CalculatorError

        st = self.states[idx]
        atoms = self.atoms(idx)
        # states may be unwrapped along the band (see _unwrap_images); the
        # codes are happier with atoms inside the cell, and energies/forces do
        # not depend on which periodic image an atom sits in.
        atoms.wrap()
        last_exc = None
        for attempt in range(retries + 1):
            directory = self.image_dir(idx)
            if attempt > 0:
                directory = f"{directory}_retry{attempt}"
                ensure_dir(directory)
            else:
                ensure_dir(directory)
            try:
                calc = self._calculator(atoms, directory)
                atoms.calc = calc
                energy = float(atoms.get_potential_energy())
                forces = np.array(atoms.get_forces())
                if not np.all(np.isfinite(forces)) or not math.isfinite(energy):
                    raise VcnebError(
                        f"non-finite energy/forces for image {idx}")
                stress6 = atoms.get_stress(voigt=True) if self.variable_cell else None
                results = getattr(calc, "results", None)
                mulliken = results.get("mulliken") if isinstance(results, dict) else None
                break
            except Exception as exc:  # noqa: BLE001 - retried below
                last_exc = exc
                if attempt >= retries:
                    raise VcnebError(
                        f"image {idx}: SCF/calculator failed after "
                        f"{retries + 1} attempts.\n"
                        f"  last error: {exc}\n"
                        f"  files: {self.image_dir(idx)}") from exc
                self.log.warn(
                    f"image {idx}: attempt {attempt + 1} failed "
                    f"({type(exc).__name__}: {str(exc).splitlines()[0][:160]}),"
                    f" retrying in {directory}")
        st.energy = energy
        st.forces = forces
        st.stress = None if stress6 is None else voigt_to_matrix(stress6)
        st.mulliken = mulliken
        if last_exc is not None and attempt > 0:
            self.log(f"  image {idx}: recovered on attempt {attempt + 1}")

    def evaluate_all(self, indices=None) -> None:
        indices = list(range(self.n_images())) if indices is None else list(indices)
        workers = max(1, int(self.cfg.neb.parallel_images))
        if workers == 1 or len(indices) <= 1:
            for idx in indices:
                self.evaluate(idx)
            return
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(self.evaluate, idx) for idx in indices]
            errors = []
            for fut, idx in zip(futures, indices):
                try:
                    fut.result()
                except Exception as exc:  # noqa: BLE001
                    errors.append(f"image {idx}: {exc}")
            if errors:
                raise VcnebError("some images failed:\n  " + "\n  ".join(errors))

    # ------------------------------------------------------------- endpoints
    def relax_endpoints(self) -> None:
        if not self.cfg.neb.relax_endpoints:
            return
        if self.dry_run:
            self.log("dry run: endpoints are used exactly as given "
                     "(no DFT relaxation, no electronic-structure call)")
            return
        from ase.filters import ExpCellFilter
        from ase.optimize import FIRE

        fmax = float(self.cfg.neb.endpoint_fmax)
        steps = int(self.cfg.neb.endpoint_steps)
        for idx, tag in ((0, "initial"), (self.n_images() - 1, "final")):
            atoms = self.atoms(idx)
            directory = self.image_dir(idx)
            ensure_dir(directory)
            atoms.calc = self._calculator(atoms, directory)
            self.log(f"relaxing {tag} endpoint (fmax={fmax} eV/A, "
                     f"cell={'free' if self.variable_cell else 'fixed'}) ...")
            logfile = os.path.join(self.cfg.logs_dir, f"relax_{tag}.log")
            target = ExpCellFilter(atoms) if self.variable_cell else atoms
            opt = FIRE(target, logfile=logfile, maxstep=float(self.cfg.neb.max_step))
            t0 = time.time()
            try:
                opt.run(fmax=fmax, steps=steps)
                converged = bool(opt.converged())
            except Exception as exc:  # noqa: BLE001
                raise VcnebError(
                    f"relaxation of the {tag} endpoint failed: {exc}\n"
                    f"  see {logfile}") from exc
            energy = float(atoms.get_potential_energy())
            forces = np.array(atoms.get_forces())
            stress6 = atoms.get_stress(voigt=True) if self.variable_cell else None
            results = getattr(atoms.calc, "results", None)
            mulliken = (results.get("mulliken")
                        if isinstance(results, dict) else None)
            self.states[idx] = ImageState(
                np.array(atoms.get_cell()), np.array(atoms.get_positions()),
                energy, forces,
                None if stress6 is None else voigt_to_matrix(stress6),
                mulliken)
            self.log(f"  {tag}: E={energy:.6f} eV  "
                     f"max|F|={max_abs(forces):.4f} eV/A  "
                     f"{'converged' if converged else 'NOT converged'}  "
                     f"({fmt_duration(time.time() - t0)})")
            if not converged:
                self.log.warn(
                    f"endpoint '{tag}' did not reach {fmax} eV/A in {steps} "
                    f"steps; the barrier will be less accurate. "
                    f"Increase neb.endpoint_steps.")
            self._save_structure(idx, os.path.join(
                self.cfg.workdir, f"endpoint_{tag}.vasp"))

    def _save_structure(self, idx: int, path: str) -> None:
        try:
            from ase.io import write

            write(path, self.atoms(idx), format="vasp")
        except Exception:  # noqa: BLE001 - never fail because of output
            pass

    # ---------------------------------------------------------- interpolation
    def _align_final_cell(self) -> None:
        """Re-index the final cell so that it lines up with the initial one.

        The same crystal can be written with its lattice vectors in any order
        and in either direction (``a1`` <-> ``a2``, ``a1`` -> ``-a1`` ...).
        The interpolated path is built in *fractional* coordinates, so two
        structures that are physically identical but indexed differently give
        a path that swaps the atoms into each other.  Here the 48 equivalent
        descriptions of the final cell are tried and the one closest to the
        initial cell is kept: cell strain (symmetric part of the deformation)
        plus the RMS atom displacement in Angstrom.

        This is a change of *description* only -- the final structure itself
        (atom positions, lattice) is untouched.
        """
        first, last = self.states[0], self.states[-1]
        c0 = np.asarray(first.cell, dtype=float)
        c1 = np.asarray(last.cell, dtype=float)
        if abs(float(np.linalg.det(c0))) < 1e-8 or abs(
                float(np.linalg.det(c1))) < 1e-8:
            return
        frac0 = np.asarray(first.pos, dtype=float) @ np.linalg.inv(c0)
        pos1 = np.asarray(last.pos, dtype=float)
        scale = float(np.mean(np.linalg.norm(c0, axis=1))) or 1.0
        best = None
        for perm in itertools.permutations(range(3)):
            for signs in itertools.product((1.0, -1.0), repeat=3):
                c1p = c1[list(perm)] * np.array(signs)[:, None]
                # F^T = c0^-1 c1p, strain = symmetric part - identity
                ft = np.linalg.solve(c0, c1p)
                strain = float(np.linalg.norm(0.5 * (ft + ft.T) - np.eye(3)))
                d = pos1 @ np.linalg.inv(c1p) - frac0
                d -= np.round(d)                      # minimum image
                disp = float(np.sqrt(np.mean(np.sum((d @ c0) ** 2, axis=1))))
                score = strain + disp / scale
                if best is None or score < best[0]:
                    best = (score, strain, disp, perm, signs, c1p)
        _, strain, disp, perm, signs, c1p = best
        if perm != (0, 1, 2) or any(s < 0 for s in signs):
            last.cell = c1p
            order = "".join("abc"[i] for i in perm)
            flipped = "".join("+" if s > 0 else "-" for s in signs)
            self.log(f"  final cell re-indexed to line up with the initial "
                     f"one: (a,b,c) -> ({order[0]},{order[1]},{order[2]}), "
                     f"signs {flipped}")
        if strain > 0.15 or disp > 2.5:
            self.log.warn(
                "the two endpoint cells do not line up (strain "
                f"{100 * strain:.0f} %, atom RMS displacement {disp:.2f} A) "
                "and no re-indexing can fix that.  The two structures are "
                "probably not the same supercell (e.g. 2x1x1 vs 3x1x1) or the "
                "atom order does not match, so the interpolated images may "
                "contain overlapping atoms; see neb.min_distance / neb.prune.")

    def interpolate(self) -> None:
        mode = str(self.cfg.neb.interpolation).lower()
        if bool(self.cfg.neb.get("align_cell", True)):
            self._align_final_cell()
        first, last = self.states[0], self.states[-1]
        n_free = self.n_free
        total = n_free + 2
        frac0 = np.linalg.solve(np.array(first.cell).T, first.pos.T).T
        frac1 = np.linalg.solve(np.array(last.cell).T, last.pos.T).T
        # minimum-image difference in fractional coordinates
        dfrac = frac1 - frac0
        dfrac -= np.round(dfrac)
        new_states = [first.copy()]
        for k in range(1, n_free + 1):
            t = k / (n_free + 1)
            cell = (1 - t) * np.array(first.cell) + t * np.array(last.cell)
            frac = frac0 + t * dfrac
            pos = frac @ cell
            new_states.append(ImageState(cell, pos))
        new_states.append(last.copy())
        self.states = new_states
        if mode == "idpp" and n_free > 0 and self.n_atoms > 1:
            try:
                self._idpp_relax()
            except Exception as exc:  # noqa: BLE001
                self.log.warn(f"IDPP interpolation failed ({exc}); "
                              f"keeping the linear path")
        elif mode == "idpp":
            self.log("  only one atom: IDPP has nothing to do, "
                     "keeping the linear path")
        self._unwrap_images()
        self._prune_close_images()
        total = self.n_images()
        self.log(f"initial path built: {total} images "
                 f"({self.n_free} movable), interpolation={mode}")
        dist = self.min_pair_distance()
        self.log(f"smallest interatomic distance on the path: {dist:.3f} A")
        try:
            n_cell = self.n_cell_components()
            steps = []
            for k in range(total - 1):
                d = self.flat(k + 1)[n_cell:] - self.flat(k)[n_cell:]
                steps.append(float(np.sqrt(np.mean(d ** 2))))
            self.log("path spacing (RMS atom displacement between neighbouring "
                     "images): "
                     + ", ".join(f"{s:.3f}" for s in steps) + " A")
            uneven = bool(steps) and max(steps) > 3.0 * max(min(steps), 1e-6)
            if uneven and self.pruned:
                self.log("  (uneven on purpose: images "
                         + ", ".join(str(i) for i in self.pruned)
                         + " were dropped, see neb.min_distance)")
            elif uneven:
                self.log.warn(
                    "the path spacing is very uneven (see the list above); "
                    "usually the atom order of the two endpoint structures "
                    "does not match. Check `endpoint_initial.vasp` and "
                    "`endpoint_final.vasp`.")
        except Exception:  # noqa: BLE001
            pass
        if dist < 1.2:
            self.log.warn(
                "atoms are very close on the initial path (< 1.2 A). "
                "The path may be wrong; check the atom correspondence or "
                "your endpoint structures.")

    def _image_min_distance(self, idx: int) -> float:
        d = self.atoms(idx).get_all_distances(mic=True)
        np.fill_diagonal(d, np.inf)
        return float(d.min())

    def _prune_close_images(self) -> None:
        """Throw away interpolated images in which two atoms overlap.

        An image with a 0.1 A contact cannot be evaluated (and no amount of
        IDPP fixes a cell that is being interpolated through a wrong shape),
        so such images are simply removed and the band continues with the
        remaining ones -- the same "discard bad candidates" idea as USPEX.
        The endpoints are never touched.
        """
        dmin = float(self.cfg.neb.get("min_distance") or 0.0)
        if dmin <= 0 or self.n_free < 1:
            return
        mode = self.cfg.neb.get("prune")
        mode = "auto" if mode is True else mode   # validate() normalises it
        if mode is False:
            return
        bad = [(i, self._image_min_distance(i))
               for i in range(1, self.n_images() - 1)]
        bad = [(i, d) for i, d in bad if d < dmin]
        if not bad:
            return
        detail = ", ".join(f"{d:.2f}" for _, d in bad)
        if mode == "error":
            raise VcnebError(
                f"{len(bad)} interpolated image(s) put atoms closer than "
                f"neb.min_distance={dmin} A (images "
                + ", ".join(str(i) for i, _ in bad)
                + f": {detail} A).\n"
                f"  Fix the endpoints (atom order / supercell must match) or "
                f"set neb.prune: auto to drop those images and continue.")
        bad_idx = {i for i, _ in bad}
        keep = [st for i, st in enumerate(self.states) if i not in bad_idx]
        if len(keep) - 2 < 1:
            raise VcnebError(
                f"every interpolated image is unphysical (closest contact "
                f"{detail} A < neb.min_distance={dmin} A), so there is no band "
                f"left.\n  The two endpoint structures almost certainly do not "
                f"correspond to each other (different supercell, atom order or "
                f"species order).  Increase neb.images only helps if the "
                f"endpoints are right.")
        self.pruned = [i for i, _ in bad]
        self.states = keep
        self.n_free = len(keep) - 2
        self.log.warn(
            "dropped images " + ", ".join(str(i) for i in self.pruned)
            + f" (closest contact {detail} A < neb.min_distance={dmin} A); "
            f"band now has {self.n_free} movable images")
        if self.n_free < 3:
            self.log.warn(
                f"only {self.n_free} movable image(s) left: a band that coarse "
                f"usually overestimates the barrier badly.  The real fix is to "
                f"make the endpoints commensurate (same supercell, same atom "
                f"order), not to prune images.")

    def _idpp_relax(self) -> None:
        """Image dependent pair potential, one image at a time.

        ASE's ``NEB``-based helper refuses cells that change along the band,
        which is exactly what a vc-NEB has.  The IDPP objective itself only
        needs distances, so each image can be optimised on its own with a
        fixed cell -- that is what this does.
        """
        from ase.optimize import MDMin

        idpp_cls = _get_idpp_class()
        images = [self.atoms(i) for i in range(self.n_images())]
        d1 = images[0].get_all_distances(mic=True)
        d2 = images[-1].get_all_distances(mic=True)
        n = self.n_images() - 1
        steps = int(self.cfg.neb.get("idpp_steps") or 200)
        fmax = float(self.cfg.neb.get("idpp_fmax") or 0.02)
        for i in range(1, self.n_images() - 1):
            atoms = images[i]
            atoms.calc = idpp_cls(d1 + i * (d2 - d1) / n, mic=True)
            MDMin(atoms, logfile=None).run(fmax=fmax, steps=steps)
            st = self.states[i]
            st.pos = np.array(atoms.get_positions())
        self.log("  IDPP pre-relaxation of the movable images done")

    def _unwrap_images(self) -> None:
        """Put neighbouring images into the same periodic image of every atom.

        IDPP (and any other minimum-image based pre-relaxation) is free to move
        an atom across a cell boundary.  The NEB tangent, the spring forces and
        the path spacing all use plain Cartesian differences, so a wrapped atom
        looks like a jump of one cell length (~5 A) and the band blows up.  Here
        every atom is unwrapped along the band,
        ``x_i = x_{i-1} + mic(x_i - x_{i-1})``, which describes the very same
        structures (periodic images) but makes the coordinate continuous.
        """
        from ase.geometry import find_mic

        if self.n_atoms < 2:
            return
        moved = 0
        for i in range(1, self.n_images()):
            prev, cur = self.states[i - 1], self.states[i]
            d = cur.pos - prev.pos
            mic, _ = find_mic(d.reshape(-1, 3), np.array(cur.cell), True)
            mic = np.asarray(mic, dtype=float).reshape(d.shape)
            if float(np.abs(d - mic).max()) < 1e-8:
                continue
            cur.pos = prev.pos + mic
            moved += 1
        if moved:
            self.log(f"  unwrapped {moved} image(s) so the path is continuous "
                     f"across the periodic boundary")

    def min_pair_distance(self) -> float:
        smallest = float("inf")
        for i in range(self.n_images()):
            atoms = self.atoms(i)
            d = atoms.get_all_distances(mic=True)
            np.fill_diagonal(d, np.inf)
            smallest = min(smallest, float(d.min()))
        return smallest

    # --------------------------------------------------------------- forces
    def neb_forces(self, ci_index: int | None = None
                   ) -> tuple[list[np.ndarray], dict]:
        """Return the NEB force for every image plus diagnostics."""
        n_img = self.n_images()
        xs = [self.flat(i) for i in range(n_img)]
        grad = []
        for i in range(n_img):
            st = self.states[i]
            g = -np.asarray(st.forces, dtype=float).ravel()
            if self.variable_cell:
                f_cell = cell_force_from_stress(st.stress, st.cell)
                # y = w*cell, so the force conjugate to y is F_cell / w
                g = np.concatenate([-f_cell.ravel() / self.cell_weight, g])
            grad.append(g)

        forces = [np.zeros_like(g) for g in grad]
        diag = {"f_par": [], "f_perp_atom": [], "f_perp_cell": [], "spring": []}
        n_cell = self.n_cell_components()
        k_spring = float(self.cfg.neb.spring_k)
        ci = ci_index
        for i in range(1, n_img - 1):
            prev, nxt = xs[i - 1], xs[i + 1]
            tau = nxt - prev
            tnorm = float(np.linalg.norm(tau))
            if tnorm < 1e-12:
                tau = xs[i] - xs[i - 1]
                tnorm = float(np.linalg.norm(tau))
            if tnorm < 1e-12:
                continue
            tau = tau / tnorm
            f = -grad[i]
            f_par = float(np.dot(f, tau))
            f_perp = f - f_par * tau
            spring = k_spring * (float(np.linalg.norm(nxt - xs[i]))
                                 - float(np.linalg.norm(xs[i] - prev)))
            if ci is not None and i == ci:
                forces[i] = f_perp - f_par * tau  # climb
            else:
                forces[i] = f_perp + spring * tau
            perp_atom = f_perp[n_cell:]
            perp_cell = f_perp[:n_cell] if n_cell else np.zeros(1)
            diag["f_par"].append(f_par)
            diag["f_perp_atom"].append(max_abs(perp_atom))
            diag["f_perp_cell"].append(max_abs(perp_cell))
            diag["spring"].append(abs(spring))
        forces[0] = np.zeros_like(grad[0])
        forces[-1] = np.zeros_like(grad[-1])
        return forces, diag

    # ------------------------------------------------------------------ run
    def _make_optimizers(self) -> dict:
        """One optimizer per movable image (sharing one would mix velocities)."""
        from ase.optimize import FIRE

        neb = self.cfg.neb
        n_cell = self.n_cell_components()
        n = n_cell + 3 * self.n_atoms
        self._points, self._opts = {}, {}
        for i in range(1, self.n_images() - 1):
            point = JointPoint(np.zeros(n))
            opt = FIRE(point, logfile=None, dt=float(neb.dt),
                       maxstep=float(neb.max_step),
                       dtmax=float(neb.dt_max))
            if getattr(self, "_optimizer_dt", None) and i in self._optimizer_dt:
                opt.dt = max(float(self._optimizer_dt[i]), 1e-4)
            opt.initialize()
            self._points[i] = point
            self._opts[i] = opt
        return self._opts

    def _optimizer_step(self, i: int, y: np.ndarray, f: np.ndarray) -> np.ndarray:
        point = self._points[i]
        point.x = np.array(y, dtype=float)
        point.forces = np.array(f, dtype=float)
        point.energy = self.states[i].energy
        self._opts[i].step(f=point.get_forces())
        return point.x

    def run(self) -> dict:
        neb = self.cfg.neb
        f_th = float(neb.f_th)
        max_steps = int(neb.max_steps)
        if self.step >= max_steps:
            raise VcnebError(
                f"this band is already at step {self.step} while max_steps="
                f"{max_steps}: there is nothing left to do.\n"
                f"  Raise neb.max_steps in the config (or pass e.g. "
                f"`--max-steps {self.step + 500}`) and run with --restart "
                f"again.")
        optimizers = self._make_optimizers()
        self.log.rule()
        self.log(f"starting vc-NEB: {self.n_free} movable images, "
                 f"f_th={f_th} eV/A, max_steps={max_steps}")
        self.log.rule()
        t_start = time.time()
        converged = False
        ci_index = None
        f_max_atom, f_cell_per_atom = getattr(
            self, "converged_forces", (float("nan"), float("nan")))
        for step in range(self.step + 1, max_steps + 1):
            self.step = step
            t0 = time.time()
            self.evaluate_all()
            energies = [st.energy for st in self.states]
            if step == 1 and max(energies) - min(energies) < 1e-9:
                self.log.warn(
                    "every image returned the same energy (spread < 1e-9 eV). "
                    "Either the potential energy surface really is flat or the "
                    "calculator is not doing anything -- check the image logs.")
            if neb.climbing_image:
                free_e = energies[1:-1]
                cand = 1 + int(np.argmax(free_e))
                if ci_index is None:
                    ci_index = cand
                else:
                    # Hysteresis: with nearly degenerate images the climbing
                    # image jumps back and forth every step and the band never
                    # settles.  Only move it when another image is clearly
                    # higher.
                    tol = max(1e-5, 1e-3 * abs(energies[cand] - energies[0]))
                    if energies[cand] > energies[ci_index] + tol:
                        ci_index = cand
            forces, diag = self.neb_forces(ci_index if neb.climbing_image else None)
            f_max_atom = max(diag["f_perp_atom"]) if diag["f_perp_atom"] else 0.0
            f_max_cell = max(diag["f_perp_cell"]) if diag["f_perp_cell"] else 0.0
            f_cell_per_atom = f_max_cell / max(1, self.n_atoms)
            barrier = max(energies[1:-1] or energies) - energies[0]
            reverse = max(energies[1:-1] or energies) - energies[-1]
            self.barrier = barrier
            self.step_times.append(time.time() - t0)
            self._record(step, energies, f_max_atom, f_cell_per_atom, ci_index)
            msg = (f"step {step:4d}  E[max]={max(energies):.5f}  "
                   f"barrier={barrier:.4f}  max|F_atom|={f_max_atom:.4f}  "
                   f"max|F_cell|/atom={f_cell_per_atom:.4f}  "
                   f"dt={min(o.dt for o in optimizers.values()):.3f}  "
                   f"({time.time() - t0:.0f}s)")
            if ci_index is not None:
                msg += f"  CI={ci_index}"
            self.log(msg)
            if f_max_atom < f_th and f_cell_per_atom < f_th:
                converged = True
                self.log(f"converged: max|F_perp| = {f_max_atom:.4f} eV/A "
                         f"< {f_th}")
                self.converged_forces = (f_max_atom, f_cell_per_atom)
                break
            self.converged_forces = (f_max_atom, f_cell_per_atom)
            xs = [self.flat(i) for i in range(self.n_images())]
            new_ys = []
            for i in range(1, self.n_images() - 1):
                new_ys.append((i, self._optimizer_step(i, xs[i], forces[i])))
            for i, new_y in new_ys:
                self.unpack(i, new_y)
            self._guard_geometry()
            self._optimizer_dt = {i: o.dt for i, o in optimizers.items()}
            self.save_state()
        else:
            self.log.warn(
                f"reached max_steps={max_steps} without converging "
                f"(last max|F_perp|={f_max_atom:.4f} eV/A, target {f_th}). "
                f"Re-run with --restart to continue.")

        self.converged = converged
        elapsed = time.time() - t_start
        self.log.rule()
        self.log(f"finished after {self.step} steps in {fmt_duration(elapsed)}")
        summary = self.write_outputs(converged, elapsed)
        return summary

    def _guard_geometry(self) -> None:
        for i in range(1, self.n_images() - 1):
            st = self.states[i]
            vol = float(np.linalg.det(st.cell))
            if vol <= 1e-3:
                raise VcnebError(
                    f"image {i} cell collapsed or became left-handed "
                    f"(volume {vol:.4f} A^3).  Reduce neb.cell_max_step "
                    f"or check the initial/final cells.")
        dist = self.min_pair_distance()
        if dist < 0.6:
            raise VcnebError(
                f"two atoms are only {dist:.2f} A apart in the band; the path "
                f"is unphysical.  Check the endpoint structures / the atom "
                f"order, and lower neb.max_step.")

    # ------------------------------------------------------------------ I/O
    def _record(self, step, energies, f_atom, f_cell, ci_index) -> None:
        self.history.append({
            "step": step,
            "energies": [float(e) for e in energies],
            "max_force_atom": float(f_atom),
            "max_force_cell_per_atom": float(f_cell),
            "climbing_image": ci_index,
        })
        path = os.path.join(self.cfg.workdir, "history.jsonl")
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(self.history[-1]) + "\n")
        if self.cfg.output.get("save_trajectory", True):
            self.write_trajectory()

    def write_trajectory(self) -> None:
        from ase.io import Trajectory

        path = os.path.join(self.cfg.workdir, "band.traj")
        traj = Trajectory(path, "w")
        for i in range(self.n_images()):
            atoms = self.atoms(i)
            atoms.wrap()  # states are unwrapped along the band, see above
            atoms.info["image"] = i
            atoms.info["energy"] = self.states[i].energy
            traj.write(atoms)
        traj.close()

    def save_state(self) -> None:
        payload = {
            "step": self.step,
            "variable_cell": self.variable_cell,
            "n_atoms": self.n_atoms,
            "n_movable": self.n_free,
            "pruned_images": list(self.pruned),
            "optimizer_dt": {str(i): float(dt) for i, dt in
                             (getattr(self, "_optimizer_dt", None) or {}).items()},
            "images": [{"cell": st.cell.tolist(), "pos": st.pos.tolist(),
                        "energy": st.energy}
                       for st in self.states],
        }
        jdump(payload, os.path.join(self.cfg.workdir, "state.json"), indent=0)

    def load_state(self, path: str) -> None:
        payload = jload(path)
        if not payload:
            raise VcnebError(f"no restart state at {path}")
        if payload.get("n_atoms") != self.n_atoms:
            raise VcnebError(
                f"restart state in {path} has {payload.get('n_atoms')} atoms, "
                f"the current setup has {self.n_atoms}")
        states = []
        for item in payload["images"]:
            states.append(ImageState(np.array(item["cell"], float),
                                     np.array(item["pos"], float),
                                     item.get("energy")))
        n_mov = payload.get("n_movable")
        n_mov = len(states) - 2 if n_mov is None else int(n_mov)
        if len(states) != n_mov + 2:
            raise VcnebError(
                f"restart state {path} is inconsistent: {len(states)} images "
                f"but n_movable={n_mov}.  Delete it and rebuild the band.")
        if n_mov != self.n_free:
            self.log.warn(
                f"the saved band has {n_mov} movable images, the config asks "
                f"for {self.n_free}; keeping the saved band.  The image count "
                f"cannot change in the middle of a run -- if you really want a "
                f"different path, delete {os.path.join(self.cfg.workdir, 'state.json')} "
                f"and start again.")
        self.n_free = n_mov
        self.pruned = [int(i) for i in payload.get("pruned_images") or []]
        self.states = states
        self.step = int(payload.get("step", 0))
        dts = payload.get("optimizer_dt") or {}
        self._optimizer_dt = {int(k): float(v) for k, v in dts.items()}
        self.log(f"restarting from step {self.step} "
                 f"({len(states)} images)")

    # -------------------------------------------------------------- outputs
    def energies(self) -> list[float]:
        return [st.energy for st in self.states]

    def write_outputs(self, converged: bool, elapsed: float) -> dict:
        energies = self.energies()
        e0, eN = energies[0], energies[-1]
        emax = max(energies)
        idx_max = int(np.argmax(energies))
        barrier = emax - e0
        reverse = emax - eN
        forces = [max_abs(st.forces) for st in self.states]
        rows = []
        for i, st in enumerate(self.states):
            total = (st.mulliken or {}).get("total") or {}
            rows.append({
                "image": i,
                "energy_eV": st.energy,
                "delta_eV": st.energy - e0,
                "max_force_eV_per_A": forces[i],
                "volume_A3": abs(float(np.linalg.det(st.cell))),
                # Mulliken totals; None when the code prints no population
                # analysis (EMT/LJ) or mulliken: false.
                "total_spin": total.get("spin"),
                "net_charge": total.get("charge"),
                "mulliken": st.mulliken,
            })
        # the population analysis is optional: only add its columns when at
        # least one image produced one, so old barrier.csv consumers keep working
        has_pop = any(r["total_spin"] is not None or r["mulliken"] for r in rows)

        def _num(value):
            return "" if value is None else f"{float(value):.6f}"

        csv_path = os.path.join(self.cfg.workdir, "barrier.csv")
        with open(csv_path, "w", encoding="utf-8") as fh:
            fh.write("image,energy_eV,delta_eV,max_force_eV_per_A,volume_A3")
            if has_pop:
                fh.write(",total_spin,net_charge")
            fh.write("\n")
            for row in rows:
                line = (f"{row['image']},{row['energy_eV']:.8f},"
                        f"{row['delta_eV']:.8f},"
                        f"{row['max_force_eV_per_A']:.6f},"
                        f"{row['volume_A3']:.6f}")
                if has_pop:
                    line += f",{_num(row['total_spin'])},{_num(row['net_charge'])}"
                fh.write(line + "\n")
        summary = {
            "engine": self.engine.name,
            "converged": bool(converged),
            "steps": self.step,
            "kpts": (list(self.engine.kpoints()) if self.engine.kpoints()
                     else "gamma"),
            "barrier_eV": barrier,
            "reverse_barrier_eV": reverse,
            "energy_initial_eV": e0,
            "energy_final_eV": eN,
            "energy_max_eV": emax,
            "image_max": idx_max,
            "f_th_eV_per_A": float(self.cfg.neb.f_th),
            "variable_cell": self.variable_cell,
            "n_atoms": self.n_atoms,
            "n_movable_images": self.n_free,
            "seconds": elapsed,
            "images": rows,
        }
        if getattr(self, "converged_forces", None):
            summary["final_max_force_atom_eV_per_A"] = float(
                self.converged_forces[0])
            summary["final_max_force_cell_per_atom_eV_per_A"] = float(
                self.converged_forces[1])
        jdump(summary, os.path.join(self.cfg.workdir, "summary.json"))
        self.log(f"barrier (forward) = {barrier:.4f} eV at image {idx_max}")
        self.log(f"barrier (reverse) = {reverse:.4f} eV")
        if has_pop:
            spins = ["?" if r["total_spin"] is None else f"{float(r['total_spin']):.2f}"
                     for r in rows]
            self.log("Mulliken total spin moment along the path "
                     f"(images 0..{len(rows) - 1}): " + "  ".join(spins))
            # Mulliken spin moments carry ~0.01 muB of basis-set noise, so
            # only a change beyond that is worth a warning
            vals = [float(r["total_spin"]) for r in rows
                    if r["total_spin"] is not None]
            if vals and max(vals) - min(vals) > 0.05:
                self.log.warn(
                    "the total spin moment changes along the path. With a "
                    "magnetization forced by the input (CP2K MULTIPLICITY "
                    "without smearing, QE tot_magnetization, ABACUS nupdown) "
                    "part of the barrier can be an occupation constraint "
                    "rather than the geometry. Compare the Mulliken table of "
                    "the endpoint and the climbing image.")
        if self.cfg.output.get("plot", True):
            try:
                plot_band(self.cfg.workdir, summary)
                self.log(f"plot written: {os.path.join(self.cfg.workdir,
                                                       'barrier.png')}")
            except Exception as exc:  # noqa: BLE001
                self.log.warn(f"could not draw the barrier plot: {exc}")
        return summary


# ---------------------------------------------------------------------- plot
def plot_band(workdir: str, summary: dict | None = None,
              out: str | None = None) -> str:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if summary is None:
        summary = jload(os.path.join(workdir, "summary.json"))
        if summary is None:
            raise VcnebError(f"no summary.json in {workdir}; run the NEB first")
    rows = summary["images"]
    idx = [r["image"] for r in rows]
    energies = [r["energy_eV"] for r in rows]
    spins = [r.get("total_spin") for r in rows]
    show_spin = any(s is not None for s in spins)
    e0 = energies[0]
    rel = [e - e0 for e in energies]
    out = out or os.path.join(workdir, "barrier.png")
    if show_spin:
        fig, (ax, ax_spin) = plt.subplots(
            2, 1, figsize=(6.0, 6.4), dpi=160, sharex=True,
            gridspec_kw={"height_ratios": [2.4, 1.0]})
    else:
        fig, ax = plt.subplots(figsize=(6.0, 4.2), dpi=160)
    ax.plot(idx, rel, "-o", color="#1f77b4", lw=2, ms=6,
            markeredgecolor="white", zorder=3)
    imax = int(np.argmax(rel))
    ax.plot([idx[imax]], [rel[imax]], "o", ms=11, mfc="none", mec="#d62728",
            mew=2, zorder=4)
    ax.annotate(f"{rel[imax]:.3f} eV", (idx[imax], rel[imax]),
                textcoords="offset points", xytext=(8, 6), fontsize=9,
                color="#d62728")
    ax.axhline(0.0, color="0.7", lw=0.8, ls="--", zorder=1)
    if not show_spin:
        ax.set_xlabel("image index")
    ax.set_ylabel("energy - E(initial)  [eV]")
    title = (f"VC-NEB  ({summary.get('engine', '?')}, "
             f"{summary.get('n_atoms', '?')} atoms"
             f"{', variable cell' if summary.get('variable_cell') else ''})")
    ax.set_title(title, fontsize=10)
    ax.set_xticks(idx)
    ax.grid(alpha=0.25)
    if show_spin:
        xs = [i for i, s in zip(idx, spins) if s is not None]
        ys = [float(s) for s in spins if s is not None]
        ax_spin.plot(xs, ys, "-o", color="#2ca02c", lw=1.6, ms=4)
        ax_spin.set_ylabel("Mulliken total spin")
        ax_spin.set_xlabel("image index")
        ax_spin.set_xticks(idx)
        ax_spin.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out
