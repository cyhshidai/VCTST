"""One thin wrapper per electronic-structure code.

Every engine exposes the same three things:

    engine.check()                 -> list[str] of problems (empty == fine)
    engine.calculator(atoms, dir)  -> an ASE calculator writing into `dir`
    engine.describe()              -> one-line human summary

QE and ABACUS reuse the ASE calculators; CP2K is driven directly (file based)
because ASE's CP2K interface needs ``cp2k_shell``, which most binary
installations do not ship.
"""

from __future__ import annotations

import copy
import glob
import os
import re
import shlex
import shutil

import numpy as np

from .util import VcnebError, ensure_dir, run_cmd, tail

ENGINE_NAMES = ("qe", "cp2k", "abacus", "emt", "lj")


def _first_token(cmd: str) -> str:
    return shlex.split(str(cmd))[0] if cmd else ""


def _need_file(path: str, what: str, problems: list[str]) -> None:
    if not path:
        problems.append(f"{what} is not configured")
    elif not os.path.isfile(path):
        problems.append(f"{what} not found: {path}")


class BaseEngine:
    name = "base"

    def __init__(self, cfg, workdir: str):
        self.cfg = cfg
        self.calc_cfg = cfg.calc
        self.workdir = workdir
        self.command = self.calc_cfg.get("command")
        self.params = copy.deepcopy(self.calc_cfg.get("params") or {})
        self.timeout = float(self.calc_cfg.get("timeout") or 3600)
        # calc.env: environment variables the code needs (LD_LIBRARY_PATH,
        # OMP_NUM_THREADS, ...).  Subprocesses inherit os.environ, so this
        # reaches pw.x / cp2k / abacus without touching the shell.
        extra_env = self.calc_cfg.get("env") or {}
        if isinstance(extra_env, dict):
            for key, value in extra_env.items():
                os.environ[str(key)] = str(value)

    # override ----------------------------------------------------------
    def check(self) -> list[str]:
        return []

    def warnings(self, symbols=None) -> list[str]:
        """Non fatal advice shown by `vcneb doctor` as WARN."""
        return []

    def calculator(self, atoms, image_dir: str):
        raise NotImplementedError

    def describe(self) -> str:
        return self.name

    # k-points ----------------------------------------------------------
    # ``calc.kpts: {density: X}`` is converted into ONE Monkhorst-Pack grid
    # (spacing <= X, 1/A, i.e. the ASE/CASTEP convention) which is then
    # frozen for the whole run.  Regenerating the grid per image would change
    # the number of k-points in integer steps along the path; each change of
    # the grid puts a small step in the total energy, i.e. a fake force and a
    # fake contribution to the barrier.  A frozen grid keeps the PES smooth
    # (the sampling error then varies continuously along the path); how large
    # that residual error is must be checked with a k-point convergence test.
    def _kspacing(self):
        sp = self.calc_cfg.get("kspacing")
        try:
            sp = float(sp) if sp else None
        except (TypeError, ValueError):
            return None
        return sp if sp and sp > 0 else None

    @staticmethod
    def grid_for_cell(cell, spacing):
        """MP grid with reciprocal spacing <= `spacing` (1/A, no 2*pi)."""
        recip = np.asarray(np.linalg.inv(np.asarray(cell, dtype=float)).T)
        norms = np.linalg.norm(recip, axis=1)
        return tuple(max(1, int(n / spacing) + 1) for n in norms)

    def plan_kpoints(self, cells):
        """Fix the k-mesh once, from every cell that can occur on the path."""
        self._kpts = tuple(int(k) for k in self.calc_cfg["kpts"]) \
            if self.calc_cfg.get("kpts") else None
        spacing = self._kspacing()
        if self._kpts or not spacing:
            return self._kpts
        grids = [self.grid_for_cell(c, spacing) for c in cells if c is not None]
        if not grids:
            return None
        # element-wise maximum: no image is sampled more coarsely than asked
        self._kpts = tuple(max(g[i] for g in grids) for i in range(3))
        return self._kpts

    def kpoints(self):
        """The frozen grid, or None for Gamma-only."""
        if not hasattr(self, "_kpts"):
            self.plan_kpoints([])
        return self._kpts

    def kpoints_too_coarse(self, atoms):
        """Axes where this cell would need more k-points than the frozen grid."""
        grid = self.kpoints()
        spacing = self._kspacing()
        if not grid or not spacing:
            return []
        need = self.grid_for_cell(atoms.get_cell(), spacing)
        return [i for i in range(3) if need[i] > grid[i]]

    def kpts_text(self) -> str:
        if self.calc_cfg.get("gamma"):
            return "Gamma"
        grid, spacing = self.kpoints(), self._kspacing()
        if grid and spacing:
            return f"{list(grid)} (frozen from kspacing {spacing})"
        if grid:
            return str(list(grid))
        if spacing:
            return f"MP grid from kspacing {spacing} (fixed at run start)"
        return "Gamma"


# --------------------------------------------------------------------- EMT/LJ
class EMTEngine(BaseEngine):
    name = "emt"

    def calculator(self, atoms, image_dir: str):
        from ase.calculators.emt import EMT

        return EMT()

    def describe(self) -> str:
        return "EMT (built-in, no input files, for workflow tests only)"


class LJEngine(BaseEngine):
    name = "lj"

    def calculator(self, atoms, image_dir: str):
        from ase.calculators.lj import LennardJones

        return LennardJones(**self.params)

    def describe(self) -> str:
        return "Lennard-Jones (built-in, for workflow tests only)"


# ------------------------------------------------------------------------ QE
class QEEngine(BaseEngine):
    name = "qe"

    def __init__(self, cfg, workdir):
        super().__init__(cfg, workdir)
        self.pseudo_dir = self.calc_cfg.get("pseudo_dir")
        self.pseudos = dict(self.calc_cfg.get("pseudopotentials") or {})

    def check(self) -> list[str]:
        problems: list[str] = []
        exe = _first_token(self.command)
        if not exe:
            problems.append("calc.command is empty, e.g. 'mpirun -np 4 pw.x'")
        elif shutil.which(exe) is None and not os.path.isfile(exe):
            problems.append(f"executable '{exe}' not found in PATH")
        if not self.pseudo_dir:
            problems.append("calc.pseudo_dir is not set (folder with the UPF files)")
        elif not os.path.isdir(self.pseudo_dir):
            problems.append(f"calc.pseudo_dir does not exist: {self.pseudo_dir}")
        if not self.pseudos:
            problems.append(
                "calc.pseudopotentials is empty; add e.g.\n"
                "    pseudopotentials:\n"
                "      Li: Li.upf")
        else:
            for sym, fname in self.pseudos.items():
                cand = os.path.join(self.pseudo_dir or "", fname)
                if self.pseudo_dir and not os.path.isfile(cand):
                    problems.append(
                        f"pseudopotential for {sym} not found: {cand}")
        for must in ("ecutwfc",):
            if must not in (self.params.get("input_data") or {}) and \
                    must not in self.params and "ecutwfc" not in self._flat_params():
                problems.append(f"calc.params.{must} is not set")
        if self.calc_cfg.get("kspacing") and self.calc_cfg.get("kpts"):
            problems.append("set either calc.kpts or kspacing, not both")
        if not self.calc_cfg.get("kpts") and not self.calc_cfg.get("kspacing") \
                and not self.calc_cfg.get("gamma"):
            problems.append(
                "no k-points given: add calc.kpts: [4, 4, 4] "
                "(or {gamma: true} for a molecule)")
        return problems

    def _flat_params(self) -> dict:
        flat = dict(self.params)
        for value in self.params.values():
            if isinstance(value, dict):
                flat.update(value)
        return flat

    def calculator(self, atoms, image_dir: str):
        from ase.calculators.espresso import Espresso, EspressoProfile

        ensure_dir(image_dir)
        profile = EspressoProfile(command=self.command,
                                 pseudo_dir=str(self.pseudo_dir))
        kwargs = copy.deepcopy(self.params)
        input_data = kwargs.pop("input_data", {}) or {}
        control = dict(input_data.get("control", {}) or {})
        control.setdefault("calculation", "scf")
        control.setdefault("tprnfor", True)
        control.setdefault("tstress", True)
        control.setdefault("disk_io", "minimal")
        control.setdefault("verbosity", "high")
        control.setdefault("outdir", os.path.abspath(image_dir))
        control.setdefault("prefix", "vcneb")
        input_data["control"] = control
        kwargs["input_data"] = input_data

        # Frozen integer grid (never `kspacing`: ASE would rebuild it from the
        # current cell for every image, see BaseEngine.plan_kpoints).
        kpts = self.kpoints()
        kwargs["kpts"] = tuple(int(k) for k in kpts) if kpts else None
        kwargs["pseudopotentials"] = self.pseudos
        return Espresso(profile=profile, directory=str(image_dir), **kwargs)

    def describe(self) -> str:
        return (f"Quantum ESPRESSO pw.x  [{self.command}]  "
                f"kpts={self.kpts_text()}")


# ---------------------------------------------------------------------- CP2K
# last (possibly signed) float on a line, Fortran D exponents allowed
_CP2K_FLOAT_RE = re.compile(
    r"([-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?|[-+]?\d+\.\d*(?:[EeDd][-+]?\d+)?)\s*$")


def _parse_cp2k_output(text: str) -> dict:
    from ase.units import Bohr, GPa, Hartree

    energy = None
    for line in text.splitlines():
        if "ENERGY|" in line and "Total FORCE_EVAL" in line:
            # formats seen in the wild:
            #   " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]   -7.4677"
            #   " ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):     -7.4677"
            m = re.search(_CP2K_FLOAT_RE, line)
            if m is None:
                continue
            try:
                energy = float(m.group(1).replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
    if energy is None:
        raise VcnebError(
            "CP2K finished but no 'ENERGY| Total FORCE_EVAL' line was found; "
            "the run probably failed early")

    lines = text.splitlines()
    results = {"energy": energy * Hartree}

    forces = _parse_cp2k_forces(lines)
    if forces is not None and len(forces):
        results["forces"] = forces * Hartree / Bohr

    stress = _parse_cp2k_stress(lines)
    if stress is not None:
        # CP2K prints the opposite sign of the ASE convention
        results["stress"] = -stress
    mulliken = _parse_cp2k_mulliken(lines)
    if mulliken is not None:
        results["mulliken"] = mulliken
    return results


def _cp2k_token_to_float(tok):
    try:
        return float(tok.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def _parse_cp2k_forces(lines):
    """Return (n_atoms, 3) forces in hartree/bohr, or None.

    Handles both the classic block::

         ATOMIC FORCES in [a.u.]
         # Atom   Kind   Element      X        Y        Z
              1      1  Li       0.0      0.0      0.0
         SUM OF ATOMIC FORCES ...

    and the (CP2K >= 2024) block::

         FORCES| Atomic forces [hartree/bohr]
         FORCES|   Atom     x        y        z        |f|
         FORCES|      1  0.0      0.0      0.0   0.0
         FORCES| Sum ...
    """
    start = None
    for idx, line in enumerate(lines):
        low = line.lower()
        if "atomic forces" in low and ("forces|" in low or "atomic forces in" in low):
            start = idx
    if start is None:
        return None
    rows = []
    for line in lines[start + 1:]:
        body = line
        if "FORCES|" in body:
            body = body.split("FORCES|", 1)[1]
        body = body.strip()
        if not body:
            if rows:
                break
            continue
        low = body.lower()
        if low.startswith("sum") or "total atomic force" in low or "stress" in low:
            break
        if re.search(r"[a-z]{2,}", low):        # header / text line
            if rows:
                break
            continue
        parts = body.split()
        if not parts[0].isdigit():
            if rows:
                break
            continue
        if len(parts) >= 6 and _cp2k_token_to_float(parts[2]) is None:
            xyz = parts[3:6]                    # old: idx kind element x y z
        else:
            nums = [p for p in parts[1:] if _cp2k_token_to_float(p) is not None]
            if len(nums) < 3:
                continue
            xyz = nums[:3]                      # new: idx x y z |f|
        vals = [_cp2k_token_to_float(v) for v in xyz]
        if any(v is None for v in vals):
            break
        rows.append(vals)
    return np.array(rows) if rows else None


_CP2K_STRESS_UNITS = {"bar": 1.0e-4, "gpa": 1.0, "kbar": 0.1, "mpa": 1.0e-3,
                      "atm": 1.01325e-4, "pa": 1.0e-9}


def _parse_cp2k_stress(lines):
    """Return the analytic stress tensor in eV/A^3 with ASE's sign convention."""
    from ase.units import GPa

    marker = None
    factor = None
    for idx, line in enumerate(lines):
        low = line.lower()
        if "stress tensor" in low and "eigen" not in low:
            m = re.search(r"\[\s*([a-zA-Z]+)\s*\]|\(\s*([a-zA-Z]+)\s*\)", line)
            unit = None
            if m:
                unit = (m.group(1) or m.group(2)).lower()
            marker = idx
            factor = _CP2K_STRESS_UNITS.get(unit, 1.0)
    if marker is None:
        return None
    rows = []
    for line in lines[marker + 1:]:
        body = line.replace("STRESS|", " ").strip()
        if not body:
            continue
        if body.lower().startswith("analytic") or body.startswith("-"):
            continue
        parts = body.split()
        axis = parts[0].lower()
        if axis not in ("x", "y", "z"):
            if rows:
                break
            continue
        nums = [_cp2k_token_to_float(p) for p in parts[1:]]
        nums = [n for n in nums if n is not None]
        if len(nums) != 3:
            if rows:
                break
            continue
        rows.append(nums)
        if len(rows) == 3:
            break
    if len(rows) != 3:
        return None
    return np.array(rows) * factor * GPa


def _parse_cp2k_mulliken(lines):
    """Last "Mulliken Population Analysis" table of a CP2K log.

    Written by population_analyses.F:write_orbpop for DFT%PRINT%MULLIKEN::

        #  Atom  Element  Kind  Atomic population (alpha,beta) Net charge  Spin moment
              1     Ni       1         9.500000     8.500000     0.000000     1.000000
        # Total charge and spin        9.500000     8.500000     0.000000     1.000000

    without UKS there are only two numbers per atom (population, net charge)
    and the total line is ``# Total charge``.  Only atoms that carry an orbital
    basis are listed, hence the explicit atom index.  With PRINT_GOP the rows
    are per orbital instead of per atom, so the atom table may be missing while
    the total line is still there.
    """
    start = None
    for idx, line in enumerate(lines):
        if "Mulliken Population Analysis" in line:
            start = idx
    if start is None:
        return None

    def floats(tokens):
        out = []
        for tok in tokens:
            val = _cp2k_token_to_float(tok)
            if val is None:
                return None
            out.append(val)
        return out

    atoms = []
    total = None
    for line in lines[start + 1:]:
        body = line.strip()
        if body.startswith("!"):
            break
        if not body:
            if atoms or total is not None:
                break
            continue
        if body.startswith("#"):
            low = body.lower()
            if "total charge" not in low:
                continue                    # column header
            nums = floats([t for t in body.split()
                           if _cp2k_token_to_float(t) is not None])
            if not nums:
                continue
            if len(nums) >= 4:              # UKS: alpha beta charge spin
                total = {"alpha": nums[0], "beta": nums[1],
                         "charge": nums[2], "spin": nums[3]}
            else:                           # RKS: population, net charge
                total = {"alpha": nums[0], "beta": None,
                         "charge": nums[-1], "spin": None}
            continue
        parts = body.split()
        if len(parts) not in (5, 7) or not parts[0].isdigit() \
                or not parts[2].isdigit():
            continue                        # orbital row or junk
        nums = floats(parts[3:])
        if nums is None:
            continue
        row = {"index": int(parts[0]) - 1,   # CP2K is 1 based
               "element": parts[1], "kind": int(parts[2])}
        if len(nums) == 2:                  # RKS
            row.update({"alpha": nums[0], "beta": None,
                        "charge": nums[1], "spin": None})
        else:                               # UKS
            row.update({"alpha": nums[0], "beta": nums[1],
                        "charge": nums[2], "spin": nums[3]})
        atoms.append(row)
    if not atoms and total is None:
        return None
    return {"atoms": atoms, "total": total}


def _cp2k_calculator_class():
    """Build the calculator class lazily so `import vcneb` stays cheap."""
    from ase.calculators.calculator import Calculator, all_changes

    class Cp2kCalculator(Calculator):
        """Minimal file-based ASE calculator: cp2k.popt -i in -o out."""

        implemented_properties = ("energy", "free_energy", "forces", "stress")
        name = "cp2k-vcneb"

        def __init__(self, directory, command, params, timeout=3600,
                     data_dir=None, kpts=None, label="cp2k"):
            super().__init__()
            # NOTE: ASE's ``label`` setter recomputes ``directory``/``prefix``,
            # so it must be set *before* self.directory, otherwise the working
            # directory silently falls back to "." and files land in the CWD.
            self.label = label
            self.directory = str(directory)
            self.prefix = label
            self.command = command
            self.params = dict(params or {})
            self.timeout = timeout
            self.data_dir = data_dir
            self.kpts = list(kpts) if kpts else None
            self._key = None

        # ase API -------------------------------------------------------
        def calculate(self, atoms=None, properties=("energy",),
                      system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            atoms = self.atoms
            if atoms is None:
                raise VcnebError("Cp2kCalculator.calculate called without atoms")
            key = (atoms.get_positions().tobytes(),
                   np.asarray(atoms.get_cell()).tobytes(),
                   tuple(sorted(atoms.get_chemical_symbols())))
            if key != self._key:
                self.results = self._run(atoms)
                self._key = key

        def _run(self, atoms):
            ensure_dir(self.directory)
            inp_name = f"{self.label}.inp"
            out_name = f"{self.label}.out"
            # CP2K appends to an existing -o file, which would mix two runs
            # (and confuse the "last block" logic of the parser): start clean.
            for name in (inp_name, out_name):
                path = os.path.join(self.directory, name)
                if os.path.isfile(path):
                    os.remove(path)
            with open(os.path.join(self.directory, inp_name), "w") as fh:
                fh.write(self._render(atoms))
            cmd = shlex.split(self.command) + ["-i", inp_name, "-o", out_name]
            env = dict(os.environ)
            if self.data_dir:
                env["CP2K_DATA_DIR"] = str(self.data_dir)
            env.setdefault("OMP_NUM_THREADS", "1")
            res = run_cmd(cmd, cwd=self.directory, timeout=self.timeout,
                          env=env,
                          log_path=os.path.join(self.directory, "stdout.log"))
            out_path = os.path.join(self.directory, out_name)
            text = ""
            if os.path.isfile(out_path):
                with open(out_path, errors="replace") as fh:
                    text = fh.read()
            if res["rc"] != 0 or "ABORT" in text[-4000:]:
                raise VcnebError(
                    f"CP2K failed (rc={res['rc']}) in {self.directory}\n"
                    f"---- last lines of {out_name} ----\n{tail(out_path, 25)}")
            parsed = _parse_cp2k_output(text)
            parsed.setdefault("free_energy", parsed["energy"])
            return parsed

        # input rendering -----------------------------------------------
        def _render(self, atoms) -> str:
            return _render_cp2k_input(atoms, self.params, self.kpts)

    return Cp2kCalculator


#: accepted spellings -> the literal CP2K &SMEAR METHOD value
_CP2K_SMEAR_METHODS = {
    "fermi_dirac": "FERMI_DIRAC",
    "fermi-dirac": "FERMI_DIRAC",
    "fermi": "FERMI_DIRAC",
    "fd": "FERMI_DIRAC",
    "energy_window": "ENERGY_WINDOW",
    "energy-window": "ENERGY_WINDOW",
    "ew": "ENERGY_WINDOW",
    "list": "LIST",
    "gauss": "GAUSSIAN",
    "gaussian": "GAUSSIAN",
    "mp": "METHFESSEL_PAXTON",
    "methfessel_paxton": "METHFESSEL_PAXTON",
    "methfessel-paxton": "METHFESSEL_PAXTON",
    "mv": "MARZARI_VANDERBILT",
    "cold": "MARZARI_VANDERBILT",
    "marzari_vanderbilt": "MARZARI_VANDERBILT",
    "marzari-vanderbilt": "MARZARI_VANDERBILT",
}

#: these take SIGMA [eV] instead of ELECTRONIC_TEMPERATURE [K]
_CP2K_SIGMA_METHODS = ("GAUSSIAN", "METHFESSEL_PAXTON", "MARZARI_VANDERBILT")

#: every key _render_cp2k_input actually reads; anything else is ignored by
#: the hand written input writer, so CP2KEngine.warnings() flags it.
_CP2K_PARAM_KEYS = frozenset((
    "basis_set", "potential", "xc", "cutoff", "rel_cutoff", "ngrids",
    "eps_scf", "max_scf", "eps_default", "smearing", "smearing_method",
    "smearing_sigma", "electronic_temperature", "window_size", "added_mos",
    "poisson_solver", "vdw", "extra", "print_level", "basis_set_file",
    "potential_file",
    # spin state, DFT+U and population analysis
    "uks", "multiplicity", "charge", "relax_multiplicity", "plus_u_method",
    "hubbard_u", "hubbard_l", "mulliken", "mulliken_orbitals",
))


def _cp2k_smear_method(value) -> str:
    """Map a user supplied smearing name onto a CP2K &SMEAR METHOD value."""
    key = str(value if value is not None else "fermi_dirac")
    key = key.strip().lower().replace(" ", "_")
    try:
        return _CP2K_SMEAR_METHODS[key]
    except KeyError:
        raise VcnebError(
            f"unknown CP2K smearing_method '{value}'. Use one of: "
            f"fermi_dirac, energy_window, list, gaussian, "
            f"methfessel_paxton, marzari_vanderbilt") from None


_CP2K_PLUS_U_METHODS = {
    "lowdin": "LOWDIN",
    "mulliken": "MULLIKEN",
    "mulliken_charges": "MULLIKEN_CHARGES",
    "mulliken_gop": "MULLIKEN_CHARGES",
    "gop": "MULLIKEN_CHARGES",
}


def _cp2k_plus_u_method(value) -> str:
    """Map a user supplied plus_u_method onto a CP2K PLUS_U_METHOD value."""
    key = str(value if value is not None else "mulliken")
    key = key.strip().lower().replace(" ", "_")
    try:
        return _CP2K_PLUS_U_METHODS[key]
    except KeyError:
        raise VcnebError(
            f"unknown CP2K plus_u_method '{value}'. Use one of: "
            f"lowdin, mulliken, mulliken_charges") from None


def _as_bool(value, default: bool = False) -> bool:
    """Accept YAML booleans as well as the usual true/false spellings."""
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in ("1", "true", "yes", "y", "on", "t"):
        return True
    if text in ("0", "false", "no", "n", "off", "f"):
        return False
    raise VcnebError(f"expected a yes/no value, got '{value}'")


def _cp2k_hubbard_entries(params, symbols=None) -> dict:
    """element -> (U_MINUS_J [eV], L) from ``params.hubbard_u``/``hubbard_l``.

    ``hubbard_u`` is the Dudarev effective U-J of ``&KIND/&DFT_PLUS_U``:

        hubbard_u: {Ni: 4.0}        # per element, eV
        hubbard_l: 2                # 2 = d, 3 = f (default 2)

    A bare number is allowed only when the structure has a single species.
    """
    raw_u = params.get("hubbard_u")
    if raw_u is None or raw_u == {} or raw_u == []:
        return {}
    if isinstance(raw_u, dict):
        table = raw_u
    elif isinstance(raw_u, (int, float)) and not isinstance(raw_u, bool):
        unique = sorted(set(symbols or []))
        if len(unique) != 1:
            raise VcnebError(
                "calc.params.hubbard_u is a single number but the structure "
                f"has {unique or 'unknown'} species; give it per element: "
                "hubbard_u: {Ni: 4.0}")
        table = {unique[0]: raw_u}
    else:
        raise VcnebError(
            "calc.params.hubbard_u must be a mapping element -> U_MINUS_J "
            f"[eV], e.g. {{Ni: 4.0}}, got {raw_u!r}")

    raw_l = params.get("hubbard_l", 2)

    def l_of(sym):
        value = raw_l.get(sym, 2) if isinstance(raw_l, dict) else raw_l
        try:
            value = int(value)
        except (TypeError, ValueError):
            raise VcnebError(
                f"calc.params.hubbard_l for {sym} must be an integer "
                f"(2 = d, 3 = f), got {value!r}") from None
        if not 0 <= value <= 3:
            raise VcnebError(
                f"calc.params.hubbard_l for {sym} must be 0..3 (2 = d, 3 = f)")
        return value

    entries = {}
    for sym, u in table.items():
        sym = str(sym)
        if symbols is not None and sym not in set(symbols):
            raise VcnebError(
                f"calc.params.hubbard_u has an entry for '{sym}' but the "
                f"structure only contains {sorted(set(symbols))}")
        try:
            u_val = float(u)
        except (TypeError, ValueError):
            raise VcnebError(
                f"calc.params.hubbard_u[{sym}] must be a number (eV), "
                f"got {u!r}") from None
        entries[sym] = (u_val, l_of(sym))
    return entries


_CP2K_DATA_NAMES: dict[str, dict] = {}
_ELEMENT_SYMBOL_RE = re.compile(r"^[A-Z][a-z]?$")


def _cp2k_data_names(path: str) -> dict:
    """element -> {name, ...} defined in a CP2K data file (cached).

    In ``GTH_POTENTIALS`` / ``BASIS_MOLOPT`` every entry starts with a header
    line ``Ni GTH-PBE-q18 GTH-PBE`` or
    ``Ni DZVP-MOLOPT-SR-GTH DZVP-MOLOPT-SR-GTH-q18`` followed by numeric data.
    Building the per-element table lets the doctor catch both a misspelt name
    and a name that exists -- but for a different element (e.g. ``GTH-PBE-q10``
    is fine for Pd but not for Ni, which only has q18).
    """
    if path in _CP2K_DATA_NAMES:
        return _CP2K_DATA_NAMES[path]
    table: dict = {}
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                body = line.strip()
                if not body or body.startswith("#"):
                    continue
                parts = body.split()
                if len(parts) < 2 or not _ELEMENT_SYMBOL_RE.match(parts[0]):
                    continue
                names = parts[1:]
                if any(_cp2k_token_to_float(tok) is not None
                       for tok in names):
                    continue                    # numeric data line
                table.setdefault(parts[0], set()).update(names)
    except OSError:
        table = {}
    _CP2K_DATA_NAMES[path] = table
    return table


def _render_cp2k_input(atoms, params, kpts) -> str:
    p = dict(params or {})
    cell = atoms.get_cell()
    a, b, c, alpha, beta, gamma = cell.cellpar()
    basis = p.get("basis_set", "DZVP-MOLOPT-SR-GTH")
    potential = p.get("potential", "GTH-PBE-q3")
    xc = p.get("xc", "PBE")
    cutoff = p.get("cutoff", 400)
    rel_cutoff = p.get("rel_cutoff", 50)
    n_grids = p.get("ngrids", 5)
    eps_scf = p.get("eps_scf", 1.0e-7)
    max_scf = p.get("max_scf", 50)
    eps_default = p.get("eps_default", 1.0e-12)
    smear = bool(p.get("smearing", False))
    e_temp = p.get("electronic_temperature", 1000)
    smear_method = _cp2k_smear_method(p.get("smearing_method"))
    smear_sigma = p.get("smearing_sigma", None)   # eV, keyword SIGMA
    window_size = p.get("window_size", None)      # eV, keyword WINDOW_SIZE
    added_mos = p.get("added_mos", None)
    poisson = p.get("poisson_solver", None)
    extra = p.get("extra", "") or ""
    vdw = p.get("vdw", None)
    symbols = atoms.get_chemical_symbols()
    # spin state / occupancy / DFT+U / population analysis
    uks = _as_bool(p.get("uks"), False)
    multiplicity = p.get("multiplicity", None)
    if multiplicity in (None, "", 0, "0"):
        multiplicity = None
    charge = int(p.get("charge", 0) or 0)
    relax_mult = p.get("relax_multiplicity", None)
    plus_u_method = p.get("plus_u_method", None)
    mulliken = _as_bool(p.get("mulliken"), True)
    mulliken_gop = _as_bool(p.get("mulliken_orbitals"), False)
    hubbard = _cp2k_hubbard_entries(p, symbols)

    def per_element(value, symbol):
        if isinstance(value, dict):
            if symbol not in value:
                raise VcnebError(
                    f"no entry for element {symbol} in the CP2K parameter "
                    f"{value}; add it as a mapping element->name")
            return value[symbol]
        return value

    out = []
    w = out.append
    w("&GLOBAL")
    w("  PROJECT vcneb")
    w("  RUN_TYPE ENERGY_FORCE")
    # PRINT_LEVEL must be >= MEDIUM, otherwise CP2K does not print the
    # "ATOMIC FORCES in [a.u.]" block nor the stress tensor and the parser
    # ends up with energy only.
    w(f"  PRINT_LEVEL {p.get('print_level', 'MEDIUM')}")
    w("&END GLOBAL")
    w("&FORCE_EVAL")
    w("  STRESS_TENSOR ANALYTICAL")  # needed by the variable-cell part
    # CP2K >= 2024 only prints the force/stress block when the run is asked
    # for it explicitly (PRINT_LEVEL alone is not enough at MEDIUM).
    w("  &PRINT")
    w("    &FORCES ON")
    w("    &END FORCES")
    w("    &STRESS_TENSOR ON")
    w("    &END STRESS_TENSOR")
    w("  &END PRINT")
    w("  METHOD QS")
    w("  &DFT")
    # CP2K keywords are read from the &DFT section, so they must come before
    # the first subsection.
    if uks:
        # unrestricted (spin polarised) Kohn-Sham; needed for any magnetic
        # state and for a meaningful spin moment in the population analysis
        w("    UKS T")
    if multiplicity:
        # 2S+1; CP2K's own default is 1 for an even, 2 for an odd electron count
        w(f"    MULTIPLICITY {int(multiplicity)}")
    if charge:
        w(f"    CHARGE {charge}")
    if relax_mult:
        w(f"    RELAX_MULTIPLICITY {float(relax_mult)}")
    if plus_u_method or hubbard:
        # MULLIKEN is the CP2K default (cheap, needs no S^{1/2}); the +U term
        # is built from this very population analysis, so keep them together.
        w(f"    PLUS_U_METHOD {_cp2k_plus_u_method(plus_u_method or 'mulliken')}")
    w(f"    BASIS_SET_FILE_NAME {p.get('basis_set_file', 'BASIS_MOLOPT')}")
    w(f"    POTENTIAL_FILE_NAME {p.get('potential_file', 'GTH_POTENTIALS')}")
    w("    &MGRID")
    w(f"      CUTOFF {cutoff}")
    w(f"      NGRIDS {n_grids}")
    w(f"      REL_CUTOFF {rel_cutoff}")
    w("    &END MGRID")
    w("    &QS")
    w(f"      EPS_DEFAULT {eps_default}")
    w("    &END QS")
    if kpts:
        w("    &KPOINTS")
        w("      SCHEME MONKHORST-PACK "
          + " ".join(str(int(k)) for k in kpts))
        w("    &END KPOINTS")
    if poisson:
        w("    &POISSON")
        w(f"      POISSON_SOLVER {poisson}")
        w("    &END POISSON")
    w("    &SCF")
    w("      SCF_GUESS ATOMIC")
    w(f"      EPS_SCF {eps_scf}")
    w(f"      MAX_SCF {max_scf}")
    if smear:
        w("      &SMEAR ON")
        # CP2K >= 2024 accepts FERMI_DIRAC, ENERGY_WINDOW, LIST, GAUSSIAN,
        # METHFESSEL_PAXTON and MARZARI_VANDERBILT (see input_cp2k_scf.F).
        # The last three take a width in energy units (SIGMA), not a temperature.
        w(f"        METHOD {smear_method}")
        if smear_method in _CP2K_SIGMA_METHODS:
            if smear_sigma is not None:
                w(f"        SIGMA [eV] {float(smear_sigma)}")
        elif smear_method == "FERMI_DIRAC":
            w(f"        ELECTRONIC_TEMPERATURE [K] {e_temp}")
        elif smear_method == "ENERGY_WINDOW" and window_size is not None:
            w(f"        WINDOW_SIZE [eV] {float(window_size)}")
        w("      &END SMEAR")
        if added_mos:
            w(f"      ADDED_MOS {added_mos}")
        w("      &DIAGONALIZATION")
        w("        ALGORITHM STANDARD")
        w("      &END DIAGONALIZATION")
        w("      &MIXING")
        w("        METHOD BROYDEN_MIXING")
        w("        ALPHA 0.2")
        w("        NBROYDEN 8")
        w("      &END MIXING")
    else:
        w("      &OT")
        w("        MINIMIZER DIIS")
        w("        PRECONDITIONER FULL_SINGLE_INVERSE")
        w("      &END OT")
        w("      &OUTER_SCF")
        w(f"        EPS_SCF {eps_scf}")
        w("        MAX_SCF 20")
        w("      &END OUTER_SCF")
    w("    &END SCF")
    if mulliken:
        # $DFT/&PRINT/&MULLIKEN -> "Mulliken Population Analysis" table in the
        # log (per atom: population alpha/beta, net charge, spin moment).
        # Cheap: no diagonalisation of S, unlike PLUS_U_METHOD LOWDIN.
        w("    &PRINT")
        w("      &MULLIKEN ON")
        if mulliken_gop:
            w("        PRINT_GOP ON")   # per-orbital (gross orbital) breakdown
        w("      &END MULLIKEN")
        w("    &END PRINT")
    w("    &XC")
    w(f"      &XC_FUNCTIONAL {xc}")
    w("      &END XC_FUNCTIONAL")
    if vdw:
        w("      &VDW_POTENTIAL")
        w(f"        POTENTIAL_TYPE {vdw}")
        w("      &END VDW_POTENTIAL")
    w("    &END XC")
    if extra.strip():
        for line in extra.strip("\n").splitlines():
            w("    " + line.rstrip())
    w("  &END DFT")
    w("  &SUBSYS")
    w("    &CELL")
    w(f"      ABC {a:.10f} {b:.10f} {c:.10f}")
    w(f"      ALPHA_BETA_GAMMA {alpha:.6f} {beta:.6f} {gamma:.6f}")
    w("      PERIODIC XYZ")
    w("    &END CELL")
    w("    &COORD")
    for sym, pos in zip(symbols, atoms.get_positions()):
        w(f"      {sym} {pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}")
    w("    &END COORD")
    seen = []
    for sym in symbols:
        if sym in seen:
            continue
        seen.append(sym)
        w(f"    &KIND {sym}")
        w(f"      ELEMENT {sym}")
        w(f"      BASIS_SET {per_element(basis, sym)}")
        w(f"      POTENTIAL {per_element(potential, sym)}")
        if sym in hubbard:
            u_val, l_val = hubbard[sym]
            # Dudarev DFT+U in the QS (Gaussian basis) code path: only
            # U_MINUS_J and L are used; U, J, alpha, beta, N belong to pwdft.
            w("      &DFT_PLUS_U ON")
            w(f"        L {l_val}")
            w(f"        U_MINUS_J [eV] {float(u_val)}")
            w("      &END DFT_PLUS_U")
        w("    &END KIND")
    w("  &END SUBSYS")
    w("&END FORCE_EVAL")
    return "\n".join(out) + "\n"


class CP2KEngine(BaseEngine):
    name = "cp2k"

    def __init__(self, cfg, workdir):
        super().__init__(cfg, workdir)
        self.data_dir = self.calc_cfg.get("data_dir")

    def check(self) -> list[str]:
        problems: list[str] = []
        exe = _first_token(self.command)
        if not exe:
            problems.append("calc.command is empty, e.g. 'mpirun -np 4 cp2k.popt'")
        elif shutil.which(exe) is None and not os.path.isfile(exe):
            problems.append(f"executable '{exe}' not found in PATH")
        else:
            for token in shlex.split(self.command):
                if token.startswith("cp2k") and shutil.which(token) is None:
                    problems.append(
                        f"'{token}' not found in PATH (this installation only "
                        f"has cp2k.popt / cp2k.psmp)")
        datadir = self.data_dir or os.environ.get("CP2K_DATA_DIR")
        if not datadir:
            problems.append(
                "calc.data_dir is not set; point it at the CP2K data folder "
                "(the one containing BASIS_MOLOPT and GTH_POTENTIALS)")
        elif not os.path.isdir(datadir):
            problems.append(f"calc.data_dir does not exist: {datadir}")
        else:
            bs_file = self.params.get("basis_set_file", "BASIS_MOLOPT")
            pot_file = self.params.get("potential_file", "GTH_POTENTIALS")
            for fname in (bs_file, pot_file):
                if not os.path.isfile(os.path.join(datadir, fname)):
                    problems.append(f"{fname} not found in {datadir}")
        if isinstance(self.params.get("basis_set"), dict) or \
                isinstance(self.params.get("potential"), dict):
            pass  # per-element tables are validated at run time
        try:
            _cp2k_smear_method(self.params.get("smearing_method"))
        except VcnebError as exc:
            problems.append(str(exc))
        return problems

    def warnings(self, symbols=None) -> list[str]:
        """Non fatal advice; `vcneb doctor` prints these as WARN."""
        notes: list[str] = []
        if not self.kpoints() and not self._kspacing():
            notes.append(
                "no k-points configured: CP2K will run Gamma-point only. "
                "Add calc.kpts: [4, 4, 4] (or {density: 0.03}) for a periodic "
                "system")
        unknown = sorted(set(self.params) - _CP2K_PARAM_KEYS)
        if unknown:
            notes.append(
                "calc.params keys that the CP2K writer does not render and "
                "will silently ignore: " + ", ".join(unknown) +
                "   (use params.extra to inject raw lines into &DFT)")
        notes.extend(self._warnings_spin(symbols))
        notes.extend(self._warnings_kind(symbols))
        try:
            method = _cp2k_smear_method(self.params.get("smearing_method"))
        except VcnebError:
            return notes
        if not self.params.get("smearing"):
            if self.params.get("smearing_method") is not None:
                notes.append(
                    "calc.params.smearing_method is set but smearing is "
                    "false; no &SMEAR block will be written")
            return notes
        if method in _CP2K_SIGMA_METHODS:
            if self.params.get("smearing_sigma") is None:
                notes.append(
                    f"{method} smearing uses SIGMA, not "
                    f"electronic_temperature; add e.g. smearing_sigma: 0.05 "
                    f"(eV), otherwise the CP2K default (0.002 Ha = 0.054 eV) "
                    f"is used")
            if self.params.get("electronic_temperature") is not None:
                notes.append(
                    f"electronic_temperature is ignored by {method} smearing; "
                    f"it only applies to fermi_dirac")
        elif method == "FERMI_DIRAC" and \
                self.params.get("smearing_sigma") is not None:
            notes.append(
                "smearing_sigma is ignored by fermi_dirac smearing; it only "
                "applies to gaussian / methfessel_paxton / "
                "marzari_vanderbilt")
        elif method == "ENERGY_WINDOW" and \
                self.params.get("window_size") is None:
            notes.append(
                "energy_window smearing needs window_size [eV], e.g. 0.3")
        return notes

    def _warnings_kind(self, symbols=None) -> list[str]:
        """Basis set / pseudopotential names per element."""
        notes: list[str] = []
        unique = sorted(set(symbols or []))
        if not unique:
            return notes
        basis = self.params.get("basis_set", "DZVP-MOLOPT-SR-GTH")
        potential = self.params.get("potential", "GTH-PBE-q3")
        if len(unique) > 1:
            for key, value in (("basis_set", basis), ("potential", potential)):
                if not isinstance(value, dict):
                    notes.append(
                        f"calc.params.{key}='{value}' is used for every "
                        f"species, but the structure has {unique}; give it "
                        f"per element, e.g. {key}: {{{unique[0]}: ...}}")
        datadir = self.data_dir or os.environ.get("CP2K_DATA_DIR")
        if not datadir or not os.path.isdir(datadir):
            return notes
        for key, value, fname in (
                ("basis_set", basis, self.params.get("basis_set_file",
                                                     "BASIS_MOLOPT")),
                ("potential", potential, self.params.get("potential_file",
                                                         "GTH_POTENTIALS"))):
            table = _cp2k_data_names(os.path.join(datadir, str(fname)))
            if not table:
                continue
            for sym in unique:
                name = value.get(sym) if isinstance(value, dict) else value
                if not name:
                    continue
                available = table.get(sym)
                if available is None:
                    notes.append(
                        f"{fname} has no entry for {sym}; check the species "
                        f"name or use another {key}_file")
                elif str(name) not in available:
                    hint = ""
                    if key == "potential":
                        same = sorted(n for n in available if "PBE" in n)
                        hint = ("; available PBE potentials for " + sym + ": "
                                + ", ".join(same)) if same else ""
                    notes.append(
                        f"calc.params.{key} '{name}' does not exist for {sym} "
                        f"in {fname}{hint}")
        return notes

    def _warnings_spin(self, symbols=None) -> list[str]:
        """Advice about UKS / multiplicity / DFT+U / Mulliken output."""
        notes: list[str] = []
        uks = _as_bool(self.params.get("uks"), False)
        multiplicity = self.params.get("multiplicity")
        if multiplicity and not uks:
            notes.append(
                "calc.params.multiplicity is ignored without uks: true "
                "(CP2K then runs a closed-shell calculation); set "
                "uks: true for a spin polarised / magnetic system")
        relax = self.params.get("relax_multiplicity")
        if relax:
            if not uks:
                notes.append(
                    "RELAX_MULTIPLICITY requires uks: true, CP2K aborts "
                    "otherwise")
            if not self.params.get("added_mos"):
                notes.append(
                    "RELAX_MULTIPLICITY needs added_mos > 0 (empty bands); "
                    "without them it has no effect")
        try:
            hubbard = _cp2k_hubbard_entries(self.params, symbols)
        except VcnebError as exc:
            notes.append(str(exc))
            hubbard = {}
        if self.params.get("plus_u_method") is not None:
            try:
                _cp2k_plus_u_method(self.params.get("plus_u_method"))
            except VcnebError as exc:
                notes.append(str(exc))
        if hubbard:
            spec = ", ".join(f"{sym}: U-J={u:g} eV, L={l}"
                             for sym, (u, l) in sorted(hubbard.items()))
            if not uks:
                notes.append(
                    "DFT+U is applied to a closed-shell density (uks is not "
                    f"set); {spec} still shifts the energy, but there is no "
                    "spin moment and FM/AFM cannot be distinguished")
            if not self.params.get("smearing"):
                notes.append(
                    "DFT+U without smearing (OT/diagonalisation) can converge "
                    "to a metastable occupation; for magnetic metals set "
                    "smearing: true and check the Mulliken spin moments")
            notes.append(
                "DFT+U (" + spec + ") makes the SCF much stiffer: converge "
                "cutoff/rel_cutoff again and expect 'DFT+U energy "
                "contribution is negative' warnings -- those are normal, not "
                "errors")
        if uks and self.params.get("multiplicity"):
            if self.params.get("smearing"):
                # measured: MULTIPLICITY 3 + &SMEAR on Ni gave M=0.81, not 2
                notes.append(
                    "with &SMEAR the MULTIPLICITY only fixes the initial "
                    "guess: alpha/beta occupations follow one Fermi level, so "
                    "the final spin moment floats. Read it from the Mulliken "
                    "table / barrier.csv instead of trusting the input")
            else:
                # qs_environment.F: nelectron_spin = (N +- (mult-1))/2, i.e.
                # each spin channel is filled with a fixed electron count
                notes.append(
                    "without smearing MULTIPLICITY pins the alpha/beta "
                    "electron counts (CP2K fills the two spin channels "
                    "separately), which biases the SCF towards that magnetic "
                    "state. If the state may change along the path, use "
                    "smearing or relax_multiplicity: <tol>")
        if not _as_bool(self.params.get("mulliken"), True):
            notes.append(
                "mulliken: false -> no charge/spin table in the log and no "
                "spin moments in barrier.csv")
        elif not uks:
            notes.append(
                "Mulliken output is on but the run is closed shell (uks not "
                "set): the table shows population and net charge only, no "
                "spin moment. Add uks: true if the system is magnetic")
        return notes

    def calculator(self, atoms, image_dir: str):
        ensure_dir(image_dir)
        cls = _cp2k_calculator_class()
        return cls(directory=image_dir, command=self.command,
                   params=self.params, timeout=self.timeout,
                   data_dir=self.data_dir, kpts=self.kpoints())

    def describe(self) -> str:
        spin = "UKS" if _as_bool(self.params.get("uks"), False) else "RKS"
        if self.params.get("multiplicity"):
            spin += f" mult={int(self.params['multiplicity'])}"
        try:
            hubbard = _cp2k_hubbard_entries(self.params)
        except VcnebError:
            hubbard = {}
        if hubbard:
            spin += "  DFT+U(" + ",".join(
                f"{sym} U-J={u:g}eV L={l}"
                for sym, (u, l) in sorted(hubbard.items())) + ")"
        if _as_bool(self.params.get("mulliken"), True):
            spin += "  mulliken=on"
        return (f"CP2K  [{self.command}]  "
                f"cutoff={self.params.get('cutoff', 400)} Ry  "
                f"basis={self.params.get('basis_set', 'DZVP-MOLOPT-SR-GTH')}  "
                f"kpts={self.kpts_text()}  {spin}")


# -------------------------------------------------------------------- ABACUS
class AbacusEngine(BaseEngine):
    name = "abacus"

    def __init__(self, cfg, workdir):
        super().__init__(cfg, workdir)
        self.pseudo_dir = self.calc_cfg.get("pseudo_dir")
        self.orb_dir = self.calc_cfg.get("data_dir")
        self.pp = dict(self.calc_cfg.get("pseudopotentials") or {})
        self.orb = dict(self.calc_cfg.get("orbitals") or {})

    def check(self) -> list[str]:
        problems: list[str] = []
        exe = _first_token(self.command)
        if not exe:
            problems.append("calc.command is empty, e.g. 'mpirun -np 4 abacus'")
        elif shutil.which(exe) is None and not os.path.isfile(exe):
            problems.append(f"executable '{exe}' not found in PATH")
        if not self.pp:
            problems.append(
                "calc.pseudopotentials is empty; add e.g.\n"
                "    pseudopotentials:\n"
                "      Cu: Cu_ONCV_PBE-1.0.upf")
        for sym, fname in self.pp.items():
            if self.pseudo_dir and not os.path.isfile(
                    os.path.join(self.pseudo_dir, fname)):
                problems.append(
                    f"UPF for {sym} not found: {os.path.join(self.pseudo_dir, fname)}")
        if self.params.get("basis_type", "lcao") == "lcao":
            if not self.orb:
                problems.append(
                    "calc.orbitals is empty; LCAO needs orbital files, e.g.\n"
                    "    orbitals:\n"
                    "      Cu: Cu_gga_9au_100Ry_4s2p2d1f.orb")
            for sym, fname in self.orb.items():
                if self.orb_dir and not os.path.isfile(
                        os.path.join(self.orb_dir, fname)):
                    problems.append(
                        f"orbital file for {sym} not found: "
                        f"{os.path.join(self.orb_dir, fname)}")
        if "ecutwfc" not in self.params:
            problems.append("calc.params.ecutwfc is not set (Ry)")
        if self.params.get("ks_solver") == "cusolver":
            problems.append(
                "calc.params.ks_solver='cusolver' needs a GPU; "
                "use 'genelpa' or 'scalapack_gvx' on CPU")
        if self.params.get("device") == "gpu" and \
                self.params.get("ks_solver") != "cusolver":
            problems.append("device: gpu is set but ks_solver is not cusolver")
        return problems

    def calculator(self, atoms, image_dir: str):
        from ase.calculators.abacus import Abacus, AbacusProfile

        ensure_dir(image_dir)
        lcao = (self.params.get("basis_type", "lcao") or "lcao").lower() == "lcao"
        # ASE's Abacus calculator falls back to $ABACUS_ORBITAL_PATH (KeyError!)
        # whenever no orbital information is given, which happens for pw basis.
        # It only writes `orbital_dir` into INPUT, which ABACUS ignores for
        # basis_type=pw, so any existing directory is fine there.
        profile = AbacusProfile(command=self.command,
                                pseudo_dir=self.pseudo_dir,
                                basis_dir=self.orb_dir if lcao
                                else (self.orb_dir or self.pseudo_dir))
        kwargs = copy.deepcopy(self.params)
        kwargs.setdefault("basis_type", "lcao")
        kwargs.setdefault("ks_solver", "genelpa")
        kwargs.setdefault("symmetry", 0)
        kwargs.setdefault("scf_thr", 1.0e-7)
        kwargs.setdefault("scf_nmax", 100)
        kwargs.setdefault("cal_force", 1)
        kwargs.setdefault("cal_stress", 1)
        kwargs.setdefault("device", "cpu")
        if self.pp:
            kwargs["pp"] = dict(self.pp)
        if self.orb:
            kwargs["basis"] = dict(self.orb)
        if self.calc_cfg.get("gamma"):
            kwargs["gamma_only"] = 1
        elif self.kpoints():
            # integer KPT file, same grid for every image (frozen, see
            # BaseEngine.plan_kpoints).  Do not hand `kspacing` to ABACUS:
            # its unit is 1/Bohr, ASE's is 1/Angstrom (factor 3.32 apart).
            kwargs["kpts"] = [int(k) for k in self.kpoints()]
        calc = Abacus(profile=profile, directory=str(image_dir), **kwargs)
        return calc

    def describe(self) -> str:
        return (f"ABACUS  [{self.command}]  ecutwfc="
                f"{self.params.get('ecutwfc')} Ry  kpts={self.kpts_text()}")


_REGISTRY = {
    "qe": QEEngine,
    "cp2k": CP2KEngine,
    "abacus": AbacusEngine,
    "emt": EMTEngine,
    "lj": LJEngine,
}


def make_engine(cfg, workdir: str | None = None) -> BaseEngine:
    name = (cfg.calc.engine or "").lower()
    if name not in _REGISTRY:
        raise VcnebError(
            f"unknown engine '{name}'; choose one of {', '.join(ENGINE_NAMES)}")
    return _REGISTRY[name](cfg, workdir or cfg.workdir)
