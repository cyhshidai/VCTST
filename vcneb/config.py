"""Configuration loading / validation for vcneb.

The user only ever edits one YAML file.  Everything the driver needs
(structure, engine, DFT parameters, NEB parameters) lives there.
"""

from __future__ import annotations

import copy
import os
import shutil

from .util import VcnebError, ensure_dir

try:  # PyYAML is the only non-core dependency besides numpy/ase
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


ENGINE_NAMES = ("qe", "cp2k", "abacus", "emt", "lj")


DEFAULTS = {
    "project": {
        "name": "vcneb",
        "workdir": "./run",
    },
    "structure": {
        "initial": None,
        "final": None,
        "build_initial": None,
        "build_final": None,
        "vacuum": None,
    },
    "calc": {
        "engine": None,
        "command": None,
        "mpi": 1,
        "env": {},
        "pseudo_dir": None,
        "data_dir": None,
        "timeout": 3600,
        "params": {},
        "kpts": None,
        "gamma": None,
        "kspacing": None,
        "pseudopotentials": {},
        "orbitals": {},
    },
    "neb": {
        "images": 5,
        "interpolation": "idpp",
        "align_cell": True,
        "tangent": "improved",  # improved (Henkelman-Jonsson) or central
        "min_distance": 2.0,
        "prune": "auto",
        "variable_cell": True,
        "pressure": 0.0,  # GPa; the band works on the enthalpy H = E + P*V
        "spring_k": 0.1,
        "f_th": 0.03,
        "max_steps": 300,
        "dt": 0.1,
        "dt_max": 0.3,
        "max_step": 0.2,
        "cell_max_step": 0.02,
        "cell_weight": None,  # null = derive from max_step/cell_max_step
        "climbing_image": True,
        "relax_endpoints": True,
        "endpoint_fmax": 0.02,
        "endpoint_steps": 200,
        "parallel_images": 1,
        "keep_intermediate_images": True,
    },
    "output": {
        "plot": True,
        "save_trajectory": True,
    },
}


def _deep_merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for key, val in (over or {}).items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], val)
        else:
            out[key] = copy.deepcopy(val)
    return out


class Section(dict):
    """dict with attribute access, so cfg.neb.f_th works and is typo-checked."""

    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError as exc:
            raise AttributeError(
                f"unknown config key '{item}' in section "
                f"{sorted(self)}") from exc

    def __setattr__(self, key, value):
        self[key] = value


def _wrap(obj):
    if isinstance(obj, dict):
        return Section({k: _wrap(v) for k, v in obj.items()})
    if isinstance(obj, list):
        return [_wrap(v) for v in obj]
    return obj


class Config:
    """Resolved configuration.

    Paths inside the YAML are interpreted relative to the YAML file, which is
    what people expect when they move a folder around.
    """

    def __init__(self, data: dict, source: str | None = None):
        self.source = os.path.abspath(source) if source else None
        self.base_dir = os.path.dirname(self.source) if source else os.getcwd()
        merged = _deep_merge(DEFAULTS, data or {})
        for key in ("project", "structure", "calc", "neb", "output"):
            if not isinstance(merged.get(key), dict):
                raise VcnebError(f"section '{key}' in {source} must be a mapping")
        self.data = merged
        self.project = _wrap(merged["project"])
        self.structure = _wrap(merged["structure"])
        self.calc = _wrap(merged["calc"])
        self.neb = _wrap(merged["neb"])
        self.output = _wrap(merged["output"])
        self._resolve_paths()

    # ------------------------------------------------------------------ paths
    def abspath(self, path: str) -> str:
        if path is None:
            return None
        path = os.path.expanduser(str(path))
        if os.path.isabs(path):
            return path
        return os.path.normpath(os.path.join(self.base_dir, path))

    def _resolve_paths(self) -> None:
        self.workdir = self.abspath(self.project.workdir or "./run")
        for key in ("initial", "final"):
            val = self.structure.get(key)
            if val:
                self.structure[key] = self.abspath(val)
        for key in ("pseudo_dir", "data_dir"):
            val = self.calc.get(key)
            if val:
                self.calc[key] = self.abspath(val)
        self.logs_dir = os.path.join(self.workdir, "logs")
        self.images_dir = os.path.join(self.workdir, "images")

    # -------------------------------------------------------------- validate
    def validate(self, for_run: bool = True) -> list[str]:
        problems: list[str] = []
        engine = (self.calc.engine or "").lower()
        if not engine:
            problems.append(
                "calc.engine is missing; pick one of " + ", ".join(ENGINE_NAMES))
        elif engine not in ENGINE_NAMES:
            problems.append(
                f"calc.engine='{engine}' is not supported; use one of "
                + ", ".join(ENGINE_NAMES))
        self.calc["engine"] = engine

        init, final = self.structure.get("initial"), self.structure.get("final")
        if not init and not self.structure.get("build_initial"):
            problems.append("structure.initial (a file ASE can read) or "
                            "structure.build_initial is required")
        if not final and not self.structure.get("build_final"):
            problems.append("structure.final (a file ASE can read) or "
                            "structure.build_final is required")
        for label, path in (("initial", init), ("final", final)):
            if path and not os.path.isfile(path):
                problems.append(f"structure.{label} file not found: {path}")

        n_img = int(self.neb.images)
        if n_img < 1:
            problems.append("neb.images must be >= 1 (it counts movable images)")
        if float(self.neb.f_th) <= 0:
            problems.append("neb.f_th must be > 0 (eV/Angstrom)")
        if int(self.neb.max_steps) < 1:
            problems.append("neb.max_steps must be >= 1")
        if str(self.neb.interpolation).lower() not in ("linear", "idpp"):
            problems.append("neb.interpolation must be 'linear' or 'idpp'")
        tangent = str(self.neb.get("tangent") or "improved").strip().lower()
        if tangent not in ("improved", "central"):
            problems.append("neb.tangent must be 'improved' "
                            "(Henkelman-Jonsson) or 'central'")
        self.neb["tangent"] = tangent
        try:
            dmin = float(self.neb.min_distance or 0.0)
        except (TypeError, ValueError):
            problems.append("neb.min_distance must be a number in Angstrom "
                            "(0 disables the check)")
            dmin = 0.0
        if dmin < 0:
            problems.append("neb.min_distance must be >= 0")
        prune = self.neb.get("prune")
        if isinstance(prune, str):
            prune = prune.strip().lower()
            if prune in ("false", "no", "off", "none"):
                prune = False
            elif prune in ("true", "yes", "on"):
                prune = "auto"
            self.neb["prune"] = prune
        elif prune is True:
            prune = "auto"
            self.neb["prune"] = prune
        if prune not in ("auto", "error", False):
            problems.append("neb.prune must be 'auto' (drop the images), "
                            "'error' (stop and tell me) or false (keep them)")

        try:
            press = float(self.neb.get("pressure") or 0.0)
        except (TypeError, ValueError):
            problems.append("neb.pressure must be a number in GPa "
                            "(0 = no pressure)")
            press = 0.0
        self.neb["pressure"] = press

        if engine in ("qe", "cp2k", "abacus"):
            if not self.calc.command:
                problems.append(
                    f"calc.command is required for engine='{engine}', "
                    f"e.g. 'mpirun -np 4 pw.x'")
            elif not for_run:
                exe = str(self.calc.command).split()[0]
                if shutil.which(exe) is None:
                    problems.append(f"executable '{exe}' not found in PATH")
        self._validate_kpts(problems)
        return problems

    def _validate_kpts(self, problems: list[str]) -> None:
        kpts = self.calc.get("kpts")
        if isinstance(kpts, dict):
            if "density" in kpts:
                self.calc["kspacing"] = float(kpts["density"])
                self.calc["kpts"] = None
            elif kpts.get("gamma"):
                self.calc["gamma"] = True
                self.calc["kpts"] = None
            elif "size" in kpts:
                self.calc["kpts"] = list(kpts["size"])
            else:
                problems.append(
                    "calc.kpts must be a list like [8, 8, 8], or "
                    "{gamma: true}, or {density: 0.03}")
        elif isinstance(kpts, (list, tuple)):
            if len(kpts) != 3 or any(int(k) < 1 for k in kpts):
                problems.append("calc.kpts must be three positive integers")
            else:
                self.calc["kpts"] = [int(k) for k in kpts]
        if self.calc.get("gamma"):
            self.calc["gamma"] = True

    # --------------------------------------------------------------- helpers
    def summary_lines(self) -> list[str]:
        c, n = self.calc, self.neb
        lines = [
            f"project      : {self.project.name}",
            f"workdir      : {self.workdir}",
            f"engine       : {c.engine}   command: {c.command or '-'}",
            f"images       : {n.images} movable "
            f"(+2 endpoints), interpolation={n.interpolation}",
            f"variable cell: {n.variable_cell}   spring_k={n.spring_k}   "
            f"f_th={n.f_th} eV/A",
            f"max_steps    : {n.max_steps}   dt={n.dt}   "
            f"climbing_image={n.climbing_image}",
        ]
        if float(n.get("pressure") or 0.0):
            lines.append(
                f"pressure     : {float(n.pressure)} GPa "
                f"(every energy below is the enthalpy H = E + P*V)")
        return lines


def load_config(path: str) -> Config:
    if yaml is None:
        raise VcnebError(
            "PyYAML is not installed.  Run:  pip install pyyaml")
    if not os.path.isfile(path):
        raise VcnebError(f"config file not found: {path}")
    with open(path, "r", encoding="utf-8") as fh:
        try:
            data = yaml.safe_load(fh) or {}
        except yaml.YAMLError as exc:
            raise VcnebError(f"cannot parse {path}:\n{exc}") from exc
    if not isinstance(data, dict):
        raise VcnebError(f"{path} must contain a YAML mapping at top level")
    return Config(data, source=path)


# ------------------------------------------------------------------ templates
_HEADER = """\
# vcneb configuration  (generated by `vcneb init`)
# Run `vcneb doctor {fname}` after editing, then `vcneb run {fname}`.
# Lengths in Angstrom, energies in eV, forces in eV/Angstrom.
#
# 只需要改三处：结构（structure）、算例命令（calc.command / 赝势路径）、
# 以及 NEB 参数（neb.images / neb.f_th / neb.max_steps）。
# 注意：截断能、k 点、赝势必须自己先做收敛测试，默认值只保证能跑通。
"""

_TEMPLATES = {
    "qe": """\
project:
  name: li_bcc_to_hcp
  workdir: ./run_li

# Any format ASE can read: POSCAR, .cif, .xyz, .traj, ...
structure:
  initial: POSCAR_initial
  final: POSCAR_final

calc:
  engine: qe
  command: mpirun -np 4 pw.x       # 4 MPI ranks per image
  mpi: 4
  # environment variables the code needs (optional), e.g. a library path
  # that is not in the default loader path:
  # env:
  #   LD_LIBRARY_PATH: /opt/qe/lib
  #   OMP_NUM_THREADS: "1"
  pseudo_dir: /path/to/pseudos      # must contain the UPF files below
  timeout: 7200
  pseudopotentials:
    Li: Li.upf
  # k-points: a list, or {{gamma: true}}, or {{density: 0.03}}
  # density = max reciprocal spacing in 1/A (小=密). 变胞时网格只算一次，
  # 取路径上最"吃"k点的那个晶胞，此后所有 image 共用同一套 k 点。
  kpts: [6, 6, 6]
  params:
    ecutwfc: 50
    ecutrho: 400
    occupations: smearing
    # 变胞 NEB 建议 methfessel-paxton 或 marzari-vanderbilt（纯数值展宽），
    # 不要用 gauss/fermi-dirac：后两者是真实展宽，能量含 -TS 项。
    smearing: marzari-vanderbilt
    degauss: 0.02          # Ry
    conv_thr: 1.0e-8

neb:
  images: 5                 # movable images between the two endpoints
  interpolation: idpp      # idpp (recommended) or linear
  align_cell: true          # re-index the final cell (axis order/sign) to match the initial one
  min_distance: 2.0         # drop images whose shortest contact is below this (A); 0 = off
  prune: auto               # auto = drop them and continue, error = stop, false = keep
  variable_cell: true      # VC-NEB: cell is optimised together with atoms
  pressure: 0.0            # GPa; the whole band uses the enthalpy H = E + P*V (leave the DFT input pressure at 0)
  spring_k: 0.1            # eV/A^2
  f_th: 0.03               # convergence: max |F_perp| in eV/A
  max_steps: 300
  dt: 0.1
  max_step: 0.2            # cap on atomic displacement per step (A)
  cell_max_step: 0.02      # cap on cell change per step (A)
  climbing_image: true
  relax_endpoints: true
  endpoint_fmax: 0.02
  parallel_images: 1       # >1 runs several images at once (watch CPU count)

output:
  plot: true
  save_trajectory: true
""",
    "cp2k": """\
project:
  name: li_bcc_to_hcp
  workdir: ./run_li

structure:
  initial: POSCAR_initial
  final: POSCAR_final

calc:
  engine: cp2k
  command: mpirun -np 4 cp2k.popt
  mpi: 4
  # env:               # environment variables, optional
  #   OMP_NUM_THREADS: "1"
  # directory holding BASIS_MOLOPT / GTH_POTENTIALS
  data_dir: /path/to/cp2k/data
  timeout: 7200
  # 金属体系必须给 k 点 + smearing，否则 SCF 不收敛
  kpts: [6, 6, 6]
  params:
    basis_set: DZVP-MOLOPT-SR-GTH
    potential: GTH-PBE-q3
    xc: PBE
    cutoff: 400            # Ry  平面波截断，需自行收敛测试
    rel_cutoff: 50         # Ry
    eps_scf: 1.0e-6
    max_scf: 50
    smearing: true             # 金属打开
    # 可选：fermi_dirac / energy_window / list / gaussian /
    #       methfessel_paxton / marzari_vanderbilt
    # 变胞 NEB 建议 marzari_vanderbilt 或 methfessel_paxton：它们的展宽是
    # 纯数值工具，fermi_dirac 是真实电子温度，能量里含 -TS 项，会给势垒
    # 混进熵贡献。前三个用 SIGMA [eV]，fermi_dirac 用 ELECTRONIC_TEMPERATURE [K]。
    smearing_method: marzari_vanderbilt
    smearing_sigma: 0.05       # eV，金属 0.02-0.1
    # electronic_temperature: 300   # K，仅 fermi_dirac 使用
    added_mos: 8               # 空带数，smearing 时给够
    # ngrids: 5
    # poisson_solver: MT       # 孤立/带电体系可加
    # ---- 自旋 / DFT+U / 布居分析（可选）----
    # uks: true                # 非限制（自旋极化）KS，磁性体系必须打开
    # multiplicity: 3          # 2S+1（此处=三重态）；只有 uks 时才有意义
    # charge: 0                # 体系总电荷
    # hubbard_u: {Ni: 4.0}     # DFT+U 的 U_eff = U-J，单位 eV，按元素给
    # hubbard_l: 2             # 加 U 的角量子数：2 = d，3 = f
    # plus_u_method: mulliken  # mulliken（默认，便宜）/ lowdin / mulliken_charges
    # mulliken: true           # 默认打印 Mulliken 布居（原子 + 净电荷 + 自旋）
    # mulliken_orbitals: false # true 时额外打印每个轨道的贡献（PRINT_GOP）

neb:
  images: 5
  interpolation: idpp
  align_cell: true          # re-index the final cell (axis order/sign) to match the initial one
  min_distance: 2.0         # drop images whose shortest contact is below this (A); 0 = off
  prune: auto               # auto = drop them and continue, error = stop, false = keep
  variable_cell: true
  pressure: 0.0            # GPa; the whole band uses the enthalpy H = E + P*V (leave the DFT input pressure at 0)
  spring_k: 0.1
  f_th: 0.03
  max_steps: 300
  dt: 0.1
  max_step: 0.2
  cell_max_step: 0.02
  climbing_image: true
  relax_endpoints: true
  parallel_images: 1

output:
  plot: true
""",
    "abacus": """\
project:
  name: cu_fcc_to_hcp
  workdir: ./run_cu

structure:
  initial: POSCAR_initial
  final: POSCAR_final

calc:
  engine: abacus
  command: mpirun -np 4 abacus
  mpi: 4
  pseudo_dir: /path/to/PP_ORB     # holds *.upf
  data_dir: /path/to/PP_ORB       # holds *.orb
  timeout: 7200
  pseudopotentials:
    Cu: Cu_ONCV_PBE-1.0.upf
  orbitals:
    Cu: Cu_gga_9au_100Ry_4s2p2d1f.orb
  kpts: [4, 4, 4]      # 或 {density: 0.03}（1/A，小=密），或 {gamma: true}
  params:
    ecutwfc: 100          # Ry，需自行收敛测试
    basis_type: lcao      # lcao（需要 orbitals）或 pw（不需要 orbitals）
    ks_solver: genelpa    # 只有 GPU 才用 'cusolver'
    dft_functional: pbe
    smearing_method: mv   # gauss / fd / mp / mv / fixed；金属用 mp 或 mv
    smearing_sigma: 0.01  # Ry，金属记得给 smearing
    scf_thr: 1.0e-7
    scf_nmax: 100
    symmetry: 0

neb:
  images: 5
  interpolation: idpp
  align_cell: true          # re-index the final cell (axis order/sign) to match the initial one
  min_distance: 2.0         # drop images whose shortest contact is below this (A); 0 = off
  prune: auto               # auto = drop them and continue, error = stop, false = keep
  variable_cell: true
  pressure: 0.0            # GPa; the whole band uses the enthalpy H = E + P*V (leave the DFT input pressure at 0)
  spring_k: 0.1
  f_th: 0.03
  max_steps: 300
  dt: 0.1
  max_step: 0.2
  cell_max_step: 0.02
  climbing_image: true
  relax_endpoints: true
  parallel_images: 1

output:
  plot: true
""",
    "emt": """\
project:
  name: cu_fcc_to_hcp_demo
  workdir: ./run_demo

# Built-in builder: no input files needed, handy to check the workflow.
structure:
  build_initial:
    formula: Cu
    crystal: fcc
    a: 3.60
  build_final:
    formula: Cu
    crystal: hcp
    a: 2.55
    c: 4.17

calc:
  engine: emt

neb:
  images: 5
  interpolation: idpp
  align_cell: true
  min_distance: 2.0
  prune: auto
  variable_cell: true
  pressure: 0.0            # GPa; the whole band uses the enthalpy H = E + P*V (leave the DFT input pressure at 0)
  spring_k: 0.1
  f_th: 0.02
  max_steps: 200
  dt: 0.05
  max_step: 0.1
  cell_max_step: 0.01
  climbing_image: true
  relax_endpoints: true

output:
  plot: true
""",
    "lj": """\
project:
  name: ar_fcc_to_hcp_demo
  workdir: ./run_demo

# Built-in Lennard-Jones: no DFT code, no input files.  Solid argon is a
# textbook fcc -> hcp test case and the whole run takes a few seconds.
structure:
  build_initial:
    formula: Ar
    crystal: fcc
    a: 5.26
    cubic: true                 # 4-atom conventional cell
  build_final:
    formula: Ar
    crystal: hcp
    a: 3.72
    c: 6.08
    orthorhombic: true          # also 4 atoms, so the atom counts match

calc:
  engine: lj
  params:
    epsilon: 0.0104             # eV   (Ar)
    sigma: 3.405                # Angstrom

neb:
  images: 5
  interpolation: idpp
  align_cell: true
  min_distance: 2.0
  prune: auto
  variable_cell: true
  pressure: 0.0            # GPa; the whole band uses the enthalpy H = E + P*V (leave the DFT input pressure at 0)
  spring_k: 0.1
  f_th: 0.002
  max_steps: 400
  dt: 0.1
  max_step: 0.1
  cell_max_step: 0.01
  climbing_image: true
  relax_endpoints: true
  endpoint_fmax: 0.002

output:
  plot: true
""",
}


def default_template(engine: str, fname: str = "neb.yaml") -> str:
    engine = (engine or "lj").lower()
    if engine not in _TEMPLATES:
        raise VcnebError(f"no template for engine '{engine}'; "
                         f"choose from {', '.join(sorted(_TEMPLATES))}")
    return _HEADER.format(fname=fname) + _TEMPLATES[engine]


def write_default_config(path: str, engine: str = "lj") -> str:
    text = default_template(engine, os.path.basename(path))
    ensure_dir(os.path.dirname(os.path.abspath(path)) or ".")
    if os.path.exists(path):
        raise VcnebError(f"{path} already exists; delete it first or pick "
                         f"another output file")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path
