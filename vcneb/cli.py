"""Command line interface: init / doctor / run / plot / status."""

from __future__ import annotations

import argparse
import os
import shutil
import sys

from . import __version__
from .config import DEFAULTS, ENGINE_NAMES, load_config, write_default_config
from .util import Logger, VcnebError, ensure_dir, jload, print_exc_brief


RUN_OUTPUTS = ("state.json", "history.jsonl", "summary.json", "barrier.csv",
               "band.traj", "barrier.png", "config.resolved.yaml")


#: engine keys the tool writes itself, so users do not have to (and should
#: not) repeat them.  Derived from the calculator()/render code paths.
_TOOL_KEYS = {
    "qe": ("calculation", "outdir", "prefix", "pseudo_dir", "disk_io",
           "tprnfor", "tstress", "verbosity", "kpts", "pseudopotentials"),
    "cp2k": ("&GLOBAL project / run_type / print_level",
             "&DFT basis_set_file_name / potential_file_name",
             "&KPOINTS (from calc.kpts)", "&SCF scf_guess",
             "&PRINT &FORCES / &STRESS_TENSOR",
             "&SUBSYS &CELL / &COORD / &KIND"),
    "abacus": ("pp", "basis", "kpts", "basis_type", "ks_solver", "symmetry",
               "scf_thr", "scf_nmax", "cal_force", "cal_stress", "device",
               "ntype"),
}


def _print_problems(problems: list[str], prefix: str = "problem") -> None:
    for i, problem in enumerate(problems, 1):
        print(f"  [{i}] {problem}")


# -------------------------------------------------------------------- init
def cmd_init(args) -> int:
    engine = args.engine.lower()
    if engine not in ENGINE_NAMES:
        raise VcnebError(f"unknown engine '{engine}'; "
                         f"choose from {', '.join(ENGINE_NAMES)}")
    path = write_default_config(args.output, engine)
    print(f"wrote {path}")
    print("\nnext steps:")
    print(f"  1) edit {path} (structures, engine command, pseudopotentials)")
    print(f"  2) vcneb doctor {path}")
    print(f"  3) vcneb run {path}")
    return 0


# ----------------------------------------------------------------- doctor
def cmd_doctor(args) -> int:
    cfg = load_config(args.config)
    ok = True
    print(f"vcneb {__version__} doctor: {args.config}\n")

    print("[dependencies]")
    for mod, why in (("numpy", "arrays"), ("ase", "atoms/calculators"),
                     ("yaml", "config file"), ("matplotlib", "barrier plot")):
        try:
            module = __import__(mod)
            version = getattr(module, "__version__", "?")
            print(f"  OK    {mod:<10} {version}")
        except ImportError:
            level = "WARN" if mod == "matplotlib" else "FAIL"
            ok = ok and level != "FAIL"
            print(f"  {level}  {mod:<10} missing (needed for {why}) -> "
                  f"pip install {mod if mod != 'yaml' else 'pyyaml'}")

    print("\n[config]")
    problems = cfg.validate(for_run=False)
    if problems:
        ok = False
        _print_problems(problems)
    else:
        print("  OK    config looks consistent")

    print("\n[structures]")
    try:
        from .structs import get_endpoints

        initial, final = get_endpoints(cfg)
        print(f"  OK    initial: {len(initial)} atoms, "
              f"V={initial.get_volume():.4f} A^3, "
              f"formula={initial.get_chemical_formula()}")
        print(f"  OK    final  : {len(final)} atoms, "
              f"V={final.get_volume():.4f} A^3, "
              f"formula={final.get_chemical_formula()}")
        for key in ("initial", "final"):
            value = cfg.structure.get(key)
            if value:
                print(f"        {key}: {value}")
            else:
                print(f"        {key}: built from config "
                      f"({cfg.structure.get('build_' + key)})")
    except VcnebError as exc:
        ok = False
        print(f"  FAIL  {exc}")

    print("\n[engine]")
    symbols = None
    try:
        from .structs import get_endpoints

        symbols = sorted(set(get_endpoints(cfg)[0].get_chemical_symbols()))
    except Exception:  # noqa: BLE001 - structures may be broken, reported above
        pass
    try:
        from .engines import make_engine

        engine = make_engine(cfg)
        print(f"  {engine.describe()}")
        engine_problems = engine.check()
        if engine_problems:
            ok = False
            _print_problems(engine_problems)
        else:
            print("  OK    engine configuration looks usable")
        for note in getattr(engine, "warnings", list)(symbols):
            print(f"  WARN  {note}")
    except VcnebError as exc:
        ok = False
        print(f"  FAIL  {exc}")

    print()
    if ok:
        print("all checks passed -> vcneb run " + args.config)
    else:
        print("some checks failed, fix them before running")
    return 0 if ok else 1


# -------------------------------------------------------------------- run
# ------------------------------------------------------------------- keys
def _print_engine_keys(engine: str) -> int:
    """List the calc.params keys an engine understands."""
    if engine == "qe":
        from ase.io.espresso_namelist.keys import pw_keys

        tool = set(_TOOL_KEYS["qe"])
        print("QE pw.x namelist keywords known to ASE.")
        print("calc.params keys are sorted into these namelists for you;")
        print("anything NOT listed here is dropped silently (typos!)")
        for namelist, keys in pw_keys.items():
            print(f"\n&{namelist.upper()}")
            for key in keys:
                mark = "   <- default set by vcneb" if key in tool else ""
                print(f"  {key}{mark}")
        print("\nauthoritative source: Modules/input_parameters.f90")
        print("  (NAMELIST / control / system / electrons / ions / cell / "
              "ions / press ...)")
    elif engine == "cp2k":
        from .engines import _CP2K_PARAM_KEYS

        print("CP2K keys that the built-in input writer renders:")
        for key in sorted(_CP2K_PARAM_KEYS):
            print(f"  {key}")
        print("\nevery other calc.params key is ignored (doctor warns about "
              "it).")
        print("inject raw lines into &DFT instead:")
        print("  params:\n    extra: |\n      &SMEAR ON\n"
              "        METHOD FERMI_DIRAC\n"
              "        ELECTRONIC_TEMPERATURE [K] 300\n      &END SMEAR")
        print("\nskeleton written by vcneb (do not repeat in params):")
        for item in _TOOL_KEYS["cp2k"]:
            print(f"  {item}")
    elif engine == "abacus":
        print("ABACUS has no whitelist: every calc.params key is written "
              "verbatim\ninto INPUT (ase/io/abacus.py: write_input), so "
              "any INPUT keyword works.")
        print("\nwritten by vcneb, do not repeat: "
              + ", ".join(_TOOL_KEYS["abacus"]))
        print("\nkeys that take several values are space separated inside one "
              "yaml string:")
        print('  orbital_corr: "2 2 -1"\n  hubbard_u: "4.0 4.0 0.0"')
        print("\nfull list from your ABACUS source tree:")
        print("  grep -oP '^### \\K.*' "
              "docs/advanced/input_files/input-main.md")
    else:
        print(f"{engine}: no calc.params keywords (no DFT input file)")
    return 0


def cmd_keys(args) -> int:
    target = (args.target or "").strip().lower()
    print(f"vcneb {__version__} keys\n")
    if target in ENGINE_NAMES:
        return _print_engine_keys(target)
    if target and target not in DEFAULTS:
        raise VcnebError(
            f"unknown target '{args.target}'; use one of "
            + ", ".join(list(DEFAULTS) + list(ENGINE_NAMES)))

    import yaml

    for name in ([target] if target else list(DEFAULTS)):
        print(f"[{name}]")
        print(yaml.safe_dump(DEFAULTS[name], sort_keys=False,
                             default_flow_style=False).rstrip())
        print()
    if not target:
        print("calc.params is engine specific:  vcneb keys qe | cp2k | abacus")
        print("commented template:              vcneb init -e cp2k -o neb.yaml")
        print("what actually ran:               <workdir>/config.resolved.yaml")
    return 0


# -------------------------------------------------------------------- run
def _clean_workdir(cfg) -> None:
    """Remove only the files this tool owns; never wipe an unknown folder."""
    target = os.path.abspath(cfg.workdir)
    guard = os.path.dirname(target) or target
    if not target.startswith(os.path.abspath(guard)):
        raise VcnebError(f"refusing to clean suspicious workdir {target}")
    for name in ("images", "logs"):
        path = os.path.join(target, name)
        if os.path.isdir(path):
            shutil.rmtree(path)
    for name in RUN_OUTPUTS + ("endpoint_initial.vasp", "endpoint_final.vasp"):
        path = os.path.join(target, name)
        if os.path.isfile(path):
            os.remove(path)


def cmd_run(args) -> int:
    cfg = load_config(args.config)
    if args.max_steps:
        cfg.neb["max_steps"] = int(args.max_steps)
    problems = cfg.validate()
    if problems:
        print("cannot run: the configuration has problems\n", file=sys.stderr)
        _print_problems(problems)
        print("\ntip: `vcneb doctor " + args.config + "` explains each item",
              file=sys.stderr)
        return 1

    from .band import Band
    from .engines import make_engine
    from .structs import get_endpoints

    ensure_dir(cfg.workdir)
    ensure_dir(cfg.logs_dir)
    ensure_dir(cfg.images_dir)
    log = Logger(os.path.join(cfg.logs_dir, "vcneb.log"))
    if args.verbose:
        os.environ["VCNEB_VERBOSE"] = "1"

    log(f"vcneb {__version__}  config: {os.path.abspath(args.config)}")
    log(f"engine: {make_engine(cfg).describe()}")
    for line in cfg.summary_lines():
        log("  " + line)

    state_path = os.path.join(cfg.workdir, "state.json")
    restart = bool(args.restart)
    if restart and not os.path.isfile(state_path):
        raise VcnebError(
            f"--restart was given but {state_path} does not exist; "
            f"run once without --restart first")
    if os.path.isfile(state_path) and not restart:
        if not args.force:
            raise VcnebError(
                f"{cfg.workdir} already contains a NEB run.\n"
                f"  use `--restart` to continue it, or `--force` to start over")
        log("--force: removing the previous run in " + cfg.workdir)
        _clean_workdir(cfg)
        ensure_dir(cfg.logs_dir)
        ensure_dir(cfg.images_dir)

    initial, final = get_endpoints(cfg)
    engine = make_engine(cfg)
    # Fix the k-mesh once for every image (see BaseEngine.plan_kpoints):
    # a cell that shrinks along the path would otherwise need a denser grid,
    # and every change of the grid is a step in the energy -> fake force.
    engine.plan_kpoints([initial.get_cell(), final.get_cell()])
    log(f"k-points: {engine.kpts_text()}   (same grid for all images)")
    engine_problems = engine.check()
    if engine_problems and not args.skip_engine_check:
        if args.dry_run:
            # a dry run only builds the path: let the user inspect structures
            # and interpolation even before the DFT code is installed
            print("warning: engine is not ready (dry run continues anyway):\n",
                  file=sys.stderr)
            _print_problems(engine_problems)
            print("", file=sys.stderr)
        else:
            print("engine is not ready:\n", file=sys.stderr)
            _print_problems(engine_problems)
            print("\ntip: `vcneb doctor " + args.config + "`", file=sys.stderr)
            return 1
    log(f"endpoints: {len(initial)} atoms, "
        f"V0={initial.get_volume():.4f} -> V1={final.get_volume():.4f} A^3")
    _warn_resources(cfg, log)

    band = Band(initial, final, cfg, engine, log, dry_run=bool(args.dry_run))
    if restart:
        band.load_state(state_path)
    elif args.dry_run:
        band.relax_endpoints()
        band.interpolate()
        _dry_run_report(band, cfg, log)
        return 0
    else:
        band.relax_endpoints()
        band.interpolate()
        band.save_state()

    try:
        summary = band.run()
    except KeyboardInterrupt:
        log.warn("interrupted; progress is in state.json, "
                 "continue with --restart")
        return 130
    _write_resolved_config(cfg, args.config)
    _maybe_clean_images(cfg, band, log)
    print()
    print(f"barrier   : {summary['barrier_eV']:.4f} eV")
    print(f"reverse   : {summary['reverse_barrier_eV']:.4f} eV")
    print(f"steps     : {summary['steps']}  "
          f"converged: {summary['converged']}")
    print(f"results   : {cfg.workdir}")
    return 0 if summary["converged"] else 2


def _dry_run_report(band, cfg, log) -> None:
    log("dry run: path built, no electronic-structure call was made")
    log(f"  atoms per image: {band.n_atoms}   "
        f"variable cell: {band.variable_cell}")
    for i in range(band.n_images()):
        atoms = band.atoms(i)
        log(f"  image {i:02d}: V={atoms.get_volume():.4f} A^3  "
            f"formula={atoms.get_chemical_formula()}  "
            f"dir={band.image_dir(i)}")
    try:
        from ase.io import write

        for tag, idx in (("initial", 0), ("final", band.n_images() - 1)):
            path = os.path.join(cfg.workdir, f"endpoint_{tag}.vasp")
            # no sort=True: the atom order must stay comparable by eye
            write(path, band.atoms(idx), format="vasp", direct=True)
            log(f"  wrote {path} (look at it in VESTA/OVITO)")
    except Exception as exc:  # noqa: BLE001
        log.warn(f"could not write the endpoint POSCARs: {exc}")
    log("next: remove --dry-run (or keep it to inspect the path)")


def _maybe_clean_images(cfg, band, log) -> None:
    """Delete per-image scratch dirs when neb.keep_intermediate_images is off."""
    if cfg.neb.get("keep_intermediate_images", True):
        return
    removed = 0
    for i in range(1, band.n_images() - 1):
        path = band.image_dir(i)
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
            removed += 1
    if removed:
        log(f"neb.keep_intermediate_images=false: removed {removed} image "
            f"directories (energies/forces stay in summary.json and band.traj)")


def _warn_resources(cfg, log) -> None:
    """Warn when the MPI/parallel settings ask for more processes than cores."""
    mpi = int(cfg.calc.get("mpi") or 1)
    par = int(cfg.neb.parallel_images)
    cores = os.cpu_count() or 1
    if mpi > 1 and par > 1 and mpi * par > cores:
        log.warn(f"calc.mpi ({mpi}) x neb.parallel_images ({par}) = "
                 f"{mpi * par} processes on {cores} cores; lower one of them")


def _write_resolved_config(cfg, source: str) -> None:
    try:
        import yaml

        with open(os.path.join(cfg.workdir, "config.resolved.yaml"), "w",
                  encoding="utf-8") as fh:
            yaml.safe_dump(cfg.data, fh, sort_keys=False,
                           default_flow_style=False)
    except Exception:  # noqa: BLE001
        pass


# ------------------------------------------------------------------- plot
def cmd_plot(args) -> int:
    from .band import plot_band

    out = plot_band(args.run_dir, out=args.out)
    print(f"wrote {out}")
    return 0


# ----------------------------------------------------------------- status
def cmd_status(args) -> int:
    run_dir = args.run_dir
    summary = jload(os.path.join(run_dir, "summary.json"))
    state = jload(os.path.join(run_dir, "state.json"))
    if summary:
        print(f"run     : {run_dir}")
        print(f"engine  : {summary.get('engine')}")
        if summary.get("tangent"):
            print(f"tangent : {summary['tangent']}")
        print(f"status  : {'converged' if summary.get('converged') else 'not converged'}"
              f" after {summary.get('steps')} steps")
        print(f"barrier : {summary['barrier_eV']:.4f} eV "
              f"(reverse {summary['reverse_barrier_eV']:.4f} eV)")
        print()
        print(" image      E (eV)      dE (eV)    max|F| (eV/A)     V (A^3)")
        for row in summary["images"]:
            print(f" {row['image']:5d}  {row['energy_eV']:11.5f}  "
                  f"{row['delta_eV']:11.5f}  "
                  f"{row['max_force_eV_per_A']:13.4f}  "
                  f"{row['volume_A3']:10.4f}")
        if not summary.get("converged"):
            steps = int(summary.get("steps") or 0)
            print()
            print("not converged yet; continue with")
            print(f"  vcneb run <config> --restart --max-steps {steps + 500}")
        return 0
    if state:
        print(f"run     : {run_dir}")
        print(f"progress: step {state.get('step')} of an unfinished run")
        if state.get("tangent"):
            print(f"tangent : {state['tangent']}")
        print(f"resume  : vcneb run <config> --restart")
        return 0
    raise VcnebError(f"{run_dir} does not look like a vcneb run directory")


# ------------------------------------------------------------------ parser
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vcneb",
        description="Turnkey variable-cell NEB for QE / CP2K / ABACUS. "
                    "One YAML file in, barrier plot out.")
    parser.add_argument("--version", action="version",
                        version=f"vcneb {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_init = sub.add_parser("init", help="write a commented config template")
    p_init.add_argument("-o", "--output", default="neb.yaml")
    p_init.add_argument("-e", "--engine", default="lj",
                        choices=list(ENGINE_NAMES))
    p_init.set_defaults(func=cmd_init)

    p_doc = sub.add_parser("doctor", help="check everything before spending CPU")
    p_doc.add_argument("config")
    p_doc.set_defaults(func=cmd_doctor)

    p_keys = sub.add_parser("keys",
                            help="list yaml fields / engine keywords")
    p_keys.add_argument("target", nargs="?", default=None,
                        help="a config section (project/structure/calc/"
                             "neb/output) or an engine (qe/cp2k/abacus)")
    p_keys.set_defaults(func=cmd_keys)

    p_run = sub.add_parser("run", help="relax endpoints, build and run the NEB")
    p_run.add_argument("config")
    p_run.add_argument("--restart", action="store_true",
                       help="continue an interrupted run")
    p_run.add_argument("--force", action="store_true",
                       help="delete the previous run in workdir and start over")
    p_run.add_argument("--dry-run", action="store_true",
                       help="build and print the initial path, run nothing")
    p_run.add_argument("--max-steps", type=int, default=None,
                       help="override neb.max_steps")
    p_run.add_argument("--skip-engine-check", action="store_true",
                       help="do not stop on engine configuration problems")
    p_run.add_argument("-v", "--verbose", action="store_true")
    p_run.set_defaults(func=cmd_run)

    p_plot = sub.add_parser("plot", help="redraw barrier.png from summary.json")
    p_plot.add_argument("run_dir")
    p_plot.add_argument("-o", "--out", default=None)
    p_plot.set_defaults(func=cmd_plot)

    p_status = sub.add_parser("status", help="show the state of a run")
    p_status.add_argument("run_dir")
    p_status.set_defaults(func=cmd_status)
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args) or 0)
    except VcnebError as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130
    except Exception as exc:  # noqa: BLE001
        if os.environ.get("VCNEB_VERBOSE"):
            raise
        print(f"\nUNEXPECTED ERROR: {print_exc_brief(exc)}\n"
              f"re-run with VCNEB_VERBOSE=1 for a traceback", file=sys.stderr)
        return 1
