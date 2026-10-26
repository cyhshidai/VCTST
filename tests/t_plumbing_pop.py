"""LJ NEB with a fake population analysis bolted on.

Exercises the ImageState.mulliken -> barrier.csv -> summary.json -> plot
plumbing without spending DFT time.  Run with
``PYTHONPATH=. python3 tests/t_plumbing_pop.py``.
"""

import numpy as np

from ase.calculators.calculator import Calculator, all_changes

from vcneb import cli, engines


class PopCalc(Calculator):
    """Real ASE calculator: LJ plus a fake Mulliken table."""

    implemented_properties = ("energy", "forces", "stress")
    name = "lj-pop"

    def __init__(self, inner):
        super().__init__()
        self.inner = inner

    def calculate(self, atoms=None, properties=("energy",),
                  system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.inner.calculate(self.atoms, properties, system_changes)
        self.results = dict(self.inner.results)
        spin = 2.0 * abs(float(np.linalg.det(self.atoms.get_cell()))) / 218.6
        self.results["mulliken"] = {
            "atoms": [{"index": i, "element": s, "kind": 1,
                       "alpha": 0.0, "beta": 0.0, "charge": 0.0, "spin": 0.0}
                      for i, s in enumerate(self.atoms.get_chemical_symbols())],
            "total": {"alpha": spin, "beta": 0.0, "charge": 0.0, "spin": spin},
        }


class FakePopLJ(engines.LJEngine):
    def calculator(self, atoms, image_dir):
        return PopCalc(super().calculator(atoms, image_dir))


engines._REGISTRY["lj"] = FakePopLJ
raise SystemExit(cli.main(["run", "examples/lj_ar.yaml", "--force"]))
