"""neb.pressure offsets the energy (H = E + P*V) and the stress (+P on the
diagonal), so the band and the endpoint relaxation optimise the enthalpy."""
import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.optimize import FIRE
from ase.filters import ExpCellFilter

from vcneb.band import EnthalpyCalculator, GPA_TO_EV_A3


def test_offset_is_consistent():
    atoms = bulk("Cu", "fcc", a=3.6).repeat((2, 2, 2))
    atoms.calc = EMT()
    e0 = atoms.get_potential_energy()
    s0 = atoms.get_stress(voigt=True)
    f0 = atoms.get_forces()
    vol = atoms.get_volume()

    pressure = 2.55  # GPa
    atoms.calc = EnthalpyCalculator(EMT(), pressure * GPA_TO_EV_A3)
    e1 = atoms.get_potential_energy()
    s1 = atoms.get_stress(voigt=True)
    f1 = atoms.get_forces()

    assert abs(e1 - (e0 + pressure * GPA_TO_EV_A3 * vol)) < 1e-9
    assert np.allclose(s1[:3], s0[:3] + pressure * GPA_TO_EV_A3)
    assert np.allclose(s1[3:], s0[3:])
    assert np.allclose(f1, f0)


def test_cell_compresses_under_pressure():
    """ExpCellFilter + the offset must shrink the cell: dH/dV = 0 <=> P = -dE/dV."""
    a0 = 3.6
    atoms = bulk("Cu", "fcc", a=a0)
    atoms.calc = EMT()
    FIRE(ExpCellFilter(atoms), logfile=None).run(fmax=1e-4, steps=200)
    v0 = atoms.get_volume()

    for pressure, expect_shrink in ((10.0, True), (-10.0, False)):
        atoms = bulk("Cu", "fcc", a=a0)
        atoms.calc = EnthalpyCalculator(EMT(), pressure * GPA_TO_EV_A3)
        FIRE(ExpCellFilter(atoms), logfile=None).run(fmax=1e-4, steps=400)
        v = atoms.get_volume()
        if expect_shrink:
            assert v < v0, (pressure, v, v0)
        else:
            assert v > v0, (pressure, v, v0)
