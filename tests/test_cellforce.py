"""Cell force vs finite differences with *fractional* positions held fixed.

This is the branch where ``-V sigma h^-T`` is exact (it is what ``stress``
means).  The Cartesian branch used by the band is checked separately in
tests/test_cellforce_basis.py.

Run with:  PYTHONPATH=. python3 tests/test_cellforce.py
"""

import numpy as np

from ase.build import bulk
from ase.calculators.emt import EMT

from vcneb.band import cell_force_from_stress, voigt_to_matrix


def energy_and_stress(cell, frac, template):
    atoms = template.copy()
    atoms.set_cell(cell, scale_atoms=False)
    atoms.set_scaled_positions(frac)
    atoms.calc = EMT()
    return atoms.get_potential_energy(), atoms.get_stress(voigt=True)


def main():
    template = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = np.array(template.get_cell())
    # a generic distortion so that every component is non-zero
    cell = cell @ np.array([[1.012, 0.004, -0.003],
                            [0.002, 0.987, 0.006],
                            [-0.005, 0.003, 1.021]])
    frac = np.array(template.get_scaled_positions())
    energy, stress6 = energy_and_stress(cell, frac, template)
    stress = voigt_to_matrix(stress6)
    f_cell = cell_force_from_stress(stress, cell)

    delta = 1e-4
    fd = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            cp = cell.copy()
            cm = cell.copy()
            cp[i, j] += delta
            cm[i, j] -= delta
            ep, _ = energy_and_stress(cp, frac, template)
            em, _ = energy_and_stress(cm, frac, template)
            fd[i, j] = -(ep - em) / (2 * delta)
    print(f"energy = {energy:.8f} eV, volume = {np.linalg.det(cell):.4f}")
    print(f"{'component':>10} {'-dE/dh (FD)':>14} {'F_cell':>13} {'diff':>10}")
    for i in range(3):
        for j in range(3):
            print(f"{i}{j:>9} {fd[i, j]:14.5f} {f_cell[i, j]:13.5f} "
                  f"{f_cell[i, j] - fd[i, j]:10.5f}")

    # The strain conjugate to the stress is symmetric, so the analytic force is
    # only defined up to an antisymmetric (pure rotation) part.  Two checks:
    #  1) diagonal components (no ambiguity) must match the finite difference,
    #  2) the residual must indeed be antisymmetric,
    #  3) stepping along +F_cell must lower the energy.
    scale = max(np.abs(np.diag(fd)).max(), 1e-12)
    diag_err = np.abs(np.diag(f_cell) - np.diag(fd)).max() / scale
    resid = f_cell - fd
    anti_err = np.abs(resid + resid.T).max() / scale
    step = 1e-3
    new_cell = cell + step * f_cell
    e_new, _ = energy_and_stress(new_cell, frac, template)
    descent = e_new - energy
    print()
    print(f"diagonal error        : {diag_err:.2e}  (want < 1e-3)")
    print(f"antisymmetry of residuar: {anti_err:.2e}  (want < 1e-2)")
    print(f"E(cell + 1e-3*F) - E   : {descent:.3e} eV  (want < 0)")
    ok = diag_err < 1e-3 and anti_err < 1e-2 and descent < 0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
