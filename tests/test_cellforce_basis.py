"""Why the cell force is ``-V sigma h^-T`` and not "its exact Cartesian form".

Two charts describe the same energy E(h, s):

  (a) fractional coordinates held fixed  ->  dE/dh = -V sigma h^-T
  (b) Cartesian coordinates held fixed   ->  dE/dh = -V sigma h^-T - s^T f

vcneb stores Cartesian positions, so (b) looks like "the exact one".  It is
not what we want: with free atoms the extra term ``-s^T f`` has the magnitude
of the atom forces (it is ~6x the stress term while the atoms are not
relaxed), so the cell would follow the atoms instead of the stress.

Measured on the EMT fcc->hcp example (5 images, spring_k=0.1, f_th=0.01):

    chart (a)  376 steps, cells 88-93 deg (i.e. unchanged), barrier 0.022 eV
    chart (b)  502 steps, cells up to 100 deg (sheared),     barrier 0.083 eV

Concluding: chart (a) is what the band uses (see ``cell_force_from_stress``),
which is also the usual variable-cell NEB convention.  This script re-checks
that decision against finite differences.

Run with:  PYTHONPATH=. python3 tests/test_cellforce_basis.py
"""

import numpy as np

from ase.build import bulk
from ase.calculators.emt import EMT

from vcneb.band import cell_force_from_stress, voigt_to_matrix


def evaluate(cell, pos, template):
    atoms = template.copy()
    atoms.set_cell(cell, scale_atoms=False)
    atoms.set_positions(pos)
    atoms.calc = EMT()
    return (atoms.get_potential_energy(),
            atoms.get_forces(),
            voigt_to_matrix(atoms.get_stress(voigt=True)))


def main():
    template = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = np.array(template.get_cell()) @ np.array(
        [[1.010, 0.005, -0.004],
         [0.003, 0.988, 0.007],
         [-0.006, 0.004, 1.018]])
    frac = np.array(template.get_scaled_positions())
    frac[2] += np.array([0.06, -0.03, 0.04])  # break the symmetry: large forces
    pos = frac @ cell

    e0, f, sigma = evaluate(cell, pos, template)
    f_cell = cell_force_from_stress(sigma, cell)
    cart_term = frac.T @ f  # the term chart (b) would add

    delta = 1e-4
    fd_cart = np.zeros((3, 3))
    fd_frac = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            step = np.zeros((3, 3))
            step[a, b] = delta
            # Cartesian positions held fixed
            ep = evaluate(cell + step, pos, template)[0]
            em = evaluate(cell - step, pos, template)[0]
            fd_cart[a, b] = -(ep - em) / (2 * delta)
            # fractional positions held fixed
            ep = evaluate(cell + step, frac @ (cell + step), template)[0]
            em = evaluate(cell - step, frac @ (cell - step), template)[0]
            fd_frac[a, b] = -(ep - em) / (2 * delta)

    print(f"max |f_atom|                 = {np.abs(f).max():.4f} eV/A")
    print(f"V                            = {abs(np.linalg.det(cell)):.3f} A^3")
    print(f"|F_h| = |-V sigma h^-T|      = {np.abs(f_cell).max():.4f} eV/A")
    print(f"|s^T f| (the extra term)     = {np.abs(cart_term).max():.4f} eV/A")
    print(f"stress vs frac-FD            = "
          f"{np.abs(f_cell - fd_frac).max():.3e}   <- must be small")
    print(f"stress vs Cartesian-FD       = "
          f"{np.abs(f_cell - fd_cart).max():.3e}   <- is not, by design")
    print(f"stress + s^T f vs Cart-FD    = "
          f"{np.abs(f_cell - cart_term - fd_cart).max():.3e}")

    # EMT's analytic stress and the energy derivative differ by ~2% (delta
    # independent, i.e. not a finite-difference error), so 5% is the tolerance.
    scale = max(np.abs(fd_frac).max(), 1e-12)
    err_frac = np.abs(f_cell - fd_frac).max() / scale
    scale_cart = max(np.abs(fd_cart).max(), 1e-12)
    err_cart = np.abs(f_cell - cart_term - fd_cart).max() / scale_cart
    print()
    if err_frac < 5e-2 and err_cart < 5e-2:
        print("PASS: -V sigma h^-T is the fixed-fractional force (chart a, what "
              "the band uses); adding -s^T f would give chart (b)")
    else:
        print("FAIL: neither chart matches the finite differences")
    return 0 if err_frac < 5e-2 and err_cart < 5e-2 else 1


if __name__ == "__main__":
    raise SystemExit(main())
