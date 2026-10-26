"""Build the fcc(111) / hcp(111) stacking demo for any element.

    python3 examples/make_111_stack.py Ar 5.26
    python3 examples/make_111_stack.py Cu 3.61

Six close-packed layers in one periodic cell.  fcc = A B C A B C,
hcp = A B A B A B: both structures share the same cell, so a NEB between
them is the classic stacking-fault ("shuffle") path.
"""
import sys

import numpy as np
from ase.build import fcc111
from ase.io import write


def main(argv):
    symbol = argv[0] if argv else "Ar"
    a = float(argv[1]) if len(argv) > 1 else 5.26
    layer_height = a / np.sqrt(3.0)

    fcc = fcc111(symbol, size=(1, 1, 6), a=a, periodic=True, vacuum=None)
    fcc.cell[2] = [0.0, 0.0, 6 * layer_height]
    hcp = fcc.copy()

    frac = fcc.get_scaled_positions()
    layer = np.rint(frac[:, 2] * 6).astype(int) % 6
    a_site, b_site = frac[0, :2], frac[1, :2]
    for i in range(len(frac)):
        frac[i, :2] = a_site if layer[i] % 2 == 0 else b_site
    hcp.set_scaled_positions(frac % 1.0)

    f1 = f"{symbol.lower()}_fcc111.vasp"
    f2 = f"{symbol.lower()}_hcp111.vasp"
    write(f1, fcc, format="vasp", direct=True, vasp5=True)
    write(f2, hcp, format="vasp", direct=True, vasp5=True)
    print(f"wrote {f1} and {f2} ({len(fcc)} atoms, a={a} A)")


if __name__ == "__main__":
    main(sys.argv[1:])
