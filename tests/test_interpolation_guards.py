"""Path-building guards: cell re-indexing (align_cell) and image pruning.

Both guards are about *not* evaluating unphysical images:

* ``_align_final_cell`` -- the final cell may be written with its lattice
  vectors in another order / direction.  That is the same crystal, but a
  different fractional coordinate system, so interpolating in fractional
  coordinates swaps the atoms into each other.
* ``_prune_close_images`` -- if two atoms still end up closer than
  ``neb.min_distance`` the image is dropped instead of being handed to the
  DFT code (the same "throw away bad candidates" idea as USPEX).

Run with:  PYTHONPATH=. python3 tests/test_interpolation_guards.py
"""

import os
import tempfile

import numpy as np

from ase.build import bulk

from vcneb.band import Band
from vcneb.config import Config
from vcneb.util import Logger, VcnebError


def make_band(initial, final, raw=None, workdir=None):
    data = {
        "project": {"workdir": workdir or tempfile.mkdtemp()},
        "calc": {"engine": "emt"},
        "neb": {"images": 4, "interpolation": "linear"},
    }
    for section, over in (raw or {}).items():
        data.setdefault(section, {}).update(over)
    cfg = Config(data)
    return Band(initial, final, cfg, None, Logger(None))


def clean_band(raw=None, workdir=None):
    """fcc -> slightly compressed fcc: nothing moves, nothing collides."""
    initial = bulk("Cu", "fcc", a=3.6, cubic=True)
    final = bulk("Cu", "fcc", a=3.2, cubic=True)
    final.set_cell(initial.get_cell(), scale_atoms=False)
    band = make_band(initial, final, raw=raw, workdir=workdir)
    band.interpolate()
    return band


# --------------------------------------------------------------- align_cell
def test_align_cell():
    a = 3.6
    initial = bulk("Cu", "fcc", a=a, cubic=True)
    final = initial.copy()
    # same crystal, lattice vectors in another order and one direction
    # flipped -- a different description of the identical structure
    perm = [2, 0, 1]
    signs = np.array([1.0, -1.0, 1.0])
    frac = np.array(final.get_scaled_positions())
    final.set_cell(np.array(final.get_cell())[perm] * signs[:, None],
                   scale_atoms=False)
    final.set_scaled_positions(frac[:, perm] * signs)
    pos_before = final.get_positions()

    band = make_band(initial, final)
    band.interpolate()

    np.testing.assert_allclose(band.states[-1].cell,
                               np.array(initial.get_cell()), atol=1e-8)
    # the final structure itself never moved, only its description
    np.testing.assert_allclose(band.states[-1].pos, pos_before, atol=1e-8)
    assert band.n_free == 4, "nothing should have been dropped"
    assert band.min_pair_distance() > 2.4, band.min_pair_distance()
    return band


def test_align_cell_off_is_worse():
    """Without the guard the very same case gets atoms within 2 A."""
    initial = bulk("Cu", "fcc", a=3.6, cubic=True)
    final = initial.copy()
    perm = [2, 0, 1]
    signs = np.array([1.0, -1.0, 1.0])
    frac = np.array(final.get_scaled_positions())
    final.set_cell(np.array(final.get_cell())[perm] * signs[:, None],
                   scale_atoms=False)
    final.set_scaled_positions(frac[:, perm] * signs)
    band = make_band(initial, final,
                     raw={"neb": {"align_cell": False, "min_distance": 0}})
    band.interpolate()
    assert band.min_pair_distance() < 2.0, band.min_pair_distance()
    return band


# ------------------------------------------------------------------- pruning
def test_prune_partial():
    band = clean_band(raw={"neb": {"min_distance": 2.0, "prune": "auto"}})
    clean = band.min_pair_distance()
    assert clean > 2.0, clean
    # force an overlap in the middle image only
    band.states[2].pos[1] = band.states[2].pos[0]
    band._prune_close_images()
    assert band.pruned == [2], band.pruned
    assert band.n_free == 3, band.n_free
    assert band.n_images() == 5
    assert band.min_pair_distance() >= 2.0
    return band


def test_prune_error_mode():
    band = clean_band(raw={"neb": {"min_distance": 2.0, "prune": "error"}})
    band.states[2].pos[1] = band.states[2].pos[0]
    try:
        band._prune_close_images()
    except VcnebError as exc:
        assert "min_distance" in str(exc), str(exc)
    else:
        raise AssertionError("prune: error should have raised")


def test_prune_all_is_an_error():
    band = clean_band(raw={"neb": {"min_distance": 2.0, "prune": "auto"}})
    for st in band.states[1:-1]:
        st.pos[1] = st.pos[0]
    try:
        band._prune_close_images()
    except VcnebError as exc:
        assert "endpoint" in str(exc), str(exc)
    else:
        raise AssertionError("an all-bad path should raise, not silently "
                             "collapse to the endpoints")


def test_min_distance_zero_disables():
    band = clean_band(raw={"neb": {"min_distance": 0}})
    for st in band.states[1:-1]:
        st.pos[1] = st.pos[0]
    band._prune_close_images()
    assert band.n_free == 4
    assert band.pruned == []


# -------------------------------------------------------------- restart path
def test_state_roundtrip_keeps_band_size():
    workdir = tempfile.mkdtemp()
    band = clean_band(raw={"neb": {"min_distance": 2.0}}, workdir=workdir)
    band.states[2].pos[1] = band.states[2].pos[0]
    band._prune_close_images()
    band.save_state()

    # the config still asks for 4 movable images -- the state wins, because
    # the image count cannot change in the middle of a run
    other = make_band(band.atoms(0), band.atoms(-1),
                      raw={"neb": {"images": 4}}, workdir=workdir)
    other.load_state(os.path.join(workdir, "state.json"))
    assert other.n_free == 3, other.n_free
    assert other.n_images() == 5
    assert other.pruned == [2], other.pruned


def main():
    band = test_align_cell()
    print(f"align_cell        : OK  (cell restored, min distance "
          f"{band.min_pair_distance():.3f} A)")
    band = test_align_cell_off_is_worse()
    print(f"align_cell: false : OK  (would give "
          f"{band.min_pair_distance():.2f} A contacts)")
    band = test_prune_partial()
    print(f"prune partial     : OK  (dropped {band.pruned}, "
          f"{band.n_free} movable left)")
    test_prune_error_mode()
    print("prune: error      : OK")
    test_prune_all_is_an_error()
    print("prune all-bad     : OK  (raises instead of collapsing)")
    test_min_distance_zero_disables()
    print("min_distance: 0   : OK  (guard off)")
    test_state_roundtrip_keeps_band_size()
    print("state roundtrip   : OK")
    print("\nall interpolation guards passed")


if __name__ == "__main__":
    main()
