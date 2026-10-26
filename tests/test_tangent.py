"""The band tangent: plain central difference vs the Henkelman-Jonsson one.

``neb_tangent`` is the only place where the path direction is defined, and the
perpendicular force (the one that actually drives the band onto the minimum
energy path) is built from it.  A wrong tangent does not crash anything, it
just makes the band cut corners of a curved path, so it gets pinned down here:

* straight band  -> both estimators must agree (nothing to cut)
* energy maximum -> improved uses the forward one-sided difference
* energy minimum -> improved uses the backward one-sided difference
* uneven spacing -> the better-sampled side wins
* degenerate band-> a zero vector instead of a NaN

Run with:  PYTHONPATH=. python3 tests/test_tangent.py
"""

import numpy as np

from vcneb.band import Band, neb_tangent
from vcneb.config import Config
from vcneb.util import Logger


def corner():
    """A 90 degree kink: prev -> cur -> next with next not along prev->cur."""
    xs = [np.array([0.0, 0.0]), np.array([1.0, 0.0]), np.array([1.0, 1.0])]
    return xs


# ------------------------------------------------------------------ tangents
def test_straight_band_both_modes_agree():
    xs = [np.array([0.0, 0.0]), np.array([1.0, 0.0]), np.array([2.0, 0.0])]
    energies = [0.0, 0.5, 1.0]
    t_central = neb_tangent(xs, energies, 1, "central")
    t_improved = neb_tangent(xs, energies, 1, "improved")
    assert np.allclose(t_central, [1.0, 0.0])
    assert np.allclose(t_improved, t_central)


def test_upwind_on_monotone_band():
    xs = corner()
    # climbing: only the uphill (forward) side defines the path direction
    assert np.allclose(neb_tangent(xs, [0.0, 0.5, 1.0], 1, "improved"),
                       [0.0, 1.0])
    # descending: only the backward side
    assert np.allclose(neb_tangent(xs, [1.0, 0.5, 0.0], 1, "improved"),
                       [1.0, 0.0])


def test_maximum_weights_the_higher_neighbour():
    xs = corner()
    # image 1 is the maximum; the backward neighbour is 1.0 eV away, the
    # forward one only 0.5 eV -> the backward difference gets the big weight
    t = neb_tangent(xs, [0.0, 1.0, 0.5], 1, "improved")
    expected = np.array([0.5, 1.0]) / np.linalg.norm([0.5, 1.0])
    assert np.allclose(t, expected)
    # ... and that is not the chord (central difference) any more
    t_central = neb_tangent(xs, [0.0, 1.0, 0.5], 1, "central")
    angle = np.degrees(np.arccos(np.clip(np.dot(t, t_central), -1, 1)))
    assert 18.0 < angle < 19.0          # 18.43 deg for this kink


def test_symmetric_maximum_equals_central():
    xs = corner()
    t_improved = neb_tangent(xs, [0.0, 1.0, 0.0], 1, "improved")
    t_central = neb_tangent(xs, [0.0, 1.0, 0.0], 1, "central")
    assert np.allclose(t_improved, t_central)


def test_degenerate_band_is_zero_not_nan():
    xs = [np.array([1.0, 1.0]), np.array([1.0, 1.0]), np.array([1.0, 1.0])]
    for mode in ("central", "improved"):
        t = neb_tangent(xs, [0.0, 1.0, 0.0], 1, mode)
        assert np.all(np.isfinite(t))
        assert np.linalg.norm(t) == 0.0


# -------------------------------------------------------------------- config
def test_config_validates_tangent():
    def problems(tangent):
        cfg = Config({"calc": {"engine": "emt"},
                      "neb": {"images": 2, "tangent": tangent}})
        cfg.structure["initial"] = "/dev/null"
        cfg.structure["final"] = "/dev/null"
        return [p for p in cfg.validate(for_run=False) if "tangent" in p]

    assert problems("improved") == []
    assert problems("central") == []
    assert problems("CENTER")          # typo'd: must be reported
    cfg = Config({"calc": {"engine": "emt"}, "neb": {"images": 2,
                                                     "tangent": "Central"}})
    cfg.validate(for_run=False)
    assert cfg.neb["tangent"] == "central"      # normalised for the run


# ---------------------------------------------------------------------- band
def test_band_uses_both_tangents():
    """neb_forces has to accept either mode and stay finite."""
    from ase.build import bulk

    initial = bulk("Cu", "fcc", a=3.6, cubic=True)
    final = initial.copy()
    final.positions[0] += [0.0, 0.0, 0.4]
    for mode in ("central", "improved"):
        cfg = Config({"calc": {"engine": "emt"},
                      "neb": {"images": 4, "interpolation": "linear",
                              "variable_cell": False, "tangent": mode}})
        band = Band(initial, final, cfg, None, Logger(None))
        band.interpolate()
        profile = np.linspace(0.0, 1.0, band.n_images())
        profile[len(profile) // 2] = 2.0        # a clear maximum in the middle
        for i, st in enumerate(band.states):
            st.energy = float(profile[i])
            st.forces = np.zeros((band.n_atoms, 3))
            st.stress = None
        forces, diag = band.neb_forces(ci_index=2)
        assert all(np.all(np.isfinite(f)) for f in forces)
        assert np.all(forces[0] == 0.0) and np.all(forces[-1] == 0.0)
        assert diag["f_perp_atom"] and np.all(np.isfinite(diag["f_perp_atom"]))


def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("PASS", name)


if __name__ == "__main__":
    main()
