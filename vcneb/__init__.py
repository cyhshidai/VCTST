"""vcneb - a turnkey variable-cell nudged elastic band driver.

Engines: Quantum ESPRESSO (pw.x), CP2K, ABACUS, plus built-in EMT/LJ for
quick tests.  One command does everything: relax endpoints -> interpolate ->
VC-NEB -> write trajectory / barrier plot / summary.
"""

__version__ = "0.1.0"

from .config import Config, load_config, write_default_config  # noqa: F401
from .engines import make_engine, ENGINE_NAMES  # noqa: F401
from .band import Band, plot_band  # noqa: F401

__all__ = ["Config", "load_config", "write_default_config", "make_engine",
           "ENGINE_NAMES", "Band", "plot_band", "__version__"]
