"""Independent CPD geometry and fixed-equilibrium optimization support."""
# Match the existing repository's isolated scientific dependency discovery.
from pathlib import Path
import sys

_deps = Path(__file__).resolve().parents[1] / '.local_dependencies/d3d-fixed-equilibrium'
if _deps.is_dir() and str(_deps) not in sys.path:
    sys.path.insert(0, str(_deps))
