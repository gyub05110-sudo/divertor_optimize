"""Local fixed-EFIT DIII-D optimization; independent of legacy FreeGS solvers."""
from pathlib import Path
import sys
_deps = Path(__file__).resolve().parents[1]/".local_dependencies/d3d-fixed-equilibrium"
if _deps.is_dir():
    sys.path.insert(0,str(_deps))
