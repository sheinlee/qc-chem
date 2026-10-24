import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pytest  # noqa: E402

from qcchem.integrals import build_mol, make_active_space  # noqa: E402


@pytest.fixture(scope="session")
def h2():
    return make_active_space("H2", build_mol("H 0 0 0; H 0 0 0.74", "sto-3g"), ncas=2, nelecas=2)


@pytest.fixture(scope="session")
def lih44():
    return make_active_space("LiH44", build_mol("Li 0 0 0; H 0 0 1.6", "sto-3g"), ncas=4, nelecas=4)


@pytest.fixture(scope="session")
def n2_stretched():
    # stretched N2: multi-reference regime, good test of sector-restricted ED vs CASCI
    return make_active_space("N2s", build_mol("N 0 0 0; N 0 0 2.0", "sto-3g"), ncas=6, nelecas=6)
