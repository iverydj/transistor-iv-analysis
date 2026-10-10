import warnings
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

from transistor_iv import Device, read_clarius_transfer  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "examples" / "data"
BASELINES = Path(__file__).resolve().parent / "data"


@pytest.fixture(params=["n", "p"])
def device(request):
    """W/L = 100/10 um on 100 nm SiO2, both polarities."""
    return Device.from_dielectric(100, 10, 3.9, 100, request.param)


def read_quiet(name: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return read_clarius_transfer(DATA / name)
