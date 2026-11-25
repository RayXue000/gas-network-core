import pytest

from gasnet.exceptions import DataContractError
from gasnet.models import FluidProperties, PipeRecord
from gasnet.physics.resistance import pipe_resistance


def test_explicit_zero_friction_factor_is_not_replaced_by_default() -> None:
    pipe = PipeRecord(
        "E1",
        "B",
        "R1",
        10.0,
        73.6,
        friction_factor=0.0,
    )
    with pytest.raises(DataContractError, match="friction factor"):
        pipe_resistance(pipe, FluidProperties())
