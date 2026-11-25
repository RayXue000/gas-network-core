from gasnet.models import RegulatorRecord
from gasnet.physics.regulator import regulator_capacity_ratio, regulator_outlet_pressure


def test_regulator_droop_above_rated_flow() -> None:
    regulator = RegulatorRecord(
        regulator_id="REG1",
        node_id="B",
        setpoint_pressure_gauge_pa=2800.0,
        rated_flow_nm3_h=80.0,
        droop_pa_per_nm3_h=10.0,
        minimum_outlet_pressure_gauge_pa=1500.0,
    )
    assert regulator_outlet_pressure(regulator, 60.0) == 2800.0
    assert regulator_outlet_pressure(regulator, 100.0) == 2600.0
    assert regulator_capacity_ratio(regulator, 160.0) == 2.0


def test_regulator_rejects_unverified_overload_by_default() -> None:
    regulator = RegulatorRecord(
        regulator_id="REG1",
        node_id="B",
        setpoint_pressure_gauge_pa=2800.0,
        rated_flow_nm3_h=80.0,
    )
    with pytest.raises(DataContractError, match="overloaded"):
        regulator_outlet_pressure(regulator, 160.0)
    assert (
        regulator_outlet_pressure(
            regulator,
            160.0,
            allow_unverified_overload=True,
        )
        == 2800.0
    )
import pytest

from gasnet.exceptions import DataContractError
