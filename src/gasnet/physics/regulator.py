from __future__ import annotations

import math

from gasnet.exceptions import DataContractError
from gasnet.models import RegulatorRecord


def regulator_outlet_pressure(
    regulator: RegulatorRecord,
    flow_nm3_h: float,
    *,
    allow_unverified_overload: bool = False,
) -> float:
    """Evaluate a documented baseline regulator curve.

    Below rated flow the outlet holds the setpoint. Above rated flow a linear
    droop is applied when a measured or manufacturer-supplied slope is present.
    Missing droop data are acceptable below rated flow. Above rated flow the
    default is strict: an outlet pressure cannot be inferred without a curve.
    Exploratory callers may opt into the old fixed-setpoint assumption explicitly.
    """
    _validate_regulator(regulator)
    if not math.isfinite(flow_nm3_h):
        raise DataContractError("regulator flow must be finite")
    overload = max(abs(flow_nm3_h) - regulator.rated_flow_nm3_h, 0.0)
    if overload > 0 and regulator.droop_pa_per_nm3_h is None and not allow_unverified_overload:
        raise DataContractError(
            "regulator is overloaded but no droop curve is available; "
            "supply a droop slope or explicitly allow the unverified assumption"
        )
    droop = 0.0 if regulator.droop_pa_per_nm3_h is None else overload * regulator.droop_pa_per_nm3_h
    return max(
        regulator.setpoint_pressure_gauge_pa - droop,
        regulator.minimum_outlet_pressure_gauge_pa,
    )


def regulator_capacity_ratio(regulator: RegulatorRecord, flow_nm3_h: float) -> float:
    _validate_regulator(regulator)
    if not math.isfinite(flow_nm3_h):
        raise DataContractError("regulator flow must be finite")
    return abs(flow_nm3_h) / regulator.rated_flow_nm3_h


def _validate_regulator(regulator: RegulatorRecord) -> None:
    values = {
        "setpoint pressure": regulator.setpoint_pressure_gauge_pa,
        "rated flow": regulator.rated_flow_nm3_h,
        "minimum outlet pressure": regulator.minimum_outlet_pressure_gauge_pa,
    }
    for name, value in values.items():
        if not math.isfinite(value):
            raise DataContractError(f"regulator {name} must be finite")
    if regulator.rated_flow_nm3_h <= 0:
        raise DataContractError("regulator rated flow must be positive")
    if regulator.minimum_outlet_pressure_gauge_pa > regulator.setpoint_pressure_gauge_pa:
        raise DataContractError("minimum outlet pressure cannot exceed the setpoint")
    if regulator.droop_pa_per_nm3_h is not None:
        if not math.isfinite(regulator.droop_pa_per_nm3_h) or regulator.droop_pa_per_nm3_h < 0:
            raise DataContractError("regulator droop slope must be finite and non-negative")
