from __future__ import annotations

import math

from gasnet.exceptions import DataContractError
from gasnet.models import FluidProperties, PipeRecord


def equivalent_length_m(pipe: PipeRecord, friction_factor: float) -> float:
    if not math.isfinite(friction_factor) or friction_factor <= 0:
        raise DataContractError("friction factor must be positive")
    if not math.isfinite(pipe.local_loss_coefficient) or pipe.local_loss_coefficient < 0:
        raise DataContractError("local loss coefficient must be finite and non-negative")
    local_length_m = (
        pipe.local_loss_coefficient * (pipe.inner_diameter_mm / 1000.0) / friction_factor
    )
    return pipe.length_m + local_length_m


def pipe_resistance(pipe: PipeRecord, fluid: FluidProperties) -> float:
    """Return R in delta_p[Pa] = R * Q[Nm3/h] * abs(Q[Nm3/h])."""
    friction_factor = (
        fluid.friction_factor_default if pipe.friction_factor is None else pipe.friction_factor
    )
    if not math.isfinite(pipe.length_m) or not math.isfinite(pipe.inner_diameter_mm):
        raise DataContractError("pipe length and inner diameter must be finite")
    if pipe.length_m <= 0 or pipe.inner_diameter_mm <= 0:
        raise DataContractError("pipe length and inner diameter must be positive")
    positive_fluid_values = {
        "standard temperature": fluid.standard_temperature_k,
        "operating temperature": fluid.operating_temperature_k,
        "standard density": fluid.standard_density_kg_m3,
        "GB 50028 constant": fluid.gb50028_constant,
    }
    for name, value in positive_fluid_values.items():
        if not math.isfinite(value) or value <= 0:
            raise DataContractError(f"{name} must be finite and positive")
    length_m = equivalent_length_m(pipe, friction_factor)
    resistance = (
        fluid.gb50028_constant
        * friction_factor
        * fluid.standard_density_kg_m3
        * fluid.operating_temperature_k
        * length_m
        / (fluid.standard_temperature_k * pipe.inner_diameter_mm**5)
    )
    if not math.isfinite(resistance) or resistance <= 0:
        raise DataContractError("computed pipe resistance must be finite and positive")
    return resistance


def pipe_pressure_drop_pa(flow_nm3_h: float, resistance: float) -> float:
    if not math.isfinite(flow_nm3_h):
        raise DataContractError("flow must be finite")
    if not math.isfinite(resistance) or resistance <= 0:
        raise DataContractError("resistance must be finite and positive")
    return resistance * flow_nm3_h * abs(flow_nm3_h)


def flow_from_pressure_drop_nm3_h(delta_pressure_pa: float, resistance: float) -> float:
    if not math.isfinite(delta_pressure_pa):
        raise DataContractError("pressure difference must be finite")
    if not math.isfinite(resistance) or resistance <= 0:
        raise DataContractError("resistance must be finite and positive")
    if delta_pressure_pa == 0:
        return 0.0
    return math.copysign(math.sqrt(abs(delta_pressure_pa) / resistance), delta_pressure_pa)


def flow_derivative_wrt_pressure(
    delta_pressure_pa: float,
    resistance: float,
    minimum_sqrt_pressure: float = 1.0e-6,
) -> float:
    if not math.isfinite(delta_pressure_pa):
        raise DataContractError("pressure difference must be finite")
    if not math.isfinite(resistance) or resistance <= 0:
        raise DataContractError("resistance must be finite and positive")
    if not math.isfinite(minimum_sqrt_pressure) or minimum_sqrt_pressure <= 0:
        raise DataContractError("minimum square-root pressure must be finite and positive")
    denominator = max(math.sqrt(abs(delta_pressure_pa)), minimum_sqrt_pressure)
    return 1.0 / (2.0 * math.sqrt(resistance) * denominator)
