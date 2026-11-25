from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from gasnet.exceptions import DataContractError
from gasnet.models import FluidProperties


@dataclass(frozen=True)
class RiserFloorResult:
    floor: int
    pressure_gauge_pa: float
    segment_flow_nm3_h: float
    friction_drop_pa: float
    buoyancy_gain_pa: float


def solve_riser_profile(
    base_pressure_gauge_pa: float,
    floor_withdrawals_nm3_h: Sequence[float],
    inner_diameter_mm: float,
    floor_height_m: float,
    fluid: FluidProperties,
    air_density_kg_m3: float = 1.293,
    gravity_m_s2: float = 9.81,
) -> list[RiserFloorResult]:
    """Solve one riser using one segment from the base to every floor.

    `floor_withdrawals_nm3_h[0]` is the first-floor withdrawal. The first segment
    carries all floor withdrawals; after each floor the local demand is removed.
    This convention avoids the previous ambiguity where floor 1 sometimes had
    zero vertical height and sometimes had one full segment.
    """
    if inner_diameter_mm <= 0 or floor_height_m <= 0:
        raise DataContractError("riser diameter and floor height must be positive")
    if any(flow < 0 for flow in floor_withdrawals_nm3_h):
        raise DataContractError("floor withdrawals must be non-negative")

    friction_factor = fluid.friction_factor_default
    segment_resistance = (
        fluid.gb50028_constant
        * friction_factor
        * fluid.standard_density_kg_m3
        * fluid.operating_temperature_k
        * floor_height_m
        / (fluid.standard_temperature_k * inner_diameter_mm**5)
    )
    buoyancy_gain = (
        gravity_m_s2
        * (air_density_kg_m3 - fluid.standard_density_kg_m3)
        * floor_height_m
    )

    pressure = float(base_pressure_gauge_pa)
    remaining_flow = float(sum(floor_withdrawals_nm3_h))
    results: list[RiserFloorResult] = []
    for floor, withdrawal in enumerate(floor_withdrawals_nm3_h, start=1):
        friction_drop = segment_resistance * remaining_flow * abs(remaining_flow)
        pressure = pressure - friction_drop + buoyancy_gain
        results.append(
            RiserFloorResult(
                floor=floor,
                pressure_gauge_pa=pressure,
                segment_flow_nm3_h=remaining_flow,
                friction_drop_pa=friction_drop,
                buoyancy_gain_pa=buoyancy_gain,
            )
        )
        remaining_flow -= float(withdrawal)
    return results

