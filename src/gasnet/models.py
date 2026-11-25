from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping


@dataclass(frozen=True)
class NodeRecord:
    node_id: str
    node_type: str
    elevation_m: float = 0.0
    x: float | None = None
    y: float | None = None
    source_ref: str | None = None


@dataclass(frozen=True)
class PipeRecord:
    pipe_id: str
    from_node: str
    to_node: str
    length_m: float
    inner_diameter_mm: float
    roughness_mm: float | None = None
    friction_factor: float | None = None
    local_loss_coefficient: float = 0.0
    material: str | None = None
    source_ref: str | None = None


@dataclass(frozen=True)
class DemandRecord:
    demand_id: str
    node_id: str
    flow_nm3_h: float
    building_id: str | None = None
    riser_id: str | None = None
    floor: int | None = None
    household_id: str | None = None
    scenario: str = "baseline"
    source_ref: str | None = None


@dataclass(frozen=True)
class RegulatorRecord:
    regulator_id: str
    node_id: str
    setpoint_pressure_gauge_pa: float
    rated_flow_nm3_h: float
    inlet_pressure_gauge_pa: float | None = None
    droop_pa_per_nm3_h: float | None = None
    minimum_outlet_pressure_gauge_pa: float = 0.0
    source_ref: str | None = None


@dataclass(frozen=True)
class FluidProperties:
    standard_temperature_k: float = 273.15
    operating_temperature_k: float = 293.15
    standard_density_kg_m3: float = 0.7174
    friction_factor_default: float = 0.015
    gb50028_constant: float = 6.26e7


@dataclass(frozen=True)
class SolverResult:
    node_pressure_gauge_pa: Mapping[str, float]
    pipe_flow_nm3_h: Mapping[str, float]
    mass_balance_residual_nm3_h: Mapping[str, float]
    converged: bool
    iterations: int
    metadata: Mapping[str, object] = field(default_factory=dict)

