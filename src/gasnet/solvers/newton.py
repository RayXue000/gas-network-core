from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping, Sequence

import networkx as nx
import numpy as np

from gasnet.exceptions import DataContractError, SolverConvergenceError, TopologyError
from gasnet.models import FluidProperties, PipeRecord, SolverResult
from gasnet.physics.resistance import (
    flow_derivative_wrt_pressure,
    flow_from_pressure_drop_nm3_h,
    pipe_resistance,
)


@dataclass(frozen=True)
class NewtonOptions:
    max_iterations: int = 50
    pressure_tolerance_pa: float = 1.0e-6
    residual_tolerance_nm3_h: float = 1.0e-8
    damping: float = 1.0
    minimum_sqrt_pressure: float = 1.0e-6


def solve_fixed_pressure_network(
    pipes: Sequence[PipeRecord],
    fixed_pressure_gauge_pa: Mapping[str, float],
    withdrawal_nm3_h: Mapping[str, float],
    fluid: FluidProperties,
    options: NewtonOptions = NewtonOptions(),
) -> SolverResult:
    """Solve a general steady network with prescribed nodal withdrawals.

    Regulator flow-dependent boundaries should be resolved by an outer boundary
    loop; this solver intentionally accepts only explicit fixed pressures.
    """
    if not pipes:
        raise TopologyError("network solver requires at least one pipe")
    if not fixed_pressure_gauge_pa:
        raise TopologyError("at least one fixed-pressure node is required")
    _validate_options(options)
    pipe_ids = [pipe.pipe_id for pipe in pipes]
    if len(pipe_ids) != len(set(pipe_ids)):
        raise DataContractError("pipe_id must be unique")
    if any(pipe.from_node == pipe.to_node for pipe in pipes):
        raise TopologyError("self-loop pipes are not supported")

    nodes = sorted({p.from_node for p in pipes} | {p.to_node for p in pipes})
    node_set = set(nodes)
    unknown_fixed_nodes = set(fixed_pressure_gauge_pa).difference(node_set)
    if unknown_fixed_nodes:
        raise DataContractError(
            f"fixed pressures reference unknown nodes: {sorted(unknown_fixed_nodes)}"
        )
    unknown_demand_nodes = set(withdrawal_nm3_h).difference(node_set)
    if unknown_demand_nodes:
        raise DataContractError(
            f"withdrawals reference unknown nodes: {sorted(unknown_demand_nodes)}"
        )
    for node, value in fixed_pressure_gauge_pa.items():
        if not math.isfinite(value):
            raise DataContractError(f"fixed pressure at {node} must be finite")
    for node, value in withdrawal_nm3_h.items():
        if not math.isfinite(value) or value < 0:
            raise DataContractError(
                f"withdrawal at {node} must be finite and non-negative"
            )

    graph = nx.Graph((pipe.from_node, pipe.to_node) for pipe in pipes)
    for component in nx.connected_components(graph):
        if component.isdisjoint(fixed_pressure_gauge_pa):
            raise TopologyError(
                f"network component has no fixed-pressure boundary: {sorted(component)}"
            )

    unknown_nodes = [node for node in nodes if node not in fixed_pressure_gauge_pa]
    unknown_index = {node: idx for idx, node in enumerate(unknown_nodes)}
    incident: dict[str, list[tuple[PipeRecord, str]]] = {node: [] for node in nodes}
    resistance = {pipe.pipe_id: pipe_resistance(pipe, fluid) for pipe in pipes}
    for pipe in pipes:
        incident[pipe.from_node].append((pipe, pipe.to_node))
        incident[pipe.to_node].append((pipe, pipe.from_node))

    initial_pressure = float(np.mean(list(fixed_pressure_gauge_pa.values())))
    pressure_vector = np.full(len(unknown_nodes), initial_pressure, dtype=float)

    def pressure(node: str) -> float:
        if node in fixed_pressure_gauge_pa:
            return float(fixed_pressure_gauge_pa[node])
        return float(pressure_vector[unknown_index[node]])

    iteration = 0
    maximum_unknown_residual = 0.0
    last_pressure_step = math.inf
    converged = not unknown_nodes
    for iteration in range(1, options.max_iterations + 1) if unknown_nodes else ():
        residual = np.zeros(len(unknown_nodes), dtype=float)
        jacobian = np.zeros((len(unknown_nodes), len(unknown_nodes)), dtype=float)
        for node in unknown_nodes:
            row = unknown_index[node]
            net_outflow = 0.0
            for pipe, neighbor in incident[node]:
                delta_pressure = pressure(node) - pressure(neighbor)
                edge_resistance = resistance[pipe.pipe_id]
                net_outflow += flow_from_pressure_drop_nm3_h(delta_pressure, edge_resistance)
                derivative = flow_derivative_wrt_pressure(
                    delta_pressure,
                    edge_resistance,
                    options.minimum_sqrt_pressure,
                )
                jacobian[row, row] += derivative
                if neighbor in unknown_index:
                    jacobian[row, unknown_index[neighbor]] -= derivative
            residual[row] = net_outflow + float(withdrawal_nm3_h.get(node, 0.0))

        maximum_unknown_residual = float(np.max(np.abs(residual)))
        if maximum_unknown_residual <= options.residual_tolerance_nm3_h and (
            iteration == 1 or last_pressure_step <= options.pressure_tolerance_pa
        ):
            converged = True
            break
        try:
            delta = np.linalg.solve(jacobian, -residual)
        except np.linalg.LinAlgError as exc:
            raise SolverConvergenceError("singular Newton Jacobian") from exc
        damped_delta = options.damping * delta
        pressure_vector += damped_delta
        last_pressure_step = float(np.max(np.abs(damped_delta)))
    if not converged:
        raise SolverConvergenceError("Newton solver exceeded the iteration limit")

    node_pressure = {node: pressure(node) for node in nodes}
    pipe_flow: dict[str, float] = {}
    for pipe in pipes:
        pipe_flow[pipe.pipe_id] = flow_from_pressure_drop_nm3_h(
            node_pressure[pipe.from_node] - node_pressure[pipe.to_node],
            resistance[pipe.pipe_id],
        )
    mass_residual = _node_residuals(pipes, pipe_flow, withdrawal_nm3_h)
    return SolverResult(
        node_pressure_gauge_pa=node_pressure,
        pipe_flow_nm3_h=pipe_flow,
        mass_balance_residual_nm3_h=mass_residual,
        converged=True,
        iterations=iteration,
        metadata={
            "method": "newton_fixed_pressure",
            "maximum_unknown_residual_nm3_h": maximum_unknown_residual,
        },
    )


def _validate_options(options: NewtonOptions) -> None:
    if options.max_iterations <= 0:
        raise DataContractError("max_iterations must be positive")
    positive_values = {
        "pressure tolerance": options.pressure_tolerance_pa,
        "residual tolerance": options.residual_tolerance_nm3_h,
        "minimum square-root pressure": options.minimum_sqrt_pressure,
    }
    for name, value in positive_values.items():
        if not math.isfinite(value) or value <= 0:
            raise DataContractError(f"{name} must be finite and positive")
    if not math.isfinite(options.damping) or not 0 < options.damping <= 1:
        raise DataContractError("Newton damping must be in (0, 1]")


def _node_residuals(
    pipes: Sequence[PipeRecord],
    pipe_flow: Mapping[str, float],
    withdrawal: Mapping[str, float],
) -> dict[str, float]:
    nodes = {p.from_node for p in pipes} | {p.to_node for p in pipes}
    net_outflow = {node: 0.0 for node in nodes}
    for pipe in pipes:
        flow = pipe_flow[pipe.pipe_id]
        net_outflow[pipe.from_node] += flow
        net_outflow[pipe.to_node] -= flow
    return {
        node: net_outflow[node] + float(withdrawal.get(node, 0.0))
        for node in nodes
    }
