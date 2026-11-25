from __future__ import annotations

from collections import deque
import math
from typing import Mapping, Sequence

import networkx as nx

from gasnet.exceptions import DataContractError, TopologyError
from gasnet.models import FluidProperties, PipeRecord, SolverResult
from gasnet.physics.resistance import pipe_pressure_drop_pa, pipe_resistance


def solve_tree_network(
    pipes: Sequence[PipeRecord],
    source_node: str,
    source_pressure_gauge_pa: float,
    withdrawal_nm3_h: Mapping[str, float],
    fluid: FluidProperties,
) -> SolverResult:
    """Exact O(N) baseline for a connected tree with prescribed withdrawals."""
    if not pipes:
        raise TopologyError("tree solver requires at least one pipe")
    if not math.isfinite(source_pressure_gauge_pa):
        raise DataContractError("source pressure must be finite")
    pipe_ids = [pipe.pipe_id for pipe in pipes]
    if len(pipe_ids) != len(set(pipe_ids)):
        raise DataContractError("pipe_id must be unique")
    graph = nx.Graph()
    pipe_by_pair: dict[frozenset[str], PipeRecord] = {}
    for pipe in pipes:
        if pipe.from_node == pipe.to_node:
            raise TopologyError(f"self-loop is not allowed for pipe {pipe.pipe_id}")
        pair = frozenset((pipe.from_node, pipe.to_node))
        if pair in pipe_by_pair:
            raise TopologyError(
                f"parallel pipes {pipe.from_node}-{pipe.to_node} require a MultiGraph solver"
            )
        graph.add_edge(pipe.from_node, pipe.to_node, pipe_id=pipe.pipe_id)
        pipe_by_pair[pair] = pipe
    if source_node not in graph:
        raise TopologyError(f"source node {source_node} is absent")
    if not nx.is_tree(graph):
        raise TopologyError("tree solver requires one connected acyclic graph")
    unknown_demand_nodes = set(withdrawal_nm3_h).difference(graph.nodes)
    if unknown_demand_nodes:
        raise DataContractError(
            f"withdrawals reference unknown nodes: {sorted(unknown_demand_nodes)}"
        )
    for node, value in withdrawal_nm3_h.items():
        if not math.isfinite(value) or value < 0:
            raise DataContractError(
                f"withdrawal at {node} must be finite and non-negative"
            )

    parent: dict[str, str | None] = {source_node: None}
    order: list[str] = []
    queue = deque([source_node])
    while queue:
        node = queue.popleft()
        order.append(node)
        for neighbor in graph.neighbors(node):
            if neighbor not in parent:
                parent[neighbor] = node
                queue.append(neighbor)

    subtree_flow = {node: float(withdrawal_nm3_h.get(node, 0.0)) for node in graph.nodes}
    for node in reversed(order[1:]):
        upstream = parent[node]
        assert upstream is not None
        subtree_flow[upstream] += subtree_flow[node]

    pressure = {source_node: float(source_pressure_gauge_pa)}
    pipe_flow: dict[str, float] = {}
    for node in order[1:]:
        upstream = parent[node]
        assert upstream is not None
        pipe = pipe_by_pair[frozenset((upstream, node))]
        flow = subtree_flow[node]
        pressure[node] = pressure[upstream] - pipe_pressure_drop_pa(
            flow, pipe_resistance(pipe, fluid)
        )
        pipe_flow[pipe.pipe_id] = (
            flow if pipe.from_node == upstream and pipe.to_node == node else -flow
        )

    residual = _mass_balance_residual(pipes, pipe_flow, withdrawal_nm3_h)
    return SolverResult(
        node_pressure_gauge_pa=pressure,
        pipe_flow_nm3_h=pipe_flow,
        mass_balance_residual_nm3_h=residual,
        converged=True,
        iterations=0,
        metadata={"method": "tree_analytic", "source_flow_nm3_h": subtree_flow[source_node]},
    )


def _mass_balance_residual(
    pipes: Sequence[PipeRecord],
    pipe_flow: Mapping[str, float],
    withdrawal: Mapping[str, float],
) -> dict[str, float]:
    nodes = {pipe.from_node for pipe in pipes} | {pipe.to_node for pipe in pipes}
    net_outflow = {node: 0.0 for node in nodes}
    for pipe in pipes:
        signed_flow = pipe_flow[pipe.pipe_id]
        net_outflow[pipe.from_node] += signed_flow
        net_outflow[pipe.to_node] -= signed_flow
    return {
        node: net_outflow[node] + float(withdrawal.get(node, 0.0))
        for node in nodes
    }
