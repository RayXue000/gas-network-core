import pytest

from gasnet.exceptions import DataContractError
from gasnet.models import FluidProperties, PipeRecord
from gasnet.physics.resistance import pipe_pressure_drop_pa, pipe_resistance
from gasnet.solvers.tree import solve_tree_network


def test_tree_solver_matches_manual_flow_aggregation() -> None:
    fluid = FluidProperties()
    pipes = [
        PipeRecord("E1", "B", "T1", 10.0, 73.6, friction_factor=0.015),
        PipeRecord("E2", "T1", "R1", 5.0, 51.4, friction_factor=0.015),
    ]
    result = solve_tree_network(
        pipes,
        source_node="B",
        source_pressure_gauge_pa=2800.0,
        withdrawal_nm3_h={"R1": 10.0},
        fluid=fluid,
    )
    expected_t1 = 2800.0 - pipe_pressure_drop_pa(10.0, pipe_resistance(pipes[0], fluid))
    expected_r1 = expected_t1 - pipe_pressure_drop_pa(10.0, pipe_resistance(pipes[1], fluid))
    assert result.pipe_flow_nm3_h == {"E1": 10.0, "E2": 10.0}
    assert result.node_pressure_gauge_pa["T1"] == pytest.approx(expected_t1)
    assert result.node_pressure_gauge_pa["R1"] == pytest.approx(expected_r1)
    assert abs(result.mass_balance_residual_nm3_h["T1"]) < 1.0e-12
    assert abs(result.mass_balance_residual_nm3_h["R1"]) < 1.0e-12


def test_tree_solver_reports_flow_in_declared_pipe_direction() -> None:
    pipe = PipeRecord("E1", "R1", "B", 10.0, 73.6, friction_factor=0.015)
    result = solve_tree_network(
        [pipe],
        source_node="B",
        source_pressure_gauge_pa=2800.0,
        withdrawal_nm3_h={"R1": 10.0},
        fluid=FluidProperties(),
    )
    assert result.pipe_flow_nm3_h["E1"] == -10.0
    assert abs(result.mass_balance_residual_nm3_h["R1"]) < 1.0e-12


def test_tree_solver_rejects_unknown_demand_node() -> None:
    pipe = PipeRecord("E1", "B", "R1", 10.0, 73.6)
    with pytest.raises(DataContractError, match="unknown nodes"):
        solve_tree_network(
            [pipe],
            source_node="B",
            source_pressure_gauge_pa=2800.0,
            withdrawal_nm3_h={"TYPO": 10.0},
            fluid=FluidProperties(),
        )
