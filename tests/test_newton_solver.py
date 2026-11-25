import pytest

from gasnet.exceptions import DataContractError, TopologyError
from gasnet.models import FluidProperties, PipeRecord
from gasnet.solvers.newton import solve_fixed_pressure_network


def test_newton_solver_matches_one_pipe_demand() -> None:
    result = solve_fixed_pressure_network(
        [PipeRecord("E1", "B", "R1", 10.0, 73.6, friction_factor=0.015)],
        fixed_pressure_gauge_pa={"B": 2800.0},
        withdrawal_nm3_h={"R1": 10.0},
        fluid=FluidProperties(),
    )
    assert result.pipe_flow_nm3_h["E1"] == pytest.approx(10.0)
    assert abs(result.mass_balance_residual_nm3_h["R1"]) < 1.0e-8


def test_newton_solver_rejects_unknown_demand_and_unbounded_component() -> None:
    pipe = PipeRecord("E1", "B", "R1", 10.0, 73.6)
    with pytest.raises(DataContractError, match="unknown nodes"):
        solve_fixed_pressure_network(
            [pipe],
            fixed_pressure_gauge_pa={"B": 2800.0},
            withdrawal_nm3_h={"TYPO": 1.0},
            fluid=FluidProperties(),
        )

    disconnected = [
        pipe,
        PipeRecord("E2", "X", "Y", 10.0, 73.6),
    ]
    with pytest.raises(TopologyError, match="no fixed-pressure boundary"):
        solve_fixed_pressure_network(
            disconnected,
            fixed_pressure_gauge_pa={"B": 2800.0},
            withdrawal_nm3_h={"R1": 1.0, "Y": 1.0},
            fluid=FluidProperties(),
        )
