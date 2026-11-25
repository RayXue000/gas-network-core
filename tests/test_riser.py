from gasnet.models import FluidProperties
from gasnet.physics.riser import solve_riser_profile


def test_riser_flow_decreases_after_each_floor() -> None:
    profile = solve_riser_profile(
        base_pressure_gauge_pa=2500.0,
        floor_withdrawals_nm3_h=[1.0, 2.0, 3.0],
        inner_diameter_mm=40.0,
        floor_height_m=3.0,
        fluid=FluidProperties(),
    )
    assert [row.segment_flow_nm3_h for row in profile] == [6.0, 5.0, 3.0]
    assert [row.floor for row in profile] == [1, 2, 3]

