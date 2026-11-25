"""Synthetic three-node example; no field measurements or site geometry."""
from gasnet.models import FluidProperties, PipeRecord
from gasnet.solvers.tree import solve_tree_network
from gasnet.solvers.newton import solve_fixed_pressure_network

pipes = [PipeRecord("p1", "source", "junction", 15.0, 60.0),
         PipeRecord("p2", "junction", "load", 10.0, 45.0)]
fluid = FluidProperties()
tree = solve_tree_network(pipes, "source", 2800.0, {"load": 5.0}, fluid)
newton = solve_fixed_pressure_network(pipes, {"source": 2800.0}, {"load": 5.0}, fluid)
assert tree.converged and newton.converged
assert abs(tree.node_pressure_gauge_pa["load"] - newton.node_pressure_gauge_pa["load"]) < 1e-5
print(dict(tree.node_pressure_gauge_pa))
