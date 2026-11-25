class GasNetError(Exception):
    """Base exception for the project."""


class DataContractError(GasNetError):
    """Input data violate a declared schema or unit contract."""


class TopologyError(GasNetError):
    """The network topology cannot support the requested calculation."""


class SolverConvergenceError(GasNetError):
    """A nonlinear solver failed to satisfy its convergence contract."""

