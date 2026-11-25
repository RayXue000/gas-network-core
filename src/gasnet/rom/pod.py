from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PODModel:
    mean_state: np.ndarray
    modes: np.ndarray
    singular_values: np.ndarray
    cumulative_energy: np.ndarray
    retained_rank: int


def fit_pod(snapshot_matrix, target_energy_fraction: float = 0.999) -> PODModel:
    snapshots = np.asarray(snapshot_matrix, dtype=float)
    if snapshots.ndim != 2:
        raise ValueError("snapshot matrix must have shape [state, snapshot]")
    if not 0 < target_energy_fraction <= 1:
        raise ValueError("target energy fraction must be in (0, 1]")
    mean_state = snapshots.mean(axis=1, keepdims=True)
    centered = snapshots - mean_state
    modes, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
    energy = singular_values**2
    if float(energy.sum()) == 0.0:
        raise ValueError("POD requires at least one non-constant snapshot direction")
    cumulative = np.cumsum(energy) / energy.sum()
    rank = int(np.searchsorted(cumulative, target_energy_fraction) + 1)
    return PODModel(mean_state, modes, singular_values, cumulative, rank)


def project_pod(snapshot_matrix, model: PODModel) -> np.ndarray:
    centered = np.asarray(snapshot_matrix, dtype=float) - model.mean_state
    return model.modes[:, : model.retained_rank].T @ centered


def reconstruct_pod(coefficients, model: PODModel) -> np.ndarray:
    coefficients = np.asarray(coefficients, dtype=float)
    return model.mean_state + model.modes[:, : model.retained_rank] @ coefficients
