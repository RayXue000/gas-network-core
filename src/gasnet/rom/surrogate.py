from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SurrogateResult:
    model: object
    x_scaler: object
    y_scaler: object | None
    training_indices: np.ndarray


def fit_surrogate(features, targets, training_indices, model_kind: str = "gpr") -> SurrogateResult:
    """Fit an optional surrogate using an explicit, externally defined split.

    The caller must create time- or scenario-based indices. This function does
    not perform random splitting, preventing accidental leakage by default.
    """
    try:
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel
        from sklearn.preprocessing import StandardScaler
    except ImportError as exc:
        raise RuntimeError("Install the optional 'rom' dependencies") from exc

    x = np.asarray(features, dtype=float)
    y = np.asarray(targets, dtype=float)
    indices = np.asarray(training_indices, dtype=int)
    scaler = StandardScaler().fit(x[indices])
    x_train = scaler.transform(x[indices])

    if model_kind != "gpr":
        raise ValueError("the initial refactor supports model_kind='gpr' only")
    kernel = ConstantKernel(1.0, (1.0e-3, 1.0e3)) * RBF(
        length_scale=np.ones(x.shape[1]), length_scale_bounds=(1.0e-2, 1.0e2)
    ) + WhiteKernel(noise_level=1.0e-5)
    model = GaussianProcessRegressor(kernel=kernel, normalize_y=True, random_state=42)
    model.fit(x_train, y[indices])
    return SurrogateResult(model=model, x_scaler=scaler, y_scaler=None, training_indices=indices)

