# Gas Network Core

Selected numerical source from Rui Xue's gas-network modelling project.
管网数值建模核心代码节选：水力求解、物理部件与降阶建模。

- `physics/`: pipe resistance, regulator response and riser pressure.
- `solvers/tree.py`: prescribed-demand tree solver with continuity checks.
- `solvers/newton.py`: fixed-pressure network solver.
- `rom/`: POD projection/reconstruction and an optional GPR surrogate.

## Run

Python 3.10+:
```sh
python -m pip install -e ".[test,rom]"
python example.py
python -m pytest -q
```

## Scope

This is a source-only subset of the local `gasnet` refactor. The numerical
modules and selected original unit tests are retained; the small example is
new and synthetic. No customer topology, household records, calibration
tables, SCADA data, private reports, UI or complete application is included.
Results of the example are numerical checks, not field validation or an
engineering design recommendation. Input units are indicated in the code.

[Academic homepage](https://rayxue000.github.io/)
· [Bilingual interactive demonstration](https://rayxue000.github.io/demos/gas-twin/index.html)
