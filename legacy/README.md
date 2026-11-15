# Legacy Scripts Archive

The scripts in this directory are preserved for archival, baseline reference, and historical auditing:

1. `tcr_simulation_ssot.py` — Pre-refactoring single source of truth simulation script.
2. `tcr_advanced.py` — Pre-refactoring multi-regime simulation and NSGA-II optimization script.
3. `thermally_coupled_reactor.py` — Original initial thermally coupled reactor script.

### Notice
These scripts have been completely superseded by the unified, tested, and modular `tcr/` package, configured via `config.yaml`, and executed exclusively through the single entry point `run_all.py`.

Please do not execute or modify these legacy scripts directly. All active simulations, parameter definitions, and tests are managed in `tcr/`, `config.yaml`, and `tests/`.
