# Step6 Development Utilities

This directory stores synthetic-data and test helpers for the transmission module.

## Contents

- `generate_synthetic_data.py`: builds synthetic case and covariate tables for development and smoke testing.
- `test_transmission_dynamics.py`: pytest-based checks for the step6 transmission scripts.

## Rule

Keep production entry points in `../scripts/` and keep development-only helpers here.
