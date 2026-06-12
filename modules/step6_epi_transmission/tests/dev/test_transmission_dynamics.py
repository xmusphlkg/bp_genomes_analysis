#!/usr/bin/env python3
"""Test suite for transmission dynamics module."""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

STEP6_DIR = Path(__file__).resolve().parents[2]
SYNTHETIC_DATA_SCRIPT = Path(__file__).resolve().parent / "generate_synthetic_data.py"
RE_ESTIMATION_SCRIPT = STEP6_DIR / "bin" / "step6_06_estimate_reproduction_numbers_v2.py"
TRANSMISSION_MODEL_SCRIPT = STEP6_DIR / "bin" / "step6_07_fit_transmission_models.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def run_command(cmd: list) -> bool:
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            print(f"✗ Command failed: {' '.join(cmd)}")
            print(f"STDERR: {result.stderr}")
            return False
        print(f"✓ {' '.join(cmd[:3])}... OK")
        return True
    except subprocess.TimeoutExpired:
        print(f"✗ Timeout: {' '.join(cmd)}")
        return False
    except Exception as e:
        print(f"✗ Error: {e}")
        return False


def generate_synthetic_test_data(output_dir: Path) -> None:
    cmd = [
        sys.executable,
        str(SYNTHETIC_DATA_SCRIPT),
        "--n-countries",
        "5",
        "--n-years",
        "10",
        "-o",
        str(output_dir),
    ]
    assert run_command(cmd)


def run_re_estimation(output_dir: Path) -> None:
    generate_synthetic_test_data(output_dir)
    cmd = [
        sys.executable,
        str(RE_ESTIMATION_SCRIPT),
        "-i",
        str(output_dir / "synthetic_cases.tsv"),
        "-o",
        str(output_dir / "re_output"),
        "--allow-annual-disaggregation",
    ]
    assert run_command(cmd)


def test_generation_interval():
    module = load_module(RE_ESTIMATION_SCRIPT, "step6_re_v2_generation")
    GenerationInterval = module.GenerationInterval
    gi = GenerationInterval(mean_days=17.0, std_days=6.0, max_days=60)
    pmf = gi.get_pmf()
    assert len(pmf) == 60, f"Expected 60 days, got {len(pmf)}"
    assert abs(pmf.sum() - 1.0) < 1e-6, f"PMF should sum to 1, got {pmf.sum()}"
    assert pmf.min() >= 0, "PMF values should be non-negative"
    print("✓ GenerationInterval class tests passed")


def test_renewal_model():
    import numpy as np
    module = load_module(RE_ESTIMATION_SCRIPT, "step6_re_v2_renewal")
    GenerationInterval = module.GenerationInterval
    RenewalModel = module.RenewalModel
    gi = GenerationInterval(mean_days=17.0, std_days=6.0)
    model = RenewalModel(gi)
    incidence = np.array([10, 15, 20, 25, 30, 35, 40, 45, 50, 55])
    results = model.estimate_re(incidence)
    assert len(results) == 10, f"Expected 10 time points, got {len(results)}"
    assert 're_estimate' in results.columns
    assert 're_ci_lower' in results.columns
    assert 're_ci_upper' in results.columns
    print("✓ RenewalModel class tests passed")


@pytest.fixture
def output_dir(tmp_path: Path) -> Path:
    run_dir = tmp_path / "transmission_test_run"
    run_dir.mkdir()
    return run_dir


def test_synthetic_data_generation(output_dir: Path):
    generate_synthetic_test_data(output_dir)
    assert (output_dir / "synthetic_cases.tsv").exists(), "synthetic_cases.tsv not created"
    assert (output_dir / "synthetic_covariates.tsv").exists(), "synthetic_covariates.tsv not created"
    assert (output_dir / "synthetic_data_metadata.json").exists(), "synthetic_data_metadata.json not created"
    print("✓ Synthetic data generation test passed")


def test_re_estimation_pipeline(output_dir: Path):
    run_re_estimation(output_dir)
    assert (output_dir / 're_output' / 'bp_country_year_re_trajectories.tsv').exists()
    assert (output_dir / 're_output' / 'bp_re_summary_statistics.tsv').exists()
    print("✓ Rₑ estimation pipeline test passed")


def test_re_estimation_rejects_annual_inputs_without_override(output_dir: Path):
    generate_synthetic_test_data(output_dir)
    cmd = [
        sys.executable,
        str(RE_ESTIMATION_SCRIPT),
        "-i",
        str(output_dir / "synthetic_cases.tsv"),
        "-o",
        str(output_dir / "re_output_fail"),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    assert result.returncode != 0
    assert "Annual country-year case totals cannot support renewal-based R_e estimation" in result.stderr


def test_transmission_model_fitting(output_dir: Path):
    run_re_estimation(output_dir)
    re_input = output_dir / 're_output' / 'bp_country_year_re_trajectories.tsv'
    covariates_input = output_dir / 'synthetic_covariates.tsv'
    cmd = [
        sys.executable,
        str(TRANSMISSION_MODEL_SCRIPT),
        "--re-input",
        str(re_input),
        "--covariates",
        str(covariates_input),
        "-o",
        str(output_dir / "models"),
        "--allow-development-re-input",
    ]
    assert run_command(cmd)
    assert (output_dir / 'models' / 'bp_transmission_model_results.json').exists()
    print("✓ Transmission model fitting test passed")


def test_transmission_model_rejects_unsupported_re_metadata(output_dir: Path):
    run_re_estimation(output_dir)
    re_input = output_dir / 're_output' / 'bp_country_year_re_trajectories.tsv'
    covariates_input = output_dir / 'synthetic_covariates.tsv'
    cmd = [
        sys.executable,
        str(TRANSMISSION_MODEL_SCRIPT),
        "--re-input",
        str(re_input),
        "--covariates",
        str(covariates_input),
        "-o",
        str(output_dir / "models_fail"),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    assert result.returncode != 0
    assert "Transmission models are disabled for unsupported R_e inputs" in result.stderr
