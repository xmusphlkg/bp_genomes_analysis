#!/usr/bin/env python3
"""Cross-validation utility for transmission-model supporting summaries."""

import os, sys, json, argparse, logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Any
import random

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root

sys.path.insert(0, str(Path(__file__).parent.parent))
import pandas as pd
import numpy as np
from scipy import stats

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CrossValidator:
    def __init__(self, data: pd.DataFrame, n_iterations: int = 10, random_seed: int = 42):
        self.data = data.copy()
        self.n_iterations = n_iterations
        self.random_seed = random_seed
        np.random.seed(random_seed)
        random.seed(random_seed)
        logger.info(f"Initialized CrossValidator with {len(data)} observations")
    
    def k_fold_split(self, k: int = 5) -> List[Tuple[pd.DataFrame, pd.DataFrame]]:
        n = len(self.data)
        if n < 2:
            return []
        k = max(2, min(int(k), n))
        indices = np.random.permutation(n)
        fold_size = n // k
        folds = []
        for i in range(k):
            start_idx = i * fold_size
            end_idx = start_idx + fold_size if i < k - 1 else n
            test_indices = indices[start_idx:end_idx]
            train_indices = np.concatenate([indices[:start_idx], indices[end_idx:]])
            train_df = self.data.iloc[train_indices].copy()
            test_df = self.data.iloc[test_indices].copy()
            folds.append((train_df, test_df))
        logger.info(f"Created {k} folds with average size {n/k:.1f}")
        return folds
    
    def leave_one_country_out(self) -> List[Tuple[pd.DataFrame, pd.DataFrame]]:
        countries = self.data['country'].unique()
        folds = []
        for country in countries:
            train_df = self.data[self.data['country'] != country].copy()
            test_df = self.data[self.data['country'] == country].copy()
            if len(test_df) > 0 and len(train_df) > 0:
                folds.append((train_df, test_df))
        logger.info(f"Created {len(folds)} leave-one-country-out folds")
        return folds
    
    def temporal_holdout(self, holdout_years: int = 3) -> List[Tuple[pd.DataFrame, pd.DataFrame]]:
        years = sorted(self.data['year'].unique())
        folds = []
        for i in range(holdout_years, len(years)):
            test_year = years[i]
            train_years = years[max(0, i-holdout_years*3):i]
            train_df = self.data[self.data['year'].isin(train_years)].copy()
            test_df = self.data[self.data['year'] == test_year].copy()
            if len(test_df) > 0 and len(train_df) > 0:
                folds.append((train_df, test_df))
        logger.info(f"Created {len(folds)} temporal holdout folds")
        return folds
    
    def calculate_metrics(self, predictions: np.ndarray, actuals: np.ndarray) -> Dict[str, float]:
        if len(predictions) == 0 or len(actuals) == 0:
            return {}
        if len(predictions) != len(actuals):
            raise ValueError(
                f"prediction/actual length mismatch in cross-validation utility: "
                f"{len(predictions)} predictions vs {len(actuals)} actuals"
            )
        predictions, actuals = np.array(predictions), np.array(actuals)
        mae = np.mean(np.abs(predictions - actuals))
        rmse = np.sqrt(np.mean((predictions - actuals) ** 2))
        if len(predictions) > 1 and np.std(actuals) > 0:
            correlation = np.corrcoef(predictions, actuals)[0, 1]
            ss_res = np.sum((actuals - predictions) ** 2)
            ss_tot = np.sum((actuals - np.mean(actuals)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        else:
            correlation, r_squared = np.nan, np.nan
        non_zero_mask = actuals != 0
        mape = np.mean(np.abs((actuals[non_zero_mask] - predictions[non_zero_mask]) / actuals[non_zero_mask])) * 100 if np.any(non_zero_mask) else np.nan
        pred_std = np.std(predictions)
        coverage = np.mean((actuals >= predictions - 2*pred_std) & (actuals <= predictions + 2*pred_std))
        return {'mae': float(mae), 'rmse': float(rmse), 'correlation': float(correlation) if not np.isnan(correlation) else None,
                'r_squared': float(r_squared) if not np.isnan(r_squared) else None, 'mape': float(mape) if not np.isnan(mape) else None, 'coverage_95ci': float(coverage)}

class TransmissionModelValidator:
    def __init__(self, cv: CrossValidator):
        self.cv = cv

    @staticmethod
    def _observed_subset(test_data: pd.DataFrame, target_col: str) -> pd.DataFrame:
        return test_data.loc[test_data[target_col].notna()].copy()
    
    def naive_baseline(self, train_data: pd.DataFrame, test_data: pd.DataFrame, target_col: str = 're_mean') -> Tuple[np.ndarray, np.ndarray]:
        train_observed = train_data.loc[train_data[target_col].notna()].copy()
        test_observed = self._observed_subset(test_data, target_col)
        train_re = train_observed[target_col].to_numpy(dtype=float)
        if len(train_re) == 0 or test_observed.empty:
            return np.array([]), np.array([])
        mean_re = np.mean(train_re)
        predictions = np.full(len(test_observed), mean_re)
        actuals = test_observed[target_col].to_numpy(dtype=float)
        return predictions, actuals
    
    def country_specific_baseline(self, train_data: pd.DataFrame, test_data: pd.DataFrame, target_col: str = 're_mean') -> Tuple[np.ndarray, np.ndarray]:
        train_observed = train_data.loc[train_data[target_col].notna()].copy()
        test_observed = self._observed_subset(test_data, target_col)
        if train_observed.empty or test_observed.empty:
            return np.array([]), np.array([])
        predictions, actuals = [], []
        country_means = train_observed.groupby('country')[target_col].mean().to_dict()
        global_mean = float(train_observed[target_col].mean())
        for _, row in test_observed.iterrows():
            predictions.append(country_means.get(row['country'], global_mean))
            actuals.append(float(row[target_col]))
        return np.array(predictions, dtype=float), np.array(actuals, dtype=float)
    
    def temporal_trend_baseline(self, train_data: pd.DataFrame, test_data: pd.DataFrame, target_col: str = 're_mean') -> Tuple[np.ndarray, np.ndarray]:
        train_observed = train_data.loc[train_data[target_col].notna()].copy()
        test_observed = self._observed_subset(test_data, target_col)
        if train_observed.empty or test_observed.empty:
            return np.array([]), np.array([])
        predictions, actuals = [], []
        country_trends = {}
        for country in train_observed['country'].unique():
            country_data = train_observed[train_observed['country'] == country].copy()
            if len(country_data) >= 2:
                years = country_data['year'].to_numpy(dtype=float)
                values = country_data[target_col].to_numpy(dtype=float)
                try:
                    slope, intercept, _, _, _ = stats.linregress(years, values)
                except ValueError:
                    slope, intercept = 0.0, float(np.mean(values))
                country_trends[country] = (slope, intercept)
        if len(train_observed) >= 2:
            global_slope, global_intercept, _, _, _ = stats.linregress(
                train_observed['year'].to_numpy(dtype=float),
                train_observed[target_col].to_numpy(dtype=float),
            )
        else:
            global_slope, global_intercept = 0.0, float(train_observed[target_col].mean())
        for _, row in test_observed.iterrows():
            slope, intercept = country_trends.get(row['country'], (global_slope, global_intercept))
            predictions.append(slope * row['year'] + intercept)
            actuals.append(float(row[target_col]))
        return np.array(predictions, dtype=float), np.array(actuals, dtype=float)

def generate_synthetic_re_data(n_samples: int = 500) -> pd.DataFrame:
    np.random.seed(42)
    countries = ['USA', 'UK', 'Australia', 'France', 'Germany', 'Italy', 'Spain', 'Canada', 'Japan', 'Brazil', 'India', 'China', 'South Africa', 'Mexico', 'Argentina', 'Netherlands', 'Belgium', 'Sweden', 'Norway', 'Denmark']
    years = list(range(2008, 2025))
    data = []
    for _ in range(n_samples):
        country, year = np.random.choice(countries), np.random.choice(years)
        re_mean = max(0.5, min(3.0, 1.2 + np.random.normal(0, 0.2) - 0.01 * (year - 2008) + np.random.normal(0, 0.15)))
        data.append({'country': country, 'year': year, 're_mean': re_mean, 're_lower': re_mean * 0.8, 're_upper': re_mean * 1.2, 'incidence': np.random.exponential(15), 'coverage': np.random.uniform(70, 95)})
    return pd.DataFrame(data)

def safe_numeric_aggregate(values, func):
    """Safely aggregate numeric values, ignoring non-numeric."""
    numeric_vals = [v for v in values if isinstance(v, (int, float)) and not np.isnan(v) if isinstance(v, float)]
    return func(numeric_vals) if numeric_vals else None


def summarize_metric_rows(metric_rows: List[Dict[str, Any]], exclude_keys: set[str] | None = None) -> Dict[str, Dict[str, float]]:
    exclude_keys = exclude_keys or set()
    summary: Dict[str, Dict[str, float]] = {}
    if not metric_rows:
        return summary
    for key in metric_rows[0].keys():
        if key in exclude_keys:
            continue
        values = [m[key] for m in metric_rows if m.get(key) is not None]
        if not values:
            continue
        summary[key] = {
            'mean': float(np.mean(values)),
            'std': float(np.std(values)),
            'min': float(np.min(values)),
            'max': float(np.max(values)),
        }
    return summary


def load_validation_dataset(data_path: str, allow_synthetic_data: bool = False) -> tuple[pd.DataFrame, dict[str, Any]]:
    path = Path(data_path)
    if not path.exists():
        if not allow_synthetic_data:
            raise FileNotFoundError(
                f"Cross-validation input not found: {path}. "
                "Synthetic fallback is disabled unless --allow-synthetic-data is set."
            )
        logger.warning("Data file not found: %s, generating explicit synthetic data fallback", path)
        synthetic = generate_synthetic_re_data(n_samples=500)
        return synthetic, {
            "input_mode": "synthetic_fallback",
            "input_path": str(path),
            "input_format": "synthetic",
            "synthetic_data_used": True,
        }

    suffix = path.suffix.lower()
    sep = "\t" if suffix in {".tsv", ".tab", ".txt"} else ","
    df = pd.read_csv(path, sep=sep)
    rename_map = {}
    if "country" not in df.columns and "country_iso3" in df.columns:
        rename_map["country_iso3"] = "country"
    if "re_mean" not in df.columns and "re_estimate" in df.columns:
        rename_map["re_estimate"] = "re_mean"
    if "re_lower" not in df.columns and "re_ci_lower" in df.columns:
        rename_map["re_ci_lower"] = "re_lower"
    if "re_upper" not in df.columns and "re_ci_upper" in df.columns:
        rename_map["re_ci_upper"] = "re_upper"
    if "incidence" not in df.columns and "smoothed_incidence" in df.columns:
        rename_map["smoothed_incidence"] = "incidence"
    if "coverage" not in df.columns and "dtp3_coverage" in df.columns:
        rename_map["dtp3_coverage"] = "coverage"
    if rename_map:
        df = df.rename(columns=rename_map)

    required_columns = {"country", "year", "re_mean"}
    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(
            f"Cross-validation input is missing required columns after normalization: {sorted(missing)}"
        )

    df["country"] = df["country"].astype(str)
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df["re_mean"] = pd.to_numeric(df["re_mean"], errors="coerce")
    df = df.dropna(subset=["country", "year"]).copy()
    if df.empty:
        raise ValueError("Cross-validation input has no rows with non-missing country/year metadata")

    return df, {
        "input_mode": "observed",
        "input_path": str(path),
        "input_format": "tsv" if sep == "\t" else "csv",
        "synthetic_data_used": False,
    }


def run_comprehensive_validation(
    data_path: str,
    output_dir: str,
    n_iterations: int = 10,
    k_folds: int = 5,
    allow_synthetic_data: bool = False,
) -> Dict[str, Any]:
    logger.info("="*70 + "\nCOMPREHENSIVE CROSS-VALIDATION FRAMEWORK\n" + "="*70)
    logger.info(f"Loading data from {data_path}...")
    df, input_metadata = load_validation_dataset(data_path, allow_synthetic_data=allow_synthetic_data)
    logger.info("Loaded %d observations (%s)", len(df), input_metadata["input_mode"])
    
    cv = CrossValidator(df, n_iterations=n_iterations)
    validator = TransmissionModelValidator(cv)
    results = {'metadata': {'timestamp': datetime.now().isoformat(), 'n_observations': len(df), 'n_countries': df['country'].nunique() if 'country' in df.columns else None, 'n_years': df['year'].nunique() if 'year' in df.columns else None, 'n_iterations': n_iterations, 'k_folds': k_folds, **input_metadata}, 'validation_strategies': {}}
    
    logger.info("\n" + "="*70 + "\nSTRATEGY 1: K-FOLD CROSS-VALIDATION\n" + "="*70)
    kfold_results = {'strategy': 'k_fold', 'k': k_folds, 'baselines': {}}
    baseline_methods = [('naive_mean', validator.naive_baseline), ('country_specific', validator.country_specific_baseline), ('temporal_trend', validator.temporal_trend_baseline)]
    
    for baseline_name, baseline_func in baseline_methods:
        logger.info(f"\nTesting {baseline_name} baseline...")
        all_metrics = []
        for iteration in range(n_iterations):
            logger.info(f"  Iteration {iteration + 1}/{n_iterations}")
            folds = cv.k_fold_split(k=k_folds)
            iteration_metrics = []
            for fold_idx, (train_df, test_df) in enumerate(folds):
                predictions, actuals = baseline_func(train_df, test_df)
                if len(predictions) > 0 and len(actuals) > 0:
                    iteration_metrics.append(cv.calculate_metrics(predictions, actuals))
            if iteration_metrics:
                avg_metrics = {
                    key: float(np.mean(values))
                    for key, values in (
                        (k, [m[k] for m in iteration_metrics if m.get(k) is not None])
                        for k in iteration_metrics[0].keys()
                    )
                    if values
                }
                all_metrics.append(avg_metrics)
                if 'mae' in avg_metrics:
                    logger.info(f"    MAE: {avg_metrics['mae']:.4f}, R²: {avg_metrics.get('r_squared', 'N/A')}")
        if all_metrics:
            baseline_summary = summarize_metric_rows(all_metrics)
            for metric_name, metric_values in baseline_summary.items():
                metric_values.pop('min', None)
                metric_values.pop('max', None)
            kfold_results['baselines'][baseline_name] = baseline_summary
    results['validation_strategies']['k_fold'] = kfold_results
    
    logger.info("\n" + "="*70 + "\nSTRATEGY 2: LEAVE-ONE-COUNTRY-OUT VALIDATION\n" + "="*70)
    locoo_results = {'strategy': 'leave_one_country_out', 'baselines': {}, 'country_results': []}
    locoo_folds = cv.leave_one_country_out()
    logger.info(f"Testing with {len(locoo_folds)} country holdouts")
    
    for baseline_name, baseline_func in baseline_methods:
        logger.info(f"\nTesting {baseline_name} baseline...")
        all_metrics = []
        for train_df, test_df in locoo_folds:
            country = test_df['country'].iloc[0] if 'country' in test_df.columns else 'Unknown'
            predictions, actuals = baseline_func(train_df, test_df)
            if len(predictions) > 0 and len(actuals) > 0:
                metrics = cv.calculate_metrics(predictions, actuals)
                metrics['country'] = str(country)
                all_metrics.append(metrics)
                locoo_results['country_results'].append({'country': str(country), 'baseline': baseline_name, **metrics})
                logger.info(f"  {country}: MAE={metrics['mae']:.4f}")
        if all_metrics:
            locoo_results['baselines'][baseline_name] = summarize_metric_rows(all_metrics, exclude_keys={'country'})
    results['validation_strategies']['leave_one_country_out'] = locoo_results
    
    logger.info("\n" + "="*70 + "\nSTRATEGY 3: TEMPORAL HOLDOUT VALIDATION\n" + "="*70)
    temporal_results = {'strategy': 'temporal_holdout', 'holdout_years': 3, 'baselines': {}, 'year_results': []}
    temporal_folds = cv.temporal_holdout(holdout_years=3)
    logger.info(f"Testing with {len(temporal_folds)} temporal holdouts")
    
    for baseline_name, baseline_func in baseline_methods:
        logger.info(f"\nTesting {baseline_name} baseline...")
        all_metrics = []
        for train_df, test_df in temporal_folds:
            year = test_df['year'].iloc[0] if 'year' in test_df.columns else 'Unknown'
            predictions, actuals = baseline_func(train_df, test_df)
            if len(predictions) > 0 and len(actuals) > 0:
                metrics = cv.calculate_metrics(predictions, actuals)
                metrics['year'] = int(year)
                all_metrics.append(metrics)
                temporal_results['year_results'].append({'year': int(year), 'baseline': baseline_name, **metrics})
                logger.info(f"  Year {year}: MAE={metrics['mae']:.4f}")
        if all_metrics:
            baseline_summary = summarize_metric_rows(all_metrics, exclude_keys={'year'})
            for metric_name, metric_values in baseline_summary.items():
                metric_values.pop('min', None)
                metric_values.pop('max', None)
            temporal_results['baselines'][baseline_name] = baseline_summary
    results['validation_strategies']['temporal_holdout'] = temporal_results
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    report_file = output_path / f"cross_validation_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(report_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, default=str)
    
    logger.info("\n" + "="*70 + "\nCROSS-VALIDATION COMPLETE\n" + "="*70)
    logger.info(f"✓ Report saved to: {report_file}")
    
    summary = generate_summary_interpretation(results)
    logger.info("\n" + "="*70 + "\nINTERPRETATION SUMMARY\n" + "="*70)
    for key, value in summary.items():
        logger.info(f"{key}: {value}")
    results['interpretation'] = summary
    return results

def generate_summary_interpretation(results: Dict[str, Any]) -> Dict[str, str]:
    interpretation = {}
    kfold_results = results['validation_strategies'].get('k_fold', {})
    if kfold_results.get('baselines'):
        best_baseline = min(kfold_results['baselines'].items(), key=lambda x: x[1]['mae']['mean'])
        interpretation['best_performing_method'] = best_baseline[0]
        interpretation['overall_mae'] = f"{best_baseline[1]['mae']['mean']:.4f} ± {best_baseline[1]['mae']['std']:.4f}"
        r2_val = best_baseline[1].get('r_squared', {}).get('mean')
        interpretation['model_quality'] = "GOOD" if r2_val and r2_val > 0.7 else "MODERATE" if r2_val and r2_val > 0.5 else "NEEDS IMPROVEMENT"
        mae_std = best_baseline[1]['mae']['std']
        interpretation['overfitting_risk'] = "LOW" if mae_std < 0.1 else "MODERATE" if mae_std < 0.2 else "HIGH"
    locoo_results = results['validation_strategies'].get('leave_one_country_out', {})
    if locoo_results.get('baselines') and 'best_baseline' in locals():
        locoo_mae = locoo_results['baselines'][best_baseline[0]]['mae']['mean']
        interpretation['generalizability'] = "GOOD" if locoo_mae < best_baseline[1]['mae']['mean'] * 1.2 else "LIMITED"
    if results['validation_strategies'].get('temporal_holdout'):
        interpretation['temporal_stability'] = "ASSESSED"
    return interpretation

def main():
    parser = argparse.ArgumentParser(description='Cross-validation framework for transmission models')
    step6_root = project_module_data_root("step6_epi_transmission")
    parser.add_argument('--data', type=str, default=str(step6_root / 'outputs' / 'bp_country_year_re_trajectories.tsv'), help='Input data file')
    parser.add_argument('--k', type=int, default=5, help='Number of folds (default: 5)')
    parser.add_argument('--iterations', type=int, default=10, help='Number of CV iterations (default: 10)')
    parser.add_argument('--output', type=str, default=str(step6_root / 'outputs' / 'cross_validation'), help='Output directory')
    parser.add_argument('--seed', type=int, default=42, help='Random seed (default: 42)')
    parser.add_argument(
        '--allow-synthetic-data',
        action='store_true',
        help='Allow explicit synthetic fallback if the requested input file is missing.',
    )
    args = parser.parse_args()
    
    print("="*70 + "\nCROSS-VALIDATION FRAMEWORK FOR TRANSMISSION MODELS\n" + "="*70)
    print(f"\nConfiguration:\n  Input: {args.data}\n  K-folds: {args.k}\n  Iterations: {args.iterations}\n  Output: {args.output}\n  Seed: {args.seed}\n")
    
    try:
        results = run_comprehensive_validation(
            data_path=args.data,
            output_dir=args.output,
            n_iterations=args.iterations,
            k_folds=args.k,
            allow_synthetic_data=args.allow_synthetic_data,
        )
        print("\n" + "="*70 + "\n✓ CROSS-VALIDATION COMPLETED SUCCESSFULLY\n" + "="*70)
    except Exception as e:
        logger.error(f"Cross-validation failed: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == '__main__':
    main()
