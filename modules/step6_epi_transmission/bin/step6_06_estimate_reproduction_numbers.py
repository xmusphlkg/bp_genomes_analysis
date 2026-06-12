#!/usr/bin/env python3
"""
Estimate effective reproduction number (Rₑ) trajectories from country-year pertussis case data.

Implements renewal equation model following Cori et al. (2013) and Thompson et al. (2019).
Rₑ(t) = I(t) / Σ_{s=1}^{∞} I(t-s) × w(s)
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import gamma

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class GenerationInterval:
    """Gamma-distributed generation interval for pertussis."""
    
    def __init__(self, mean_days: float = 17.0, std_days: float = 6.0, max_days: int = 60):
        self.mean_days = mean_days
        self.std_days = std_days
        self.max_days = max_days
        alpha = (mean_days / std_days) ** 2
        beta = std_days ** 2 / mean_days
        self.alpha = alpha
        self.beta = beta
        self._compute_pmf()
    
    def _compute_pmf(self):
        days = np.arange(0, self.max_days + 1)
        cdf_vals = gamma.cdf(days, a=self.alpha, scale=self.beta)
        pmf = np.diff(cdf_vals)
        self.pmf = pmf / pmf.sum()
        self.days = days[:-1]
    
    def get_pmf(self) -> np.ndarray:
        return self.pmf


class RenewalModel:
    """Renewal equation model for Rₑ estimation."""
    
    def __init__(self, generation_interval: GenerationInterval, prior_shape: float = 1.0, 
                 prior_rate: float = 1.0, window_size: Optional[int] = None):
        self.generation_interval = generation_interval
        self.prior_shape = prior_shape
        self.prior_rate = prior_rate
        self.window_size = window_size or len(generation_interval.get_pmf())
    
    def compute_infectivity(self, incidence: np.ndarray, t: int) -> float:
        if t == 0:
            return 0.0
        max_lag = min(t, self.window_size, len(self.generation_interval.pmf))
        infectivity = sum(incidence[t - s] * self.generation_interval.pmf[s - 1] for s in range(1, max_lag + 1) if t - s >= 0)
        return infectivity
    
    def estimate_re(self, incidence: np.ndarray) -> pd.DataFrame:
        n_time = len(incidence)
        re_estimates, re_lower, re_upper = [], [], []
        
        for t in range(n_time):
            if t < 1:
                re_estimates.extend([np.nan, np.nan, np.nan])
                continue
            
            infectivity = self.compute_infectivity(incidence, t)
            if infectivity > 0:
                posterior_shape = incidence[t] + self.prior_shape
                posterior_rate = infectivity + self.prior_rate
                re_mean = posterior_shape / posterior_rate
                ci_lower = gamma.ppf(0.025, posterior_shape, scale=1/posterior_rate)
                ci_upper = gamma.ppf(0.975, posterior_shape, scale=1/posterior_rate)
                re_estimates.extend([re_mean, ci_lower, ci_upper])
            else:
                re_estimates.extend([np.nan, np.nan, np.nan])
        
        return pd.DataFrame({
            'time_index': range(n_time), 'incidence': incidence,
            're_estimate': re_estimates[::3], 're_ci_lower': re_estimates[1::3], 're_ci_upper': re_estimates[2::3]
        })


def load_public_health_data(input_path: Path) -> pd.DataFrame:
    df = pd.read_csv(input_path, sep='\t')
    column_mapping = {
        'iso3': 'country',
        'iso3_code': 'country',
        'country_iso3': 'country',
        'country_name': 'country_name',
        'case_count': 'cases',
        'reported_cases': 'cases',
        'pop': 'population',
    }
    df = df.rename(columns=column_mapping)
    
    for col in ['country', 'year', 'cases']:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
    
    df['year'] = pd.to_numeric(df['year'], errors='coerce')
    df['cases'] = pd.to_numeric(df['cases'], errors='coerce')
    df = df.dropna(subset=['country', 'year', 'cases'])
    df.loc[df['cases'] < 0, 'cases'] = 0
    
    logger.info(f"Loaded data for {df['country'].nunique()} countries, {df['year'].nunique()} years")
    return df


def reconstruct_weekly_incidence(annual_cases: pd.DataFrame, method: str = 'uniform') -> pd.DataFrame:
    weekly_data = []
    for (country, year), group in annual_cases.groupby(['country', 'year']):
        annual_total = group['cases'].values[0]
        for week in range(1, 53):
            if method == 'uniform':
                weekly_cases = annual_total / 52.0
            elif method == 'seasonal':
                seasonal_factor = 1 + 0.3 * np.sin(2 * np.pi * (week - 10) / 52)
                weekly_cases = (annual_total / 52.0) * seasonal_factor
            else:
                raise ValueError(f"Unknown method: {method}")
            weekly_data.append({'country': country, 'year': year, 'week': week, 'cases': weekly_cases})
    return pd.DataFrame(weekly_data)


def estimate_country_re_trajectories(weekly_incidence: pd.DataFrame, gen_interval: GenerationInterval, output_path: Path) -> pd.DataFrame:
    model = RenewalModel(gen_interval)
    all_results = []
    
    for country in weekly_incidence['country'].unique():
        country_data = weekly_incidence[weekly_incidence['country'] == country].sort_values(['year', 'week']).copy()
        country_data['time_index'] = ((country_data['year'] - country_data['year'].min()) * 52 + country_data['week'])
        
        max_time = int(country_data['time_index'].max())
        incidence_array = np.zeros(max_time + 1)
        for _, row in country_data.iterrows():
            incidence_array[int(row['time_index'])] = row['cases']
        
        re_results = model.estimate_re(incidence_array)
        re_results['country'] = country
        re_results['year'] = re_results['time_index'] // 52 + country_data['year'].min()
        re_results['week'] = (re_results['time_index'] % 52) + 1
        all_results.append(re_results)
    
    results_df = pd.concat(all_results, ignore_index=True)
    results_df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Saved Rₑ estimates to {output_path}")
    return results_df


def generate_summary_statistics(re_estimates: pd.DataFrame, output_path: Path) -> pd.DataFrame:
    valid = re_estimates[re_estimates['re_estimate'].notna()].copy()
    summary = valid.groupby(['country', 'year']).agg(
        re_mean=('re_estimate', 'mean'), re_median=('re_estimate', 'median'),
        re_std=('re_estimate', 'std'), re_min=('re_estimate', 'min'), re_max=('re_estimate', 'max'),
        n_weeks=('re_estimate', 'count'), total_cases=('incidence', 'sum')
    ).reset_index()
    summary['n_weeks_above_1.0'] = valid.groupby(['country', 'year']).apply(lambda x: (x['re_estimate'] > 1.0).sum()).values
    summary.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Saved summary statistics to {output_path}")
    return summary


def main():
    parser = argparse.ArgumentParser(description='Estimate Rₑ from pertussis case data')
    parser.add_argument('--input', '-i', type=Path, required=True)
    parser.add_argument('--output-dir', '-o', type=Path, required=True)
    parser.add_argument('--gi-mean', type=float, default=17.0)
    parser.add_argument('--gi-sd', type=float, default=6.0)
    parser.add_argument('--disaggregation-method', type=str, choices=['uniform', 'seasonal'], default='uniform')
    parser.add_argument('--verbose', '-v', action='store_true')
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    public_health_data = load_public_health_data(args.input)
    gen_interval = GenerationInterval(mean_days=args.gi_mean, std_days=args.gi_sd)
    weekly_incidence = reconstruct_weekly_incidence(public_health_data, method=args.disaggregation_method)
    
    re_output = args.output_dir / 'bp_country_year_re_trajectories.tsv'
    re_estimates = estimate_country_re_trajectories(weekly_incidence, gen_interval, re_output)
    
    summary_output = args.output_dir / 'bp_re_summary_statistics.tsv'
    generate_summary_statistics(re_estimates, summary_output)
    
    metadata = {
        'run_timestamp': datetime.now().isoformat(), 'input_file': str(args.input),
        'gi_mean': args.gi_mean, 'gi_sd': args.gi_sd, 'disaggregation_method': args.disaggregation_method,
        'n_countries': int(public_health_data['country'].nunique()),
        'n_years': int(public_health_data['year'].nunique()),
        'year_range': f"{public_health_data['year'].min()}-{public_health_data['year'].max()}"
    }
    with open(args.output_dir / 'bp_re_run_metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    
    logger.info("Rₑ estimation completed successfully")
    return 0


if __name__ == '__main__':
    sys.exit(main())
