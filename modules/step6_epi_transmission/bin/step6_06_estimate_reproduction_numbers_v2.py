#!/usr/bin/env python3
"""
Estimate effective reproduction number (Rₑ) trajectories from country-year pertussis case data.

IMPROVED VERSION with robustness features:
- Incidence smoothing to prevent artificial spikes
- Burn-in period exclusion for early unreliable estimates  
- Biologically plausible range truncation (0 < Re < 20)
- Quality flag tracking for transparency
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import gamma
from scipy.optimize import minimize_scalar

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class GenerationInterval:
    """Represents the generation interval distribution for pertussis."""
    
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
    """Renewal equation model for estimating time-varying reproduction numbers."""
    
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
        infectivity = 0.0
        for s in range(1, max_lag + 1):
            if t - s >= 0:
                infectivity += incidence[t - s] * self.generation_interval.pmf[s - 1]
        
        return infectivity
    
    @staticmethod
    def smooth_incidence(incidence: np.ndarray, window: int = 7) -> np.ndarray:
        """Apply centered moving average to smooth incidence data."""
        if len(incidence) < window:
            return incidence.copy()
        
        smoothed = np.convolve(incidence, np.ones(window)/window, mode='same')
        
        # Handle edges with smaller windows
        half_window = window // 2
        for i in range(half_window):
            left_window = i + half_window + 1
            smoothed[i] = np.mean(incidence[:left_window])
            smoothed[-(i+1)] = np.mean(incidence[-(left_window):])
        
        return smoothed
    
    def estimate_re(self, incidence: np.ndarray, burn_in_weeks: int = 4,
                   max_plausible_re: float = 20.0) -> pd.DataFrame:
        """
        Estimate time-varying reproduction numbers with robustness improvements.
        
        Args:
            incidence: Array of incident cases (weekly)
            burn_in_weeks: Number of initial weeks to exclude (unreliable due to limited history)
            max_plausible_re: Maximum biologically plausible Re for pertussis
            
        Returns:
            DataFrame with Rₑ estimates, credible intervals, and quality flags
        """
        n_time = len(incidence)
        
        # Apply smoothing to reduce artificial spikes from uniform disaggregation
        smoothed_incidence = self.smooth_incidence(incidence, window=7)
        
        re_estimates = []
        re_lower = []
        re_upper = []
        quality_flags = []
        
        for t in range(n_time):
            flag = "OK"
            
            # Mark burn-in period as unreliable
            if t < burn_in_weeks:
                re_estimates.append(np.nan)
                re_lower.append(np.nan)
                re_upper.append(np.nan)
                quality_flags.append("BURN_IN_PERIOD")
                continue
            
            # Use smoothed incidence for estimation
            current_incidence = smoothed_incidence[t]
            infectivity = self.compute_infectivity(smoothed_incidence, t)
            
            if infectivity <= 0:
                re_estimates.append(np.nan)
                re_lower.append(np.nan)
                re_upper.append(np.nan)
                quality_flags.append("ZERO_INFECTIVITY")
                continue
            
            # Bayesian posterior estimation
            posterior_shape = current_incidence + self.prior_shape
            posterior_rate = infectivity + self.prior_rate
            re_mean = posterior_shape / posterior_rate
            
            # Check biological plausibility
            if re_mean > max_plausible_re:
                logger.debug(f"Week {t}: Re={re_mean:.2f} capped at {max_plausible_re} (implausible)")
                re_mean = max_plausible_re
                flag = "CAPPED_AT_MAX"
            elif re_mean < 0:
                re_mean = 0.0
                flag = "SET_TO_ZERO"
            
            # Calculate credible intervals
            ci_lower = gamma.ppf(0.025, posterior_shape, scale=1/posterior_rate)
            ci_upper = gamma.ppf(0.975, posterior_shape, scale=1/posterior_rate)
            
            re_estimates.append(re_mean)
            re_lower.append(ci_lower)
            re_upper.append(ci_upper)
            quality_flags.append(flag)
        
        return pd.DataFrame({
            'time_index': range(n_time),
            'raw_incidence': incidence,
            'smoothed_incidence': smoothed_incidence,
            're_estimate': re_estimates,
            're_ci_lower': re_lower,
            're_ci_upper': re_upper,
            'quality_flag': quality_flags
        })


def load_public_health_data(input_path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    df = pd.read_csv(input_path, sep='\t')
    column_mapping = {
        'iso3': 'country', 'iso3_code': 'country', 'country_iso3': 'country',
        'country_name': 'country_name', 'case_count': 'cases',
        'reported_cases': 'cases', 'pop': 'population',
        'epi_week': 'week', 'week_num': 'week', 'week_number': 'week',
    }
    df = df.rename(columns=column_mapping)
    
    for col in ['country', 'year', 'cases']:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
    
    df['year'] = pd.to_numeric(df['year'], errors='coerce')
    df['cases'] = pd.to_numeric(df['cases'], errors='coerce')
    df = df.dropna(subset=['country', 'year', 'cases'])
    df.loc[df['cases'] < 0, 'cases'] = 0

    metadata: dict[str, object] = {
        'input_temporal_resolution': 'annual_country_year',
        'observed_subannual_incidence': False,
    }

    if 'week' in df.columns:
        df['week'] = pd.to_numeric(df['week'], errors='coerce')
        weekly = df.dropna(subset=['week']).copy()
        weekly = weekly[weekly['week'].between(1, 53)].copy()
        if not weekly.empty and weekly.groupby(['country', 'year'], dropna=False).size().gt(1).any():
            weekly = weekly.sort_values(['country', 'year', 'week']).reset_index(drop=True)
            metadata = {
                'input_temporal_resolution': 'observed_weekly',
                'observed_subannual_incidence': True,
            }
            logger.info(
                "Loaded observed weekly incidence for %d countries across %d years",
                weekly['country'].nunique(),
                weekly['year'].nunique(),
            )
            return weekly, metadata

    logger.info(f"Loaded annual case totals for {df['country'].nunique()} countries, {df['year'].nunique()} years")
    annual = df.sort_values(['country', 'year']).reset_index(drop=True)
    return annual, metadata


def reconstruct_weekly_incidence(
    annual_cases: pd.DataFrame,
    method: str = 'uniform',
    *,
    allow_annual_disaggregation: bool = False,
) -> pd.DataFrame:
    if not allow_annual_disaggregation:
        raise ValueError(
            "Annual country-year case totals cannot support renewal-based R_e estimation. "
            "Use observed weekly incidence, or rerun with --allow-annual-disaggregation for "
            "development-only smoke tests."
        )
    weekly_data = []
    for (country, year), group in annual_cases.groupby(['country', 'year']):
        total_cases = group['cases'].values[0]
        population = group['population'].values[0] if 'population' in group.columns else None
        
        if method == 'uniform':
            weekly_cases = [total_cases / 52.0] * 52
        else:
            weekly_cases = [total_cases / 52.0] * 52
            logger.warning(f"Seasonal disaggregation not yet implemented, using uniform for {country} {year}")
        
        for week in range(52):
            row = {'country': country, 'year': year, 'week': week + 1, 'cases': weekly_cases[week]}
            if population is not None:
                row['population'] = population
            weekly_data.append(row)
    
    result = pd.DataFrame(weekly_data)
    logger.info(f"Generated weekly data: {len(result)} rows ({len(result)//52} country-years)")
    return result


def estimate_country_re_trajectories(weekly_incidence: pd.DataFrame, 
                                     gen_interval: GenerationInterval,
                                     output_path: Path,
                                     burn_in_weeks: int = 4,
                                     max_plausible_re: float = 20.0) -> pd.DataFrame:
    model = RenewalModel(gen_interval)
    all_trajectories = []
    
    grouped = weekly_incidence.groupby('country', dropna=False)
    logger.info(f"Estimating continuous Rₑ trajectories for {len(grouped)} countries...")
    
    for idx, (country, group) in enumerate(grouped, start=1):
        order_columns = ['year']
        if 'week' in group.columns:
            order_columns.append('week')
        group = group.sort_values(order_columns).reset_index(drop=True)
        incidence = group['cases'].values
        
        if len(incidence) <= burn_in_weeks:
            logger.warning(f"Skipping {country}: only %d weekly observations after ordering", len(incidence))
            continue
        
        trajectory = model.estimate_re(
            incidence,
            burn_in_weeks=burn_in_weeks,
            max_plausible_re=max_plausible_re,
        )
        trajectory['country'] = np.repeat(country, len(trajectory))
        trajectory['year'] = group['year'].to_numpy()
        if 'week' in group.columns:
            trajectory['week'] = group['week'].to_numpy()
        if 'population' in group.columns:
            trajectory['population'] = group['population'].to_numpy()
        trajectory['country_time_index'] = np.arange(len(group), dtype=int)
        all_trajectories.append(trajectory)
        
        if idx % 100 == 0:
            logger.info(f"Processed {idx}/{len(grouped)} countries")
    
    combined = pd.concat(all_trajectories, ignore_index=True)
    combined.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Saved Rₑ trajectories to {output_path}")
    
    abnormal = combined[(combined['re_estimate'] > 20) | (combined['re_estimate'] < 0)]
    logger.info(f"Abnormal RE values after filtering: {len(abnormal)} ({len(abnormal)/len(combined)*100:.2f}%)")
    
    return combined


def generate_summary_statistics(re_trajectories: pd.DataFrame, output_path: Path):
    summary_stats = []
    for (country, year), group in re_trajectories.groupby(['country', 'year']):
        valid_re = group['re_estimate'].dropna()
        
        stats = {
            'country': country, 'year': year,
            'mean_re': valid_re.mean() if len(valid_re) > 0 else np.nan,
            'median_re': valid_re.median() if len(valid_re) > 0 else np.nan,
            'min_re': valid_re.min() if len(valid_re) > 0 else np.nan,
            'max_re': valid_re.max() if len(valid_re) > 0 else np.nan,
            'std_re': valid_re.std() if len(valid_re) > 0 else np.nan,
            'num_weeks_valid': len(valid_re),
            'num_weeks_capped': len(group[group['quality_flag'] == 'CAPPED_AT_MAX']),
            'percent_capped': len(group[group['quality_flag'] == 'CAPPED_AT_MAX']) / len(group) * 100 if len(group) > 0 else np.nan,
        }
        summary_stats.append(stats)
    
    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Saved summary statistics to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Estimate time-varying reproduction numbers')
    parser.add_argument('--input', '-i', type=Path, required=True, help='Input TSV with country-year case counts')
    parser.add_argument('--output-dir', '-o', type=Path, required=True)
    parser.add_argument('--gi-mean', type=float, default=17.0, help='Generation interval mean (days)')
    parser.add_argument('--gi-sd', type=float, default=6.0, help='Generation interval SD (days)')
    parser.add_argument('--disaggregation-method', type=str, choices=['uniform', 'seasonal'], default='uniform')
    parser.add_argument('--burn-in-weeks', type=int, default=4, help='Weeks to exclude at start of each year')
    parser.add_argument('--max-re', type=float, default=20.0, help='Maximum biologically plausible Re')
    parser.add_argument(
        '--allow-annual-disaggregation',
        action='store_true',
        help='Development-only override that permits synthetic weekly disaggregation from annual totals',
    )
    parser.add_argument('--verbose', '-v', action='store_true')
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    public_health_data, input_metadata = load_public_health_data(args.input)
    gen_interval = GenerationInterval(mean_days=args.gi_mean, std_days=args.gi_sd)
    synthetic_disaggregation_used = not bool(input_metadata['observed_subannual_incidence'])
    if synthetic_disaggregation_used:
        weekly_incidence = reconstruct_weekly_incidence(
            public_health_data,
            method=args.disaggregation_method,
            allow_annual_disaggregation=args.allow_annual_disaggregation,
        )
    else:
        weekly_incidence = public_health_data.copy()
    
    re_output = args.output_dir / 'bp_country_year_re_trajectories.tsv'
    re_estimates = estimate_country_re_trajectories(
        weekly_incidence,
        gen_interval,
        re_output,
        burn_in_weeks=args.burn_in_weeks,
        max_plausible_re=args.max_re,
    )
    
    summary_output = args.output_dir / 'bp_re_summary_statistics.tsv'
    generate_summary_statistics(re_estimates, summary_output)

    metadata = {
        'run_timestamp': datetime.now().isoformat(),
        'input': str(args.input),
        'output_dir': str(args.output_dir),
        'n_country_years': int(re_estimates[['country', 'year']].drop_duplicates().shape[0]),
        'n_countries': int(re_estimates['country'].nunique()),
        'generation_interval_mean_days': float(args.gi_mean),
        'generation_interval_sd_days': float(args.gi_sd),
        'disaggregation_method': args.disaggregation_method,
        'burn_in_weeks': int(args.burn_in_weeks),
        'max_re': float(args.max_re),
        'input_temporal_resolution': input_metadata['input_temporal_resolution'],
        'observed_subannual_incidence': bool(input_metadata['observed_subannual_incidence']),
        'synthetic_disaggregation_used': bool(synthetic_disaggregation_used),
        'manuscript_supported': bool(input_metadata['observed_subannual_incidence']),
        'development_override_allow_annual_disaggregation': bool(args.allow_annual_disaggregation),
    }
    with open(args.output_dir / 'bp_re_run_metadata.json', 'w') as handle:
        json.dump(metadata, handle, indent=2)
    
    logger.info("Rₑ estimation complete!")


if __name__ == '__main__':
    main()
