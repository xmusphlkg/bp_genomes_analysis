#!/usr/bin/env python3
"""Generate synthetic pertussis case data for testing transmission dynamics module."""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def generate_synthetic_cases(n_countries: int = 20, n_years: int = 20, seed: int = 42) -> pd.DataFrame:
    np.random.seed(seed)
    countries = [f"Country_{chr(65+i)}" for i in range(n_countries)]
    years = list(range(2000, 2000 + n_years))
    
    data = []
    for country in countries:
        baseline_incidence = np.random.uniform(5, 50)
        trend = np.random.uniform(-0.02, 0.02)
        for year in years:
            seasonal_effect = 1 + 0.3 * np.sin(2 * np.pi * (year - 2000) / 5 + np.random.uniform(0, 2*np.pi))
            noise = np.random.normal(1, 0.2)
            cases = max(0, baseline_incidence * (1 + trend * (year - 2000)) * seasonal_effect * noise)
            data.append({'country': country, 'year': year, 'cases': int(cases)})
    
    df = pd.DataFrame(data)
    print(f"Generated {len(df)} observations for {n_countries} countries over {n_years} years")
    return df


def generate_synthetic_covariates(countries: list, years: list, seed: int = 42) -> pd.DataFrame:
    np.random.seed(seed + 1)
    data = []
    for country in countries:
        base_dtp3 = np.random.uniform(60, 95)
        trend = np.random.uniform(0, 1.5)
        for year in years:
            dtp3 = np.clip(base_dtp3 + trend * (year - years[0]) + np.random.normal(0, 3), 30, 99)
            reported_cases = max(10, int(np.exp(np.random.normal(8.5, 0.8))))
            ipw_prevalence = np.clip(np.random.beta(1.5, 2.5), 0, 1)
            naive_prevalence = np.clip(ipw_prevalence + np.random.normal(0, 0.05), 0, 1)
            genomes_per_case = np.clip(np.random.lognormal(mean=-5.3, sigma=0.6), 1e-5, None)
            data.append({
                'country': country,
                'year': year,
                'dtp3_coverage': round(dtp3, 1),
                'reported_cases': reported_cases,
                'post_covid_period': int(year >= 2024),
                'n_genomes_prn_interpretable': int(np.random.randint(6, 35)),
                'ipw_prevalence': round(ipw_prevalence, 6),
                'naive_prevalence': round(naive_prevalence, 6),
                'genomes_per_case': round(float(genomes_per_case), 6),
            })
    df = pd.DataFrame(data)
    df['ap_exposure_v1_score'] = (
        (df['dtp3_coverage'] - df['dtp3_coverage'].mean()) / max(df['dtp3_coverage'].std(), 1e-6) +
        0.7 * ((df['year'] - df['year'].mean()) / max(df['year'].std(), 1e-6)) +
        0.5 * df['post_covid_period']
    )
    return df


def main():
    parser = argparse.ArgumentParser(description='Generate synthetic test data')
    parser.add_argument('--n-countries', type=int, default=20)
    parser.add_argument('--n-years', type=int, default=20)
    parser.add_argument('--output-dir', '-o', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    cases_df = generate_synthetic_cases(args.n_countries, args.n_years, seed=args.seed)
    cases_df.to_csv(args.output_dir / 'synthetic_cases.tsv', sep='\t', index=False)
    
    countries = cases_df['country'].unique()
    years = sorted(cases_df['year'].unique())
    covariates_df = generate_synthetic_covariates(list(countries), list(years), seed=args.seed)
    covariates_df.to_csv(args.output_dir / 'synthetic_covariates.tsv', sep='\t', index=False)
    
    metadata = {
        'generated_at': datetime.now().isoformat(),
        'n_countries': args.n_countries,
        'n_years': args.n_years,
        'seed': args.seed,
        'files': ['synthetic_cases.tsv', 'synthetic_covariates.tsv']
    }
    with open(args.output_dir / 'synthetic_data_metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"✓ Generated synthetic data in {args.output_dir}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
