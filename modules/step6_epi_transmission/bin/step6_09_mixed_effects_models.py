#!/usr/bin/env python3
"""
step6_09_mixed_effects_models.py — Fit mixed-effects and clustered GLMs.

Addresses reviewer concerns:
- GLM does not account for within-country correlation
- Small denominators destabilize estimates
- Need vaccine program covariates beyond DTP3

Models fitted:
1. Binomial GLM (baseline replication)
2. GLM with country-cluster sandwich SEs
3. GEE with exchangeable correlation (clustered by country)
4. GLM excluding cells with n < 10
5. Country fixed-effects model
6. Extended model with acellular_vs_whole_cell covariate

Per-country logistic regressions for countries with ≥3 observations.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "outputs"
INPUT_PATH = OUTPUT_DIR / "bp_country_year_analysis_input.tsv"


def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df[df["n_genomes_prn_interpretable"] > 0].copy()
    for col in ["n_prn_disrupted", "n_genomes_prn_interpretable",
                "dtp3_coverage", "reported_cases", "incidence_per_100k"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["n_prn_disrupted", "n_genomes_prn_interpretable", "dtp3_coverage"])
    return df


def _prepare_X(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    X = df[[c for c in cols if c in df.columns]].copy()
    if "reported_cases" in X.columns:
        X["reported_cases"] = np.log1p(X["reported_cases"].fillna(0))
    if "incidence_per_100k" in X.columns:
        X["incidence_per_100k"] = np.log1p(X["incidence_per_100k"].fillna(0))
    return sm.add_constant(X)


def _binomial_y(df: pd.DataFrame):
    return np.column_stack([
        df["n_prn_disrupted"].values,
        (df["n_genomes_prn_interpretable"] - df["n_prn_disrupted"]).values,
    ])


def fit_glm(df, preds, label, cov_type: str | None = None, cov_kwds: dict[str, Any] | None = None):
    X = _prepare_X(df, preds)
    y = _binomial_y(df)
    m = sm.GLM(y, X.values, family=sm.families.Binomial())
    fit_kwargs: dict[str, Any] = {}
    if cov_type is not None:
        fit_kwargs["cov_type"] = cov_type
    if cov_kwds:
        fit_kwargs["cov_kwds"] = cov_kwds
    r = m.fit(**fit_kwargs)
    return {
        "model": label,
        "converged": True,
        "n_obs": int(r.nobs),
        "covariance_type": getattr(r, "cov_type", "nonrobust"),
        "params": dict(zip(X.columns, r.params.tolist())),
        "bse": dict(zip(X.columns, r.bse.tolist())),
        "pvalues": dict(zip(X.columns, r.pvalues.tolist())),
        "conf_low": dict(zip(X.columns, r.conf_int()[:, 0].tolist())),
        "conf_high": dict(zip(X.columns, r.conf_int()[:, 1].tolist())),
        "aic": float(r.aic),
        "bic": float(r.bic),
    }


def fit_cluster_robust(df, preds):
    X = _prepare_X(df, preds)
    y = _binomial_y(df)
    m = sm.GLM(y, X.values, family=sm.families.Binomial())
    n_clusters = int(df["country_iso3"].nunique())
    if n_clusters < 3:
        return {"model": "glm_country_cluster_sandwich", "converged": False, "error": "fewer than 3 country clusters"}
    try:
        r = m.fit(
            cov_type="cluster",
            cov_kwds={"groups": df["country_iso3"].astype(str).to_numpy(dtype=object)},
        )
        z = r.params / r.bse
        return {
            "model": "glm_country_cluster_sandwich",
            "converged": True,
            "n_obs": int(r.nobs),
            "n_clusters": n_clusters,
            "covariance_type": getattr(r, "cov_type", "cluster"),
            "params": dict(zip(X.columns, r.params.tolist())),
            "bse_robust": dict(zip(X.columns, r.bse.tolist())),
            "pvalues_robust": dict(zip(X.columns, r.pvalues.tolist())),
            "z_robust": dict(zip(X.columns, z.tolist())),
        }
    except Exception as e:
        return {"model": "glm_country_cluster_sandwich", "converged": False, "error": str(e)}


def fit_gee(df, preds):
    from statsmodels.genmod.generalized_estimating_equations import GEE
    from statsmodels.genmod.cov_struct import Exchangeable
    X = _prepare_X(df, preds)
    y = _binomial_y(df)
    groups = df["country_iso3"].values
    try:
        m = GEE(y, X.values, groups, family=sm.families.Binomial(),
                cov_struct=Exchangeable())
        r = m.fit()
        return {
            "model": "gee_exchangeable",
            "converged": True,
            "n_obs": int(r.nobs),
            "n_groups": int(r.n_groups),
            "params": dict(zip(X.columns, r.params.tolist())),
            "bse": dict(zip(X.columns, r.bse.tolist())),
            "pvalues": dict(zip(X.columns, r.pvalues.tolist())),
            "scale": float(r.scale),
        }
    except Exception as e:
        return {"model": "gee_exchangeable", "converged": False, "error": str(e)}


def fit_fixed_effects(df, preds):
    cc = df.groupby("country_iso3").size()
    eligible = cc[cc > 1].index.tolist()
    if len(eligible) < 2:
        return {"model": "glm_country_fixed_effects", "converged": False,
                "error": "insufficient countries with multiple years"}
    ds = df[df["country_iso3"].isin(eligible)].copy()
    dummies = pd.get_dummies(ds["country_iso3"], drop_first=True, prefix="ctry")
    X = _prepare_X(ds, preds)
    X = pd.concat([X, dummies], axis=1)
    y = _binomial_y(ds)
    try:
        m = sm.GLM(y, X.values, family=sm.families.Binomial())
        r = m.fit()
        return {
            "model": "glm_country_fixed_effects",
            "converged": True,
            "n_obs": int(r.nobs),
            "n_countries_fe": len(eligible),
            "params": dict(zip(X.columns, r.params.tolist())),
            "bse": dict(zip(X.columns, r.bse.tolist())),
            "pvalues": dict(zip(X.columns, r.pvalues.tolist())),
        }
    except Exception as e:
        return {"model": "glm_country_fixed_effects", "converged": False, "error": str(e)}


def within_country_models(df, preds):
    results = []
    cc = df.groupby("country_iso3").size()
    eligible = cc[cc >= 3].index.tolist()
    for ctry in eligible:
        dc = df[df["country_iso3"] == ctry].copy()
        X = _prepare_X(dc, preds)
        y = _binomial_y(dc)
        try:
            m = sm.GLM(y, X.values, family=sm.families.Binomial())
            r = m.fit()
            results.append({
                "country": ctry,
                "n_obs": len(dc),
                "converged": True,
                "params": dict(zip(X.columns, r.params.tolist())),
                "pvalues": dict(zip(X.columns, r.pvalues.tolist())),
            })
        except Exception as e:
            results.append({"country": ctry, "n_obs": len(dc),
                            "converged": False, "error": str(e)})
    return results


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_data(INPUT_PATH)
    print(f"Loaded {len(df)} country-year observations ({df['country_iso3'].nunique()} countries)")

    preds = ["dtp3_coverage", "reported_cases"]
    all_results = []

    # 1. Baseline GLM
    print("1. Baseline GLM...")
    all_results.append(fit_glm(df, preds, "binomial_glm_baseline"))

    # 2. Cluster-robust SEs
    print("2. Cluster-robust GLM...")
    all_results.append(fit_cluster_robust(df, preds))

    # 3. GEE
    print("3. GEE exchangeable...")
    all_results.append(fit_gee(df, preds))

    # 4. Removed invalid weighted GLM
    print("4. Removing invalid weighted GLM variant...")
    all_results.append(
        {
            "model": "glm_weighted_by_n_removed",
            "converged": False,
            "error": (
                "removed_invalid_weighting: grouped binomial responses already encode trial counts; "
                "additional normalized variance weights created pseudo-precision"
            ),
        }
    )

    # 5. Exclude n < 10
    print("5. Excluding n<10...")
    d10 = df[df["n_genomes_prn_interpretable"] >= 10].copy()
    if len(d10) >= 5:
        all_results.append(
            fit_glm(
                d10,
                preds,
                "glm_exclude_n_lt_10_country_cluster_or_hc1",
                cov_type="cluster" if d10["country_iso3"].nunique() >= 3 else "HC1",
                cov_kwds={"groups": d10["country_iso3"].astype(str).to_numpy(dtype=object)}
                if d10["country_iso3"].nunique() >= 3
                else None,
            )
        )
    else:
        all_results.append({"model": "glm_exclude_n_lt_10_country_cluster_or_hc1", "converged": False,
                            "error": f"only {len(d10)} obs after filter"})

    # 6. Fixed effects
    print("6. Country fixed effects...")
    all_results.append(fit_fixed_effects(df, preds))

    # 7. Extended with acellular
    print("7. Extended model (acellular covariate)...")
    de = df.dropna(subset=["acellular_vs_whole_cell"]).copy()
    if len(de) >= 5:
        de["is_acellular"] = (de["acellular_vs_whole_cell"] == "mixed_or_acellular").astype(int)
        epreds = ["dtp3_coverage", "reported_cases", "is_acellular"]
        all_results.append(fit_glm(de, epreds, "glm_extended_acellular"))

    # 8. Within-country
    print("8. Within-country models...")
    wc = within_country_models(df, preds)

    # --- Save coefficient table ---
    rows = []
    for res in all_results:
        mdl = res.get("model", "?")
        params = res.get("params", {})
        for pred, coef in params.items():
            row = {"model": mdl, "predictor": pred, "coef": coef,
                   "converged": res.get("converged", False)}
            bse = res.get("bse", {})
            pval = res.get("pvalues", {})
            bse_rob = res.get("bse_robust", {})
            pv_rob = res.get("pvalues_robust", {})
            cl = res.get("conf_low", {})
            ch = res.get("conf_high", {})
            if res.get("covariance_type"):
                row["covariance_type"] = res["covariance_type"]
            if pred in bse:
                row["bse"] = bse[pred]
            if pred in pval:
                row["pvalue"] = pval[pred]
            if pred in bse_rob:
                row["bse_robust"] = bse_rob[pred]
            if pred in pv_rob:
                row["pvalue_robust"] = pv_rob[pred]
            if pred in cl:
                row["conf_low"] = cl[pred]
            if pred in ch:
                row["conf_high"] = ch[pred]
            rows.append(row)

    out_df = pd.DataFrame(rows)
    out_path = OUTPUT_DIR / "bp_country_year_mixed_effects_coefficients.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nCoefficients → {out_path}")

    # --- Save diagnostics JSON ---
    diag = {
        "models": all_results,
        "within_country": wc,
        "data_summary": {
            "n_observations": len(df),
            "n_countries": int(df["country_iso3"].nunique()),
            "n_cells_n_ge_10": int((df["n_genomes_prn_interpretable"] >= 10).sum()),
            "n_cells_n_lt_10": int((df["n_genomes_prn_interpretable"] < 10).sum()),
            "country_counts": df.groupby("country_iso3").size().to_dict(),
        },
    }
    diag_path = OUTPUT_DIR / "bp_country_year_mixed_effects_diagnostics.json"
    with diag_path.open("w") as f:
        json.dump(diag, f, indent=2, default=str)
    print(f"Diagnostics → {diag_path}")

    # --- Print summary ---
    print("\n=== MODEL COMPARISON (DTP3 coverage coefficient) ===")
    for res in all_results:
        mdl = res["model"]
        ok = "OK" if res.get("converged") else "FAIL"
        params = res.get("params", {})
        pv = res.get("pvalues", {})
        pv_rob = res.get("pvalues_robust", {})
        dtp3 = str(params.get("dtp3_coverage", "—"))
        dtp3_p = str(pv.get("dtp3_coverage", "—"))
        dtp3_p_rob = str(pv_rob.get("dtp3_coverage", "—"))
        print(f"  {mdl:45s} [{ok}]  coef={dtp3:>10s}  p={str(dtp3_p):>10s}  p_robust={str(dtp3_p_rob):>10s}")

    print("\n=== WITHIN-COUNTRY MODELS ===")
    for wcr in wc:
        ok = "OK" if wcr.get("converged") else "FAIL"
        pars = wcr.get("params", {})
        pvs = wcr.get("pvalues", {})
        dtp3 = str(pars.get("dtp3_coverage", "—"))
        dtp3_p = str(pvs.get("dtp3_coverage", "—"))
        print(f"  {wcr['country']:5s} (n={wcr['n_obs']:2d}) [{ok}]  DTP3_coef={dtp3}  p={dtp3_p}")


if __name__ == "__main__":
    main()
