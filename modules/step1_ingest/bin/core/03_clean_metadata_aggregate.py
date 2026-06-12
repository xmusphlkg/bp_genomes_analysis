import argparse
import re
from datetime import datetime
from pathlib import Path

import pandas as pd

def find_col(df: pd.DataFrame, patterns: list[str]) -> str | None:
    cols = list(df.columns)
    low = [c.lower() for c in cols]
    for p in patterns:
        p = p.lower()
        for c, cl in zip(cols, low):
            if p in cl:
                return c
    return None

def parse_date(s: str):
    if s is None:
        return (pd.NaT, pd.NA, pd.NA, pd.NA, "missing")
    s = str(s).strip()
    if s == "" or s.lower() in {"missing", "na", "n/a", "unknown"}:
        return (pd.NaT, pd.NA, pd.NA, pd.NA, "missing")
    s = s.split("T")[0].strip()

    if re.fullmatch(r"\d{4}", s):
        return (pd.NaT, int(s), pd.NA, pd.NA, "year")
    if re.fullmatch(r"\d{4}-\d{2}", s):
        y, m = s.split("-")
        return (pd.NaT, int(y), int(m), pd.NA, "year_month")
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", s):
        dt = pd.to_datetime(s, errors="coerce")
        if pd.isna(dt):
            return (pd.NaT, pd.NA, pd.NA, pd.NA, "unparsed")
        iso = dt.isocalendar()
        return (dt, int(dt.year), int(dt.month), int(iso.week), "full_date")

    # fallback
    dt = pd.to_datetime(s, errors="coerce", dayfirst=False)
    if pd.isna(dt):
        return (pd.NaT, pd.NA, pd.NA, pd.NA, "unparsed")
    iso = dt.isocalendar()
    return (dt, int(dt.year), int(dt.month), int(iso.week), "parsed_fallback")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--extended", required=True, help="Extended TSV produced by step2")
    ap.add_argument("--prefix", default="bp")
    ap.add_argument("--readme", default="readme.md")
    args = ap.parse_args()

    ext = Path(args.extended)
    if not ext.exists():
        raise SystemExit(f"ERROR: extended TSV not found: {ext}")

    df = pd.read_csv(ext, sep="\t", dtype=str)

    # detect columns from human-readable headers
    col_acc  = find_col(df, ["assembly accession", "accession"])
    col_geo  = find_col(df, ["biosample geographic location", "geographic location"])
    col_date = find_col(df, ["biosample collection date", "collection date"])

    if not col_acc or not col_geo or not col_date:
        raise SystemExit("ERROR: Could not detect accession/geo/date columns from TSV header.")

    total = len(df)
    geo_raw = df[col_geo].fillna("").astype(str).str.strip()
    date_raw = df[col_date].fillna("").astype(str).str.strip()
    missing_geo = int((geo_raw == "").sum())
    missing_date = int((date_raw == "").sum())

    # country
    df["geo_raw"] = df[col_geo]
    df["country"] = geo_raw.str.split(":").str[0].str.strip()
    df.loc[df["country"] == "", "country"] = pd.NA

    # date parsing
    df["collection_date_raw"] = date_raw
    parsed = df["collection_date_raw"].apply(parse_date)
    df["collection_date"] = [x[0] for x in parsed]
    df["year"] = [x[1] for x in parsed]
    df["month"] = [x[2] for x in parsed]
    df["iso_week"] = [x[3] for x in parsed]
    df["date_resolution"] = [x[4] for x in parsed]

    df["month_key"] = df.apply(
        lambda r: f"{int(r['year']):04d}-{int(r['month']):02d}"
        if pd.notna(r["year"]) and pd.notna(r["month"]) else pd.NA,
        axis=1,
    )
    df["week_key"] = df.apply(
        lambda r: f"{int(r['year']):04d}-W{int(r['iso_week']):02d}"
        if pd.notna(r["year"]) and pd.notna(r["iso_week"]) else pd.NA,
        axis=1,
    )

    # outputs
    df.to_csv(f"{args.prefix}_metadata_clean.csv", index=False)

    date_summary = df["date_resolution"].value_counts(dropna=False).reset_index()
    date_summary.columns = ["date_resolution", "count"]
    date_summary.to_csv(f"{args.prefix}_date_resolution_summary.csv", index=False)

    cy = (df.dropna(subset=["country", "year"])
            .assign(year=lambda x: x["year"].astype(int))
            .groupby(["country", "year"]).size().reset_index(name="n_genomes"))
    cy.to_csv(f"{args.prefix}_country_year_counts.csv", index=False)

    cm = (df.dropna(subset=["country", "month_key"])
            .groupby(["country", "month_key"]).size().reset_index(name="n_genomes"))
    cm.to_csv(f"{args.prefix}_country_month_counts.csv", index=False)

    cw = (df.dropna(subset=["country", "week_key"])
            .groupby(["country", "week_key"]).size().reset_index(name="n_genomes"))
    cw.to_csv(f"{args.prefix}_country_week_counts.csv", index=False)

    # month-ready subset (monthly usable)
    hi = df[df["date_resolution"].isin(["full_date", "year_month", "parsed_fallback"])].copy()
    hi.to_csv(f"{args.prefix}_samples_month_ready.csv", index=False)

    cm_hi = (hi.dropna(subset=["country", "month_key"])
               .groupby(["country", "month_key"]).size().reset_index(name="n_genomes"))
    cm_hi.to_csv(f"{args.prefix}_country_month_counts_hires.csv", index=False)

    # accession lists
    acc_month = hi[col_acc].dropna().drop_duplicates()
    acc_month.to_csv("assembly_accessions_month_ready.txt", index=False, header=False)

    acc_all = df[col_acc].dropna().drop_duplicates()
    acc_all.to_csv("assembly_accessions_all.txt", index=False, header=False)

    # summary tables for readme
    top_year = (cy.groupby("country")["n_genomes"].sum().sort_values(ascending=False).head(20))
    top_month = (cm_hi.groupby("country")["n_genomes"].sum().sort_values(ascending=False).head(20))

    year_min = int(cy["year"].min()) if not cy.empty else None
    year_max = int(cy["year"].max()) if not cy.empty else None

    # write readme.ms
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    readme = Path(args.readme)
    with readme.open("w", encoding="utf-8") as f:
        f.write(f"# Bordetella pertussis: NCBI metadata step1 (generated {now})\n\n")
        f.write("This folder contains step-1 outputs for *Bordetella pertussis* genome assembly metadata (NCBI Datasets CLI).\n\n")

        f.write("## Pipeline steps\n")
        f.write("1) Fetch genome assembly report (JSONL) from NCBI Datasets.\n")
        f.write("2) Export min/extended TSV (including BioSample collection date & geographic location).\n")
        f.write("3) Clean metadata: derive country/year/month_key/week_key; create count tables.\n\n")

        f.write("## Key statistics\n")
        f.write(f"- Total assemblies in report: {total}\n")
        f.write(f"- Missing geographic location (raw): {missing_geo} ({missing_geo/total:.2%})\n")
        f.write(f"- Missing collection date (raw): {missing_date} ({missing_date/total:.2%})\n\n")

        f.write("### Date resolution breakdown\n")
        for _, r in date_summary.iterrows():
            f.write(f"- {r['date_resolution']}: {int(r['count'])}\n")
        f.write("\n")
        if year_min is not None:
            f.write(f"- Year range (year available): {year_min}–{year_max}\n\n")

        f.write("## Modeling note (important)\n")
        f.write("- Many samples only have **year-level** dates. For global spread + national weekly/monthly cases,\n")
        f.write("  use a **two-track strategy**:\n")
        f.write("  1) Global analysis at **annual** resolution: use `*_country_year_counts.csv`.\n")
        f.write("  2) High-resolution (monthly/weekly) analysis on a subset with full dates: use `*_samples_month_ready.csv`.\n\n")

        f.write("## Top countries by genomes (year available)\n")
        for c, v in top_year.items():
            f.write(f"- {c}: {int(v)}\n")
        f.write("\n")

        f.write("## Top countries by genomes (month-ready subset)\n")
        for c, v in top_month.items():
            f.write(f"- {c}: {int(v)}\n")
        f.write("\n")

        f.write("## Output files\n")
        f.write(f"- `{args.prefix}_genome_report.jsonl`: assembly report (JSONL)\n")
        f.write(f"- `{args.prefix}_min_metadata.tsv`: small TSV for quick inspection\n")
        f.write(f"- `{args.prefix}_extended_metadata.tsv`: extended TSV with BioSample date/geo\n")
        f.write(f"- `{args.prefix}_metadata_clean.csv`: per-assembly cleaned metadata with derived time keys\n")
        f.write(f"- `{args.prefix}_date_resolution_summary.csv`: date resolution counts\n")
        f.write(f"- `{args.prefix}_country_year_counts.csv`: country×year genome counts\n")
        f.write(f"- `{args.prefix}_country_month_counts.csv`: country×month counts (where month exists)\n")
        f.write(f"- `{args.prefix}_country_week_counts.csv`: country×ISO-week counts (where full date exists)\n")
        f.write(f"- `{args.prefix}_samples_month_ready.csv`: subset usable for monthly analysis\n")
        f.write(f"- `{args.prefix}_country_month_counts_hires.csv`: month counts from month-ready subset\n")
        f.write("- `assembly_accessions_month_ready.txt`: accessions for downloading genomes (month-ready)\n")
        f.write("- `assembly_accessions_all.txt`: accessions for downloading genomes (all)\n")

        f.write("\n## Download genomes (optional)\n")
        f.write("To download genomes referenced by the accessions lists, run:\n\n")
        f.write("- `bash run_step1.sh --download-genomes`\n\n")
        f.write("### Speed tips\n")
        f.write("- Use parallel *datasets* downloads by splitting the accession list (recommended):\n")
        f.write("  - `PARALLEL=4 bash run_step1.sh --download-genomes`\n")
        f.write("- aria2 mode requires a URL list file (one URL per line). If `ARIA2_URLS` is not set, the pipeline will automatically fall back to `datasets` download.\n")
        f.write("  - Example: `USE_ARIA2=1 ARIA2_URLS=urls.txt ARIA2_JOBS=8 bash run_step1.sh --download-genomes`\n")

    print(f"Wrote: {readme}")
    print("Wrote core outputs:")
    print(f"  - {args.prefix}_metadata_clean.csv")
    print(f"  - {args.prefix}_country_year_counts.csv")
    print(f"  - {args.prefix}_country_month_counts.csv")
    print(f"  - {args.prefix}_country_week_counts.csv")
    print(f"  - {args.prefix}_samples_month_ready.csv")
    print("  - assembly_accessions_month_ready.txt")
    print("  - assembly_accessions_all.txt")

if __name__ == "__main__":
    main()
