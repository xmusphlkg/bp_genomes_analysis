#!/usr/bin/env python3
"""Extract key statistics for manuscript."""

import os
import csv
import json
from pathlib import Path
from datetime import datetime

def main():
    print("=" * 70)
    print("EXTRACTING KEY STATISTICS FOR MANUSCRIPT")
    print("=" * 70)
    
    base_dir = Path(__file__).resolve().parents[3]
    project_data_root = Path(
        os.environ.get(
            "PERTUSSIS_PROJECT_DATA_ROOT",
            str(base_dir / "pertussis_data" / "pertussis_gene"),
        )
    )
    manuscript_dir = base_dir / "manuscript"
    # Try to find and load metadata
    metadata_candidates = [
        base_dir / "state" / "manifest" / "manifest.tsv",
        manuscript_dir / "figure_data" / "project_genome_metadata_manifest.tsv",
        project_data_root / "step1_ingest" / "bp_metadata_clean.csv",
    ]
    
    rows = None
    fieldnames = []
    loaded_source = None
    for candidate in metadata_candidates:
        if candidate.exists():
            print(f"Loading: {candidate.name}")
            try:
                delimiter = "\t" if candidate.suffix == ".tsv" else ","
                with open(candidate, newline="") as handle:
                    reader = csv.DictReader(handle, delimiter=delimiter)
                    rows = list(reader)
                    fieldnames = reader.fieldnames or []
                loaded_source = candidate
                print(f"✓ Loaded {len(rows):,} records")
                break
            except Exception as e:
                print(f"Error: {e}")
    
    stats = {"extraction_date": datetime.now().isoformat()}
    
    if rows is not None:
        stats["source_file"] = str(loaded_source.relative_to(base_dir))

        # Basic counts
        stats["total_genomes"] = len(rows)
        
        # Country counts
        country_field = "country_iso3" if "country_iso3" in fieldnames else "country"
        if country_field in fieldnames:
            invalid_country_values = {"", "na", "nan", "missing", "not provided"}
            countries = sorted(
                {
                    str(row.get(country_field, "")).strip()
                    for row in rows
                    if str(row.get(country_field, "")).strip().lower()
                    not in invalid_country_values
                }
            )
            stats["country_field"] = country_field
            stats["unique_countries"] = len(countries)
            stats["countries_list"] = countries
        
        # Temporal range
        year_field = None
        if "collection_year" in fieldnames:
            year_field = "collection_year"
            year_values = [row.get(year_field, "") for row in rows]
        elif "year" in fieldnames:
            year_field = "year"
            year_values = [row.get(year_field, "") for row in rows]
        elif "collection_date" in fieldnames:
            year_field = "collection_date"
            year_values = [str(row.get(year_field, ""))[:4] for row in rows]
        else:
            year_values = []

        years = []
        for value in year_values:
            try:
                years.append(int(float(value)))
            except (TypeError, ValueError):
                continue
        if len(years) > 0:
            stats["year_field"] = year_field
            stats["min_year"] = min(years)
            stats["max_year"] = max(years)
    
    # Save results
    json_out = manuscript_dir / "key_statistics.json"
    with open(json_out, "w") as f:
        json.dump(stats, f, indent=2, default=str)
    
    txt_out = manuscript_dir / "key_statistics.txt"
    with open(txt_out, "w") as f:
        f.write("KEY STATISTICS\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Extraction Date: {stats['extraction_date']}\n\n")
        f.write(f"Source: {stats.get('source_file', 'N/A')}\n")
        total_genomes = stats.get("total_genomes", "N/A")
        total_genomes_text = (
            f"{total_genomes:,}" if isinstance(total_genomes, int) else str(total_genomes)
        )
        f.write(f"Total Genomes: {total_genomes_text}\n")
        f.write(f"Countries/Territories: {stats.get('unique_countries', 'N/A')}\n")
        f.write(f"Time Span: {stats.get('min_year', 'N/A')} - {stats.get('max_year', 'N/A')}\n")
    
    print(f"\n✓ Saved: {json_out}")
    print(f"✓ Saved: {txt_out}")
    print("\n" + "=" * 70)

if __name__ == "__main__":
    main()
