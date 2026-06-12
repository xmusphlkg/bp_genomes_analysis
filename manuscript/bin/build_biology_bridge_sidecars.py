#!/usr/bin/env python3
"""Build manuscript-facing biology bridge sidecars.

The tables generated here are deliberately small. They connect the
genome-defined prn disruption layer to phenotype-linked published
surveillance, post-pandemic counterexamples, and AMR/lineage context
without redefining the paper's primary estimands.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FIGURE_DATA = ROOT / "manuscript" / "figure_data"
SUPPLEMENTARY = ROOT / "manuscript" / "supplementary"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def truthy(value: str) -> bool:
    return str(value).strip().lower() == "true"


def frac(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return ""
    return f"{numerator / denominator:.6f}"


def count_group(rows: list[dict[str, str]], group_variable: str, group_value: str) -> dict[str, object]:
    selected = [row for row in rows if (row.get(group_variable) or "unassigned") == group_value]
    interpretable = [row for row in selected if truthy(row.get("prn_interpretable", ""))]
    disrupted = [row for row in interpretable if truthy(row.get("prn_disrupted", ""))]
    gap1043 = [row for row in disrupted if "gap1043" in row.get("prn_event_id", "")]
    ptxp3 = [row for row in selected if row.get("ptxP_label") in {"ptxP_3", "ptxP3"}]
    fim31 = [row for row in selected if row.get("fim3_label") in {"fim3_1", "fim3-1"}]
    return {
        "analysis_layer": "internal_genome_context",
        "group_variable": group_variable,
        "group_value": group_value,
        "n_records": len(selected),
        "n_prn_interpretable": len(interpretable),
        "n_prn_disrupted": len(disrupted),
        "prn_disrupted_fraction_among_interpretable": frac(len(disrupted), len(interpretable)),
        "n_gap1043": len(gap1043),
        "n_ptxP3": len(ptxp3),
        "n_fim3_1": len(fim31),
        "interpretation": "",
    }


def build_internal_bridge() -> list[dict[str, object]]:
    annotation = read_tsv(FIGURE_DATA / "published_overlap_annotation.tsv")
    concordance = read_tsv(FIGURE_DATA / "published_overlap_concordance.tsv")

    rows: list[dict[str, object]] = []
    for item in concordance:
        if item["summary_level"] == "overall":
            rows.append(
                {
                    "analysis_layer": "published_overlap",
                    "group_variable": item["metric_name"],
                    "group_value": "overall",
                    "n_records": item["n_overlap_rows"],
                    "n_prn_interpretable": "",
                    "n_prn_disrupted": "",
                    "prn_disrupted_fraction_among_interpretable": "",
                    "n_gap1043": "",
                    "n_ptxP3": "",
                    "n_fim3_1": "",
                    "interpretation": (
                        f"{item['n_concordant']} of {item['n_compared_rows']} compared rows "
                        f"were concordant (fraction {float(item['concordance_fraction']):.3f}). "
                        f"{item.get('notes', '')}".strip()
                    ),
                }
            )

    overlap_counter = Counter(row.get("published_overlap_found") for row in annotation)
    rows.append(
        {
            "analysis_layer": "published_overlap",
            "group_variable": "published_PRN_status_rows",
            "group_value": "overall",
            "n_records": overlap_counter["True"],
            "interpretation": (
                "Rows with a published PRN-status annotation; concordance is restricted to "
                "the compared subset rather than all annotated rows."
            ),
        }
    )

    for value in ["MS", "MR_A2047G", "other_23S_allele", "unassigned"]:
        group = count_group(annotation, "marker_23s_status", value)
        if int(group["n_records"]) > 0:
            group["interpretation"] = (
                "Internal 23S-marker stratum; used as support-only AMR context and not as a "
                "primary causal adjustment."
            )
            rows.append(group)

    for value in [
        "United States Sublineage 1",
        "ptxP3/fim3-1",
        "ptxP3/fim3-2",
        "MR-MT28-others",
        "MR-MT28-PG1",
        "MR-MT28-PG2",
        "MS-MT28",
        "unassigned",
    ]:
        group = count_group(annotation, "published_sublineage", value)
        if int(group["n_records"]) > 0:
            group["interpretation"] = (
                "Published-lineage stratum; used to check whether prn disruption is simply "
                "a single AMR or lineage label."
            )
            rows.append(group)

    return rows


def build_external_context() -> list[dict[str, object]]:
    return [
        {
            "context_id": "japan_otsuka_2012_expression_bridge",
            "country_or_region": "Japan",
            "period": "1990-2010 isolate collection",
            "evidence_type": "PRN immunoblot plus prn mechanism sequencing",
            "n_total_or_frame": 121,
            "n_prn_producing": 88,
            "n_prn_deficient": 33,
            "n_amr_or_a2047g": "",
            "headline": "Thirty-three of 121 isolates were PRN-negative; 24 carried an 84-bp signal-sequence deletion and nine carried IS481 insertions with ACTAGG direct repeats.",
            "interpretive_role": "phenotype-genotype bridge for deletion and IS481 lesion classes",
            "reference": "Otsuka et al. PLoS ONE 2012;7:e31985. doi:10.1371/journal.pone.0031985",
        },
        {
            "context_id": "usa_pawloski_2014_expression_bridge",
            "country_or_region": "United States",
            "period": "1935-2012 surveillance isolate collection",
            "evidence_type": "Western blot plus prn sequencing",
            "n_total_or_frame": 1300,
            "n_prn_producing": "",
            "n_prn_deficient": 306,
            "n_amr_or_a2047g": "",
            "headline": "PRN deficiency was confirmed in 306 isolates, including 276 prn::IS481 isolates, and exceeded half of tested isolates by 2012.",
            "interpretive_role": "phenotype-genotype bridge for recurrent IS481-mediated PRN loss",
            "reference": "Pawloski et al. Clin. Vaccine Immunol. 2014;21:119-125. doi:10.1128/CVI.00717-13",
        },
        {
            "context_id": "australia_lam_2014_expression_bridge",
            "country_or_region": "Australia",
            "period": "2008-2012 rise in PRN-negative isolates",
            "evidence_type": "Western immunoblot plus prn/promoter mechanism assessment",
            "n_total_or_frame": 320,
            "n_prn_producing": 224,
            "n_prn_deficient": 96,
            "n_amr_or_a2047g": "",
            "headline": "Ninety-six of 320 isolates were PRN-negative; among 80 mechanism-resolved PRN-negative isolates, 77 had insertion-sequence lesions, while 16 PRN-negative isolates lacked an obvious prn/promoter change.",
            "interpretive_role": "phenotype-genotype bridge and caution that protein testing can reveal lesions missed by locus inspection",
            "reference": "Lam et al. Emerg. Infect. Dis. 2014;20:626-633. doi:10.3201/eid2004.131478",
        },
        {
            "context_id": "europe_eupert_1998_2015_expression_bridge",
            "country_or_region": "Europe",
            "period": "EUpert I-IV, 1998-2015",
            "evidence_type": "ELISA antigen-expression testing across European isolate collections",
            "n_total_or_frame": 661,
            "n_prn_producing": "",
            "n_prn_deficient": 66,
            "n_amr_or_a2047g": "",
            "headline": "PRN-deficient frequency increased from 1.0% in 1998-2001 to 24.9% in 2012-2015; IS481 and promoter inversion were the largest mechanism classes.",
            "interpretive_role": "multi-country phenotype-genotype bridge for recurrent PRN loss",
            "reference": "Barkoff et al. Euro Surveill. 2019;24:1700832. doi:10.2807/1560-7917.ES.2019.24.7.1700832",
        },
        {
            "context_id": "france_2023_2024_resurgence",
            "country_or_region": "France",
            "period": "June 2023-April 2024",
            "evidence_type": "national-reference-centre isolate phenotype and WGS",
            "n_total_or_frame": 67,
            "n_prn_producing": 66,
            "n_prn_deficient": "",
            "n_amr_or_a2047g": 1,
            "headline": "Post-pandemic resurgence was dominated by pertactin-producing isolates despite one macrolide-resistant isolate.",
            "interpretive_role": "post-pandemic counterexample to a one-way PRN-loss resurgence model",
            "reference": "Rodrigues et al. Euro Surveill. 2024;29:2400459. doi:10.2807/1560-7917.ES.2024.29.31.2400459",
        },
        {
            "context_id": "belgium_2022_2023_expression",
            "country_or_region": "Belgium",
            "period": "2022-2023, after pre-pandemic PRN-negative peak",
            "evidence_type": "national-reference-centre PRN expression testing and WGS",
            "n_total_or_frame": 94,
            "n_prn_producing": 94,
            "n_prn_deficient": 0,
            "n_amr_or_a2047g": "",
            "headline": "All 94 tested post-COVID isolates expressed pertactin after a pre-pandemic PRN-negative peak around 65%.",
            "interpretive_role": "post-pandemic counterexample and phenotype bridge",
            "reference": "Martini et al. Microbiol. Spectr. 2026;14:e01535-25. doi:10.1128/spectrum.01535-25",
        },
        {
            "context_id": "belgium_2000_2023_expression_bridge",
            "country_or_region": "Belgium",
            "period": "2000-2023 national-reference-centre collection",
            "evidence_type": "PRN ELISA plus WGS mechanism assessment",
            "n_total_or_frame": 432,
            "n_prn_producing": 252,
            "n_prn_deficient": 180,
            "n_amr_or_a2047g": "",
            "headline": "Among 432 isolates tested by PRN ELISA, 180 were PRN-negative; WGS on 416 isolates linked loss to IS481, deletion, premature stop, promoter inversion and promoter deletion classes.",
            "interpretive_role": "phenotype-genotype bridge and post-COVID boundary",
            "reference": "Martini et al. Microbiol. Spectr. 2026;14:e01535-25. doi:10.1128/spectrum.01535-25",
        },
        {
            "context_id": "netherlands_2023_2024_prn_positive_boundary",
            "country_or_region": "Netherlands",
            "period": "2023-2024 resurgence",
            "evidence_type": "official national surveillance report",
            "n_total_or_frame": "",
            "n_prn_producing": "",
            "n_prn_deficient": "",
            "n_amr_or_a2047g": 0,
            "headline": "The reported PRN-deficient fraction fell to 7% in 2023-2024 after 13% in 2015-2017 and 25% in 2018-2020; no 23S rRNA A2047G macrolide-resistance mutation was detected.",
            "interpretive_role": "post-pandemic counterexample and AMR-negative boundary condition",
            "reference": "RIVM. Resurgence of Bordetella pertussis in the Netherlands in 2023-2024. Report 2025-0092. https://www.rivm.nl/bibliotheek/rapporten/2025-0092.pdf",
        },
        {
            "context_id": "slovenia_2002_2020_expression",
            "country_or_region": "Slovenia",
            "period": "2002-2020; high-frequency period 2017-2020",
            "evidence_type": "ELISA antigen-expression testing plus prn WGS",
            "n_total_or_frame": 123,
            "n_prn_producing": 75,
            "n_prn_deficient": 48,
            "n_amr_or_a2047g": "",
            "headline": "Overall 48 of 123 isolates were PRN-deficient; 44 of 49 isolates from 2017-2020 did not express PRN.",
            "interpretive_role": "direct PRN-expression bridge for genome-defined prn lesions",
            "reference": "Barkoff et al. Emerg. Infect. Dis. 2024;30:2429-2432. doi:10.3201/eid3011.231393",
        },
        {
            "context_id": "slovenia_2024_resurgence_boundary",
            "country_or_region": "Slovenia",
            "period": "2024 epidemic",
            "evidence_type": "national epidemic isolate report",
            "n_total_or_frame": "",
            "n_prn_producing": "",
            "n_prn_deficient": 0,
            "n_amr_or_a2047g": 0,
            "headline": "The 2024 epidemic report found no pertactin-deficient or macrolide-resistant isolates.",
            "interpretive_role": "post-pandemic counterexample after an earlier high PRN-deficient period",
            "reference": "Kastrin et al. Clin. Microbiol. Infect. 2025. doi:10.1016/j.cmi.2025.04.016",
        },
        {
            "context_id": "beijing_2022_2023_amr_antigen_deficiency",
            "country_or_region": "Beijing, China",
            "period": "2022-2023",
            "evidence_type": "sentinel-hospital isolate WGS and antimicrobial susceptibility",
            "n_total_or_frame": 44,
            "n_prn_producing": 14,
            "n_prn_deficient": 30,
            "n_amr_or_a2047g": 44,
            "headline": "All isolates were erythromycin resistant with the same 23S mutation background; 30 were Prn deficient and two were Fha deficient.",
            "interpretive_role": "AMR and antigen-deficiency boundary condition, not a direct vaccine-programme contrast",
            "reference": "Zhou et al. China CDC Wkly. 2024;6:437-441. doi:10.46234/ccdcw2024.085",
        },
        {
            "context_id": "shanghai_2018_2024_mr_mt28",
            "country_or_region": "Shanghai, China",
            "period": "2018-2024",
            "evidence_type": "active-surveillance isolate WGS and antimicrobial susceptibility",
            "n_total_or_frame": 178,
            "n_prn_producing": "",
            "n_prn_deficient": "",
            "n_amr_or_a2047g": "",
            "headline": "MT28 increased from 16% before 2020 to 61.17% after 2020; post-2020 MT28 isolates were macrolide resistant, ptxP3-associated and carried prn150 antigen profiles.",
            "interpretive_role": "AMR and lineage-background context, especially for China-dominated post-pandemic signals",
            "reference": "Xu et al. J. Clin. Microbiol. 2025;63:e0106425. doi:10.1128/jcm.01064-25",
        },
        {
            "context_id": "global_mr_mt28_post_pandemic",
            "country_or_region": "Global, enriched for China and post-pandemic genomes",
            "period": "pre-, during and post-COVID-19 comparison",
            "evidence_type": "global genomic surveillance and AMR/antigen-profile analysis",
            "n_total_or_frame": 8117,
            "n_prn_producing": "",
            "n_prn_deficient": "",
            "n_amr_or_a2047g": "",
            "headline": "MR-MT28 expanded in eastern China and spread internationally; most non-Chinese MR-MT28 isolates were PRN-deficient.",
            "interpretive_role": "lineage and AMR context for interpreting post-pandemic public-genome signals",
            "reference": "Zhang et al. J. Infect. 2026;92:106718. doi:10.1016/j.jinf.2026.106718",
        },
        {
            "context_id": "eu_eea_2024_resurgence",
            "country_or_region": "EU/EEA",
            "period": "2024",
            "evidence_type": "regional epidemiological surveillance report",
            "n_total_or_frame": 209674,
            "n_prn_producing": "",
            "n_prn_deficient": "",
            "n_amr_or_a2047g": "",
            "headline": "EU/EEA reported 209,674 pertussis cases in 2024, with resurgence attributed to multiple immunity, behaviour and testing factors.",
            "interpretive_role": "epidemiological boundary showing that resurgence is not reducible to PRN loss",
            "reference": "ECDC. Annual epidemiological report for 2024: Pertussis. 2026.",
        },
    ]


def main() -> None:
    internal_fields = [
        "analysis_layer",
        "group_variable",
        "group_value",
        "n_records",
        "n_prn_interpretable",
        "n_prn_disrupted",
        "prn_disrupted_fraction_among_interpretable",
        "n_gap1043",
        "n_ptxP3",
        "n_fim3_1",
        "interpretation",
    ]
    external_fields = [
        "context_id",
        "country_or_region",
        "period",
        "evidence_type",
        "n_total_or_frame",
        "n_prn_producing",
        "n_prn_deficient",
        "n_amr_or_a2047g",
        "headline",
        "interpretive_role",
        "reference",
    ]

    internal_rows = build_internal_bridge()
    external_rows = build_external_context()

    write_tsv(FIGURE_DATA / "biology_bridge_internal_sensitivity.tsv", internal_rows, internal_fields)
    write_tsv(
        SUPPLEMENTARY / "Supplementary_Table_60_Biology_Bridge_Internal_Sensitivity.tsv",
        internal_rows,
        internal_fields,
    )
    write_tsv(FIGURE_DATA / "biology_bridge_external_context.tsv", external_rows, external_fields)
    write_tsv(
        SUPPLEMENTARY / "Supplementary_Table_61_External_Phenotype_AMR_Context.tsv",
        external_rows,
        external_fields,
    )


if __name__ == "__main__":
    main()
