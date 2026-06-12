#!/usr/bin/env python3

import argparse
from pathlib import Path


def read_lines(path: Path) -> list[str]:
    return [l.strip() for l in path.read_text(encoding="utf-8").splitlines() if l.strip()]


def alt_accessions(acc: str) -> list[str]:
    acc = acc.strip()
    out: list[str] = [acc]
    if acc.startswith("GCF_"):
        out.append("GCA_" + acc[len("GCF_") :])
    elif acc.startswith("GCA_"):
        out.append("GCF_" + acc[len("GCA_") :])
    return list(dict.fromkeys(out))


def pick_fasta(assembly_dir: Path) -> Path | None:
    # Prefer *_genomic.fna (datasets layout), fall back to any .fna/.fna.gz
    cands: list[Path] = []
    for p in assembly_dir.rglob("*"):
        if not p.is_file():
            continue
        n = p.name
        if n.endswith("_genomic.fna"):
            cands.append(p)
        elif n.endswith(".fna") or n.endswith(".fna.gz"):
            cands.append(p)

    if not cands:
        return None
    cands.sort(key=lambda x: (0 if x.name.endswith("_genomic.fna") else 1, len(x.name)))
    return cands[0]


def main() -> None:
    ap = argparse.ArgumentParser(description="Resolve assembly accession -> local FASTA path")
    ap.add_argument("--accessions", required=True, help="Accession list (one per line)")
    ap.add_argument("--data-root", required=True, help="datasets extracted root: ncbi_dataset/data")
    ap.add_argument("--out-tsv", required=True, help="Output TSV with resolved FASTA paths")
    ap.add_argument("--out-missing", required=True, help="Output list of missing accessions")
    args = ap.parse_args()

    acc_path = Path(args.accessions)
    if not acc_path.exists() or acc_path.stat().st_size == 0:
        raise SystemExit(f"ERROR: accession list missing or empty: {acc_path}")

    data_root = Path(args.data_root)
    if not data_root.exists() or not data_root.is_dir():
        raise SystemExit(f"ERROR: data root not found: {data_root}")

    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    out_missing = Path(args.out_missing)
    out_missing.parent.mkdir(parents=True, exist_ok=True)

    accs = read_lines(acc_path)

    rows: list[dict[str, str]] = []
    missing: list[str] = []

    for acc in accs:
        resolved_dir: Path | None = None
        resolved_acc: str | None = None

        for a in alt_accessions(acc):
            d = data_root / a
            if d.exists() and d.is_dir():
                resolved_dir = d
                resolved_acc = a
                break

        if resolved_dir is None or resolved_acc is None:
            missing.append(acc)
            rows.append(
                {
                    "input_accession": acc,
                    "resolved_accession": "",
                    "status": "missing_dir",
                    "fasta_path": "",
                    "note": "no folder under data-root (tried GCA/GCF)",
                }
            )
            continue

        fasta = pick_fasta(resolved_dir)
        if fasta is None:
            rows.append(
                {
                    "input_accession": acc,
                    "resolved_accession": resolved_acc,
                    "status": "missing_fasta",
                    "fasta_path": "",
                    "note": "folder exists but no .fna found",
                }
            )
            continue

        rows.append(
            {
                "input_accession": acc,
                "resolved_accession": resolved_acc,
                "status": "ok",
                "fasta_path": str(fasta),
                "note": "",
            }
        )

    # Write TSV
    header = ["input_accession", "resolved_accession", "status", "fasta_path", "note"]
    with out_tsv.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r.get(h, "") for h in header) + "\n")

    out_missing.write_text("\n".join(missing) + ("\n" if missing else ""), encoding="utf-8")

    ok = sum(1 for r in rows if r["status"] == "ok")
    print("[Index] input accessions:", len(accs))
    print("[Index] ok FASTA:", ok)
    print("[Index] missing dirs:", sum(1 for r in rows if r["status"] == "missing_dir"))
    print("[Index] missing fasta:", sum(1 for r in rows if r["status"] == "missing_fasta"))
    print("Wrote:")
    print(f"  - {out_tsv}")
    print(f"  - {out_missing}")


if __name__ == "__main__":
    main()
