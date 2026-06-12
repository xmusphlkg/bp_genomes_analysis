#!/usr/bin/env python3

import argparse
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
import zipfile


@dataclass
class CheckResult:
    ok: bool
    message: str


def read_lines(path: Path) -> list[str]:
    return [l.strip() for l in path.read_text(encoding="utf-8").splitlines() if l.strip()]


def is_safe_member(name: str) -> bool:
    # Prevent zip-slip: no absolute paths and no parent traversal.
    p = Path(name)
    if p.is_absolute():
        return False
    if any(part == ".." for part in p.parts):
        return False
    return True


def safe_extract_zip(zip_path: Path, outdir: Path) -> None:
    with zipfile.ZipFile(zip_path, "r") as z:
        for member in z.infolist():
            if not is_safe_member(member.filename):
                raise SystemExit(f"ERROR: unsafe zip member path: {member.filename}")
        z.extractall(outdir)


def verify_zip(zip_path: Path) -> CheckResult:
    if not zip_path.exists() or zip_path.stat().st_size == 0:
        return CheckResult(False, f"zip missing or empty: {zip_path}")

    try:
        with zipfile.ZipFile(zip_path, "r") as z:
            # Checks CRC for each member; returns first bad member name or None.
            bad = z.testzip()
            if bad is not None:
                return CheckResult(False, f"zip failed CRC check at member: {bad}")
    except zipfile.BadZipFile:
        return CheckResult(False, "zip is not a valid zip file (likely interrupted download)")
    except Exception as e:
        return CheckResult(False, f"zip validation error: {e}")

    return CheckResult(True, "zip OK")


def find_data_root(outdir: Path) -> Path | None:
    # datasets zip usually extracts into outdir/ncbi_dataset/data
    cand = outdir / "ncbi_dataset" / "data"
    if cand.exists() and cand.is_dir():
        return cand
    # Sometimes users point outdir at the extracted root already
    cand2 = outdir / "data"
    if cand2.exists() and cand2.is_dir():
        return cand2
    return None


def count_fasta_files(assembly_dir: Path) -> int:
    n = 0
    for p in assembly_dir.rglob("*"):
        if p.is_file() and (p.name.endswith(".fna") or p.name.endswith(".fna.gz")):
            n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Verify NCBI datasets download zip, extract safely, and validate extracted genome folders/files."
    )
    ap.add_argument("--zip", dest="zip_path", help="datasets zip file (optional; if provided, script can auto-extract)")
    ap.add_argument("--outdir", required=True, help="extraction directory")
    ap.add_argument("--accessions", help="accession list txt (one per line); used to validate extracted contents")
    ap.add_argument(
        "--force-extract",
        action="store_true",
        help="delete outdir contents first, then extract (requires --zip)",
    )
    ap.add_argument(
        "--no-extract",
        action="store_true",
        help="do not extract even if --zip is provided; only verify zip and/or outdir",
    )
    ap.add_argument("--max-missing-show", type=int, default=20, help="max missing accessions to print")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    zip_path = Path(args.zip_path) if args.zip_path else None

    if args.force_extract and zip_path is None:
        raise SystemExit("ERROR: --force-extract requires --zip")

    # 1) Zip checks + extraction (optional)
    if zip_path is not None:
        zr = verify_zip(zip_path)
        if not zr.ok:
            print(f"[Zip] FAIL: {zr.message}")
            print("[Hint] Delete the partial zip and re-download (datasets download does not resume).")
            print(f"       rm -f {zip_path}")
            sys.exit(2)
        print(f"[Zip] OK: {zip_path}")

        if not args.no_extract:
            if args.force_extract:
                for child in outdir.iterdir():
                    if child.is_dir():
                        shutil.rmtree(child)
                    else:
                        child.unlink()

            # Skip extraction if it looks already extracted
            data_root = find_data_root(outdir)
            if data_root is None or args.force_extract:
                print(f"[Extract] extracting to: {outdir}")
                safe_extract_zip(zip_path, outdir)
            else:
                print(f"[Extract] skip (already extracted): {data_root}")

    # 2) Outdir checks
    data_root = find_data_root(outdir)
    if data_root is None:
        print(f"[Outdir] FAIL: cannot find extracted data root under: {outdir}")
        print("[Hint] Expected: outdir/ncbi_dataset/data (datasets zip layout)")
        sys.exit(3)

    print(f"[Outdir] OK: data root = {data_root}")

    # 3) Accession checks (optional)
    if args.accessions:
        acc_path = Path(args.accessions)
        if not acc_path.exists() or acc_path.stat().st_size == 0:
            raise SystemExit(f"ERROR: accession list missing or empty: {acc_path}")

        accs = read_lines(acc_path)
        missing: list[str] = []
        no_fasta: list[str] = []

        for acc in accs:
            d = data_root / acc
            if not d.exists():
                missing.append(acc)
                continue
            if count_fasta_files(d) == 0:
                no_fasta.append(acc)

        ok = (len(missing) == 0) and (len(no_fasta) == 0)
        print(f"[Check] accessions: {len(accs)}")
        print(f"[Check] missing dirs: {len(missing)}")
        print(f"[Check] no FASTA: {len(no_fasta)}")

        if missing:
            show = missing[: args.max_missing_show]
            print("[Missing] " + ", ".join(show) + (" ..." if len(missing) > len(show) else ""))
        if no_fasta:
            show = no_fasta[: args.max_missing_show]
            print("[NoFASTA] " + ", ".join(show) + (" ..." if len(no_fasta) > len(show) else ""))

        if not ok:
            sys.exit(4)

    print("[Done] verify OK")


if __name__ == "__main__":
    main()
