#!/usr/bin/env python3

import argparse
import math
import shutil
import subprocess
import zipfile
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


def remove_if_exists(path: Path) -> None:
    try:
        path.unlink()
    except FileNotFoundError:
        return
    except OSError:
        return


def read_lines(path: Path) -> list[str]:
    return [l.strip() for l in path.read_text(encoding="utf-8").splitlines() if l.strip()]


def find_data_root(outdir: Path) -> Path | None:
    cand = outdir / "ncbi_dataset" / "data"
    if cand.exists() and cand.is_dir():
        return cand
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


def alt_accessions(acc: str) -> list[str]:
    acc = acc.strip()
    out: list[str] = [acc]
    if acc.startswith("GCF_"):
        out.append("GCA_" + acc[len("GCF_") :])
    elif acc.startswith("GCA_"):
        out.append("GCF_" + acc[len("GCA_") :])
    return list(dict.fromkeys(out))


def already_extracted_accessions(outdir: Path, accessions: list[str]) -> set[str]:
    data_root = find_data_root(outdir)
    if data_root is None:
        return set()

    present: set[str] = set()
    for acc in accessions:
        for a in alt_accessions(acc):
            d = data_root / a
            if d.exists() and d.is_dir() and count_fasta_files(d) > 0:
                present.add(acc)
                break
    return present


@dataclass
class ZipCheck:
    ok: bool
    message: str


def verify_zip(zip_path: Path) -> ZipCheck:
    if not zip_path.exists() or zip_path.stat().st_size == 0:
        return ZipCheck(False, f"zip missing or empty: {zip_path}")

    try:
        with zipfile.ZipFile(zip_path, "r") as z:
            bad = z.testzip()
            if bad is not None:
                return ZipCheck(False, f"zip failed CRC check at member: {bad}")
    except zipfile.BadZipFile:
        return ZipCheck(False, "zip is not a valid zip file (likely interrupted download)")
    except Exception as e:
        return ZipCheck(False, f"zip validation error: {e}")

    return ZipCheck(True, "zip OK")


def is_safe_member(name: str) -> bool:
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


def safe_extract_zip_best_effort(zip_path: Path, outdir: Path) -> tuple[int, list[str]]:
    """Extract as much as possible from a potentially-corrupt zip.

    Returns (failures_count, failed_members).
    """
    failures = 0
    failed: list[str] = []
    with zipfile.ZipFile(zip_path, "r") as z:
        for member in z.infolist():
            if not is_safe_member(member.filename):
                failures += 1
                failed.append(member.filename)
                continue
            try:
                z.extract(member, outdir)
            except (zipfile.BadZipFile, RuntimeError, OSError, zlib.error) as e:
                failures += 1
                failed.append(member.filename)
                print(f"[Warn] best-effort extract skipped bad member: {member.filename} ({e})")
    return failures, failed


def split_lines(lines: list[str], parts: int | None, chunk_size: int | None) -> list[list[str]]:
    if not lines:
        return []

    if chunk_size is not None and chunk_size > 0:
        return [lines[i : i + chunk_size] for i in range(0, len(lines), chunk_size)]

    if parts is None or parts <= 1:
        return [lines]

    part_size = max(1, math.ceil(len(lines) / parts))
    return [lines[i * part_size : (i + 1) * part_size] for i in range(parts) if lines[i * part_size : (i + 1) * part_size]]


def run_datasets_download(
    inputfile: Path,
    zip_out: Path,
    include: str,
    fast_zip_validation: bool,
    stream: bool,
) -> int:
    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        "--inputfile",
        str(inputfile),
        "--include",
        include,
        "--filename",
        str(zip_out),
    ]
    if fast_zip_validation:
        cmd.append("--fast-zip-validation")

    if stream:
        return subprocess.run(cmd).returncode

    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise SystemExit(f"ERROR downloading genomes:\n{p.stderr.strip()}")
    return 0


def aria2_download(urls_file: Path, outdir: Path, aria2_opts: str, jobs: int) -> None:
    aria2c = shutil.which("aria2c")
    if aria2c is None:
        raise SystemExit("ERROR: aria2c not found in PATH. Install aria2 to use this mode.")

    cmd = [aria2c, "-i", str(urls_file), "-d", str(outdir), "-j", str(jobs)]
    if aria2_opts:
        cmd += aria2_opts.split()

    p = subprocess.run(cmd)
    if p.returncode != 0:
        raise SystemExit(f"aria2 download failed (exit {p.returncode})")


def write_chunk_file(lines: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def zip_base_name(zip_arg: str) -> str:
    p = Path(zip_arg)
    if p.suffix.lower() == ".zip":
        return p.stem
    return p.name


def download_and_extract_chunk(
    chunk_idx: int,
    chunk_file: Path,
    zip_path: Path,
    outdir: Path,
    include: str,
    fast_zip_validation: bool,
    retries: int,
    resume: bool,
    cleanup_zips: bool,
) -> None:
    # Resume shortcut
    if resume and zip_path.exists() and zip_path.stat().st_size > 0:
        zc = verify_zip(zip_path)
        if zc.ok:
            print(f"[Resume] chunk {chunk_idx}: zip OK, skip download: {zip_path.name}")
        else:
            print(f"[Resume] chunk {chunk_idx}: zip invalid, re-download: {zip_path.name} ({zc.message})")
            remove_if_exists(zip_path)

    attempts = max(1, retries + 1)
    for attempt in range(1, attempts + 1):
        if not zip_path.exists():
            try:
                rc = run_datasets_download(
                    chunk_file,
                    zip_path,
                    include=include,
                    fast_zip_validation=fast_zip_validation,
                    stream=True,
                )
            except KeyboardInterrupt:
                remove_if_exists(zip_path)
                raise
            if rc != 0:
                remove_if_exists(zip_path)
                if attempt < attempts:
                    print(
                        f"[Warn] chunk {chunk_idx}: datasets download failed (exit {rc}); retrying ({attempt}/{attempts})"
                    )
                    continue
                raise SystemExit(f"ERROR: chunk {chunk_idx} download failed (exit {rc})")

        zc = verify_zip(zip_path)
        if not zc.ok:
            if attempt < attempts:
                remove_if_exists(zip_path)
                print(f"[Warn] chunk {chunk_idx}: {zc.message}; retrying ({attempt}/{attempts})")
                continue

            # Last attempt: salvage what we can and continue.
            print(f"[Warn] chunk {chunk_idx}: {zc.message}; attempting best-effort extraction")
            failures, _ = safe_extract_zip_best_effort(zip_path, outdir)
            if failures:
                print(f"[Warn] chunk {chunk_idx}: best-effort extraction had {failures} skipped member(s)")

            # Preserve corrupt zip for inspection, but do not reuse on resume.
            corrupt = zip_path.with_name(zip_path.stem + ".corrupt.zip")
            try:
                if corrupt.exists():
                    corrupt.unlink()
                zip_path.rename(corrupt)
            except OSError:
                # Fall back to delete; next run will re-download.
                remove_if_exists(zip_path)

            return

        # Extract
        safe_extract_zip(zip_path, outdir)
        print(f"[OK] chunk {chunk_idx}: extracted {zip_path.name}")

        if cleanup_zips:
            remove_if_exists(zip_path)
        return


def ensure_downloads_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def chunk_workspace(downloads_dir: Path, base: str, chunk_size: int) -> Path:
    # Keep chunk manifests/zips in a stable subfolder to avoid mixing strategies.
    d = downloads_dir / f"{base}.chunks_{chunk_size}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def chunk_manifest_paths(chunks_dir: Path, base: str, chunk_idx: int) -> tuple[Path, Path]:
    chunk_file = chunks_dir / f"{base}.part{chunk_idx:04d}.accessions.txt"
    zip_path = chunks_dir / f"{base}.part{chunk_idx:04d}.zip"
    return chunk_file, zip_path


def load_or_create_chunk_manifests(
    chunks_dir: Path,
    base: str,
    accessions: list[str],
    chunk_size: int,
) -> list[tuple[int, Path, Path]]:
    # Stable chunking: part0001 contains accessions[0:chunk_size], etc.
    chunks = split_lines(accessions, parts=None, chunk_size=chunk_size)
    jobs: list[tuple[int, Path, Path]] = []
    for idx, lines in enumerate(chunks, start=1):
        chunk_file, zip_path = chunk_manifest_paths(chunks_dir, base, idx)
        if not chunk_file.exists() or chunk_file.stat().st_size == 0:
            write_chunk_file(lines, chunk_file)
        jobs.append((idx, chunk_file, zip_path))
    return jobs


def chunk_needs_work(outdir: Path, chunk_file: Path, resume: bool) -> bool:
    if not resume:
        return True
    lines = read_lines(chunk_file)
    if not lines:
        return False
    present = already_extracted_accessions(outdir, lines)
    return len(present) < len(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", required=True, help="accession list txt (one per line)")
    ap.add_argument("--zip", required=True, help="output zip filename (or prefix when --parallel > 1)")
    ap.add_argument("--outdir", required=True, help="unzip directory")

    ap.add_argument("--parallel", type=int, default=1, help="number of concurrent datasets downloads")
    ap.add_argument("--retries", type=int, default=0, help="retry datasets download failures N times (no extract retries)")
    ap.add_argument(
        "--no-fast-zip-validation",
        action="store_true",
        help="do NOT pass --fast-zip-validation to datasets",
    )

    ap.add_argument(
        "--resume",
        action="store_true",
        help="resume by skipping already extracted accessions and reusing valid part zip files",
    )
    ap.add_argument(
        "--chunk-size",
        type=int,
        default=200,
        help="accessions per part zip (smaller chunks reduce restart cost)",
    )
    ap.add_argument(
        "--downloads-dir",
        default=None,
        help="where to store part zip files and chunk lists (default: <outdir>/_downloads)",
    )
    ap.add_argument(
        "--cleanup-zips",
        action="store_true",
        help="delete part zip files after successful extraction (disables zip-based resume)",
    )

    ap.add_argument("--use-aria2", action="store_true", help="use aria2c with provided --urls file")
    ap.add_argument("--urls", help="file with one http/ftp url per line (required if --use-aria2)")
    ap.add_argument("--aria2-opts", default="-x16 -s16", help="extra options passed to aria2c (string)")
    ap.add_argument("--aria2-jobs", type=int, default=4, help="concurrent file jobs for aria2c (-j)")

    args = ap.parse_args()

    acc = Path(args.accessions)
    if not acc.exists() or acc.stat().st_size == 0:
        raise SystemExit(f"ERROR: accession list missing or empty: {acc}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.use_aria2:
        if not args.urls:
            raise SystemExit("ERROR: --urls must be provided when --use-aria2 is used")
        urls_file = Path(args.urls)
        if not urls_file.exists() or urls_file.stat().st_size == 0:
            raise SystemExit(f"ERROR: urls file missing or empty: {urls_file}")
        print("Using aria2c to download URLs...")
        aria2_download(urls_file, outdir, args.aria2_opts, int(args.aria2_jobs))
        print("aria2 downloads completed.")
        return

    parallel = max(1, int(args.parallel))
    fast_zip_validation = not args.no_fast_zip_validation
    retries = max(0, int(args.retries))

    accessions = read_lines(acc)
    if not accessions:
        raise SystemExit(f"ERROR: accession list missing or empty: {acc}")

    resume = bool(args.resume)
    if resume:
        present = already_extracted_accessions(outdir, accessions)
        print(f"[Resume] already extracted: {len(present)} / {len(accessions)}")

    base = zip_base_name(args.zip)
    downloads_dir = Path(args.downloads_dir) if args.downloads_dir else (outdir / "_downloads")
    downloads_dir = ensure_downloads_dir(downloads_dir)
    chunks_dir = chunk_workspace(downloads_dir, base, int(args.chunk_size))

    all_chunk_jobs = load_or_create_chunk_manifests(chunks_dir, base, accessions, int(args.chunk_size))
    chunk_jobs = [(i, cf, zp) for (i, cf, zp) in all_chunk_jobs if chunk_needs_work(outdir, cf, resume)]

    print(f"[Plan] total chunks: {len(all_chunk_jobs)} (chunk_size={int(args.chunk_size)})")
    print(f"[Plan] chunks needing work: {len(chunk_jobs)}")
    print(f"[Plan] downloads dir: {chunks_dir}")

    if not chunk_jobs:
        print("[Done] nothing to download (all accessions already extracted)")
        return

    include = "genome,seq-report"
    cleanup_zips = bool(args.cleanup_zips)

    # Run with limited concurrency
    try:
        with ThreadPoolExecutor(max_workers=parallel) as ex:
            futures = []
            for chunk_idx, chunk_file, zip_path in chunk_jobs:
                futures.append(
                    ex.submit(
                        download_and_extract_chunk,
                        chunk_idx,
                        chunk_file,
                        zip_path,
                        outdir,
                        include,
                        fast_zip_validation,
                        retries,
                        resume,
                        cleanup_zips,
                    )
                )

            for fut in as_completed(futures):
                fut.result()
    except KeyboardInterrupt:
        print("[Interrupted] stopping. Re-run with --resume to continue.")
        raise SystemExit(130)

    print(f"[Done] downloaded+extracted chunks: {len(chunk_jobs)}")
    if not cleanup_zips:
        print(f"[Note] part zips kept for resume under: {chunks_dir}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        raise SystemExit("Interrupted (Ctrl-C)")
