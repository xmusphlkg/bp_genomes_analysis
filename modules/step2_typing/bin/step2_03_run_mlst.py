#!/usr/bin/env python3

import argparse
import concurrent.futures as cf
import os
import subprocess
from pathlib import Path
from typing import Optional


def read_genome_paths(tsv_path: Path) -> list[str]:
    paths: list[str] = []
    seen: set[str] = set()
    with tsv_path.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            idx_status = header.index("status")
            idx_path = header.index("fasta_path")
        except ValueError:
            raise SystemExit("ERROR: genome_paths TSV missing required columns: status, fasta_path")

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_status, idx_path):
                continue
            if parts[idx_status] != "ok":
                continue
            p = parts[idx_path].strip()
            if not p:
                continue
            if p in seen:
                continue
            seen.add(p)
            paths.append(p)
    return paths


def read_completed(out_path: Path) -> set[str]:
    if not out_path.exists() or out_path.stat().st_size == 0:
        return set()
    done: set[str] = set()
    with out_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                # Keep resume robust: also treat our own error marker lines as completed.
                # Format: "# ERROR\t<FASTA_PATH>\t..."
                if line.startswith("# ERROR\t"):
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2 and parts[1].strip():
                        done.add(parts[1].strip())
                continue
            first = line.split("\t", 1)[0].strip()
            if first:
                done.add(first)
    return done


def build_env() -> dict[str, str]:
    """Return environment for running mlst.

    We also extend PERL5LIB to handle conda perl 'site_perl' version directory mismatches.
    """
    env = os.environ.copy()
    conda_prefix = env.get("CONDA_PREFIX")
    if conda_prefix:
        site = Path(conda_prefix) / "lib" / "perl5" / "site_perl"
        extra = [
            site / "5.22.0",
            site / "5.22.0" / "x86_64-linux-thread-multi",
            site / "5.22.2",
            site / "5.22.2" / "x86_64-linux-thread-multi",
        ]
        extra_paths = [str(p) for p in extra if p.exists()]
        if extra_paths:
            old = env.get("PERL5LIB", "")
            env["PERL5LIB"] = ":".join(extra_paths + ([old] if old else []))
    return env


def run_one_mlst(fasta_path: str, mlst_cmd: str, timeout_s: Optional[int]) -> tuple[str, int, str, str]:
    """Run mlst on one genome FASTA.

    Returns: (fasta_path, returncode, stdout, stderr)
    """
    env = build_env()
    try:
        proc = subprocess.run(
            [mlst_cmd, fasta_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            env=env,
            timeout=timeout_s,
        )
        return (fasta_path, int(proc.returncode), proc.stdout or "", proc.stderr or "")
    except FileNotFoundError:
        return (fasta_path, 127, "", f"mlst not found: {mlst_cmd}")
    except subprocess.TimeoutExpired as e:
        stderr = (e.stderr or "") if isinstance(e.stderr, str) else ""
        return (fasta_path, 124, e.stdout or "", (stderr + "\nTIMEOUT").strip())


def main() -> None:
    ap = argparse.ArgumentParser(description="Run mlst on genomes listed in genome_paths.tsv (resumable)")
    ap.add_argument("--genome-paths", required=True, help="TSV from step2_02_index_genomes.py")
    ap.add_argument("--out", required=True, help="Output TSV (mlst results; appended)")
    ap.add_argument("--mlst-cmd", default="mlst", help="mlst command")
    ap.add_argument("--resume", action="store_true", help="skip FASTAs already present in --out")
    ap.add_argument("--jobs", type=int, default=1, help="number of worker processes (default: 1)")
    ap.add_argument(
        "--timeout",
        type=int,
        default=0,
        help="per-genome timeout seconds (0 = no timeout)",
    )
    ap.add_argument("--progress-every", type=int, default=50, help="print progress every N genomes")
    ap.add_argument("--stderr-log", default=None, help="optional file to append mlst stderr")
    args = ap.parse_args()

    genome_paths = Path(args.genome_paths)
    if not genome_paths.exists() or genome_paths.stat().st_size == 0:
        raise SystemExit(f"ERROR: genome paths TSV missing or empty: {genome_paths}")

    paths = read_genome_paths(genome_paths)
    if not paths:
        raise SystemExit("ERROR: no genomes with status=ok found in genome paths TSV")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    done = read_completed(out_path) if args.resume else set()
    todo = [p for p in paths if p not in done]
    if args.resume:
        print(f"[MLST] completed: {len(done)}", flush=True)
        print(f"[MLST] remaining: {len(todo)}", flush=True)

    if not todo:
        print("[MLST] nothing to do", flush=True)
        return

    jobs = int(args.jobs)
    if jobs < 1:
        raise SystemExit("ERROR: --jobs must be >= 1")
    jobs = min(jobs, len(todo))

    timeout_s: Optional[int] = int(args.timeout) if int(args.timeout) > 0 else None

    log_fh = None
    if args.stderr_log:
        log_path = Path(args.stderr_log)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_fh = log_path.open("a", encoding="utf-8")

    try:
        with out_path.open("a", encoding="utf-8") as out:
            processed = 0

            if jobs == 1:
                for fasta_path in todo:
                    fasta_path, rc, stdout, stderr = run_one_mlst(fasta_path, args.mlst_cmd, timeout_s)
                    if log_fh is not None and stderr:
                        log_fh.write(f"# {fasta_path}\n")
                        log_fh.write(stderr)
                        if not stderr.endswith("\n"):
                            log_fh.write("\n")
                        log_fh.flush()

                    if rc != 0:
                        msg = (stderr or "").strip().replace("\n", " ")
                        out.write(f"# ERROR\t{fasta_path}\texit={rc}\t{msg}\n")
                    else:
                        for line in (stdout or "").splitlines():
                            if line.strip():
                                out.write(line.rstrip("\n") + "\n")

                    out.flush()
                    processed += 1
                    if args.progress_every > 0 and processed % int(args.progress_every) == 0:
                        print(f"[MLST] processed {processed}/{len(todo)}", flush=True)
            else:
                print(f"[MLST] jobs: {jobs}", flush=True)
                with cf.ProcessPoolExecutor(max_workers=jobs) as ex:
                    futures = [ex.submit(run_one_mlst, p, args.mlst_cmd, timeout_s) for p in todo]
                    for fut in cf.as_completed(futures):
                        fasta_path, rc, stdout, stderr = fut.result()

                        if log_fh is not None and stderr:
                            log_fh.write(f"# {fasta_path}\n")
                            log_fh.write(stderr)
                            if not stderr.endswith("\n"):
                                log_fh.write("\n")
                            log_fh.flush()

                        if rc != 0:
                            msg = (stderr or "").strip().replace("\n", " ")
                            out.write(f"# ERROR\t{fasta_path}\texit={rc}\t{msg}\n")
                        else:
                            for line in (stdout or "").splitlines():
                                if line.strip():
                                    out.write(line.rstrip("\n") + "\n")

                        out.flush()
                        processed += 1
                        if args.progress_every > 0 and processed % int(args.progress_every) == 0:
                            print(f"[MLST] processed {processed}/{len(todo)}", flush=True)
    finally:
        if log_fh is not None:
            log_fh.close()

    print(f"Wrote/updated: {out_path}", flush=True)


if __name__ == "__main__":
    main()
