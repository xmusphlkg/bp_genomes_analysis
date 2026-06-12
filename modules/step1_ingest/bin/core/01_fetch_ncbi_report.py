import argparse
import subprocess
import sys
from pathlib import Path

def run(cmd: list[str]) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed:\n{' '.join(cmd)}\n\nSTDERR:\n{p.stderr.strip()}")
    sys.stdout.write(p.stdout)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxon", required=True, help='e.g. "Bordetella pertussis"')
    ap.add_argument("--out", required=True, help="Output JSONL filename")
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Try --limit all (newer versions), else fallback to a huge number
    cmd_all = ["datasets", "summary", "genome", "taxon", args.taxon, "--as-json-lines", "--limit", "all"]
    cmd_big = ["datasets", "summary", "genome", "taxon", args.taxon, "--as-json-lines", "--limit", "1000000"]

    try:
        # write directly to file by redirecting stdout ourselves
        p = subprocess.run(cmd_all, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if p.returncode != 0 or not p.stdout.strip():
            # fallback
            p = subprocess.run(cmd_big, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if p.returncode != 0:
                raise RuntimeError(p.stderr.strip() or "unknown error")
        out.write_text(p.stdout, encoding="utf-8")
    except Exception as e:
        raise SystemExit(f"ERROR fetching report JSONL: {e}")

    print(f"Wrote: {out} (lines: {sum(1 for _ in out.open('r', encoding='utf-8'))})")

if __name__ == "__main__":
    main()
