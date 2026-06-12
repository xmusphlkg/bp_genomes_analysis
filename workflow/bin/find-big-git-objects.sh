#!/usr/bin/env bash
# List the largest git objects (blobs) and their paths
# Usage: run from repo root: ./scripts/find-big-git-objects.sh 50

set -euo pipefail
N=${1:-50}

echo "Finding top $N largest git objects (blobs) in repository..."

# Produce list of object SHA, size, path
git rev-list --objects --all \
  | git cat-file --batch-check='%(objectname) %(objecttype) %(objectsize) %(rest)' \
  | awk '$2=="blob" {print $1, $3, substr($0, index($0,$4))}' \
  | sort -k2 -n \
  | tail -n $N \
  | awk '{printf("%s\t%s\t%s\n", $1, $2, substr($0, index($0,$3)))}'

echo "Done. Use the paths shown to decide which files to remove or move to Git LFS." 
