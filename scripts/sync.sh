#!/usr/bin/env bash
# Sync between the local git copy and the cluster mirror (no rsync on Windows; tar over ssh).
#   scripts/sync.sh push   -> code (qcchem/ configs/ scripts/ slurm/ tests/ *.md environment.yml) to cluster:~/qc-chem
#   scripts/sync.sh pull   -> cluster:~/qc-chem/results -> local results/  (merge, never delete locally)
#   scripts/sync.sh check  -> print remote tree size and scratch usage
set -euo pipefail
HOST=cluster
REMOTE=~/qc-chem
SCRATCH=/data/$USER/qc-chem
cd "$(dirname "$0")/.."

case "${1:-}" in
  push)
    ssh -o BatchMode=yes "$HOST" "mkdir -p $REMOTE $SCRATCH"
    tar czf - --exclude='__pycache__' --exclude='.pytest_cache' \
        qcchem configs scripts slurm tests PLAN.md TASK_BRIEF.md progress.md environment.yml 2>/dev/null \
      | ssh -o BatchMode=yes "$HOST" "tar xzf - -C $REMOTE"
    echo "push ok -> $HOST:$REMOTE"
    ;;
  pull)
    mkdir -p results
    ssh -o BatchMode=yes "$HOST" "cd $REMOTE && [ -d results ] && tar czf - results || tar czf - -T /dev/null" \
      | tar xzf - -C .
    echo "pull ok <- $HOST:$REMOTE/results"
    ;;
  check)
    ssh -o BatchMode=yes "$HOST" "du -sh $REMOTE $SCRATCH 2>/dev/null; squeue -u $USER -h | wc -l | xargs echo queued_jobs:"
    ;;
  *)
    echo "usage: $0 {push|pull|check}" >&2; exit 2;;
esac
