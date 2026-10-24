#!/bin/bash
#SBATCH -J qcchem
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 00:60:00
#SBATCH -o slurm-%j.out
# usage: sbatch slurm/casci_job.sh configs/<name>.yaml   (run from ~/qc-chem on cluster)
set -euo pipefail
source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate qchem
export OMP_NUM_THREADS=8 PYSCF_TMPDIR=/data/$USER/qc-chem/tmp PYSCF_MAX_MEMORY=90000
mkdir -p $PYSCF_TMPDIR
cd $HOME/qc-chem
export PYTHONPATH=$HOME/qc-chem
echo "host=$(hostname) start=$(date -Is) cfg=$1"
python scripts/run_system.py "$1" ${2:-}
echo "end=$(date -Is)"
