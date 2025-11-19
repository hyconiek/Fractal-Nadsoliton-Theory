#!/bin/bash
#SBATCH --job-name=toe_ensemble
#SBATCH --output=toe_ensemble_%j.out
#SBATCH --error=toe_ensemble_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

module load anaconda/3
source activate your_env_name  # dopasuj do środowiska

PYTHON=python3
SCRIPT=/home/krzysiek/Pobrane/TOE/edison/110_FIX_SELFCONSISTENCY.py

# Parametry (można podać przez sbatch --export=...)
N=${N:-24}
ANCHOR=${ANCHOR:-246.0}

# Uruchomienie
$PYTHON $SCRIPT --anchor $ANCHOR > "ensemble_N${N}.log" 2>&1
