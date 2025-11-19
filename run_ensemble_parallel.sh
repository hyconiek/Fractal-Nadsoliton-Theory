#!/usr/bin/env bash
# run_ensemble_parallel.sh
# Uruchamia ensemble runs równolegle (GNU parallel)
# Usage: ./run_ensemble_parallel.sh

PYTHON=python3
SCRIPT=110_FIX_SELFCONSISTENCY.py
# Lista N do sprawdzenia
N_LIST=(12 16 20 24 28 32)

export PYTHONPATH="$(pwd)"

# Funkcja uruchamiająca pojedynczy run (używana przez parallel)
run_one() {
  N=$1
  OUTDIR="ensemble_results_N${N}"
  mkdir -p "$OUTDIR"
  # Uruchom w trybie nie-dry-run, jedna iteracja krótsza dla szybkości -- można zmienić
  $PYTHON $SCRIPT --anchor 246.0 > "$OUTDIR/run.log" 2>&1
}

export -f run_one

# Uruchomienie równoległe: zmień -j do liczby równoległych zadań
parallel -j 6 run_one ::: "${N_LIST[@]}"
