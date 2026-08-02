#!/usr/bin/env bash
set -eu

cd "$(dirname "$0")"
nohup python3 -u run_fin_p475_p485_p486_p487_unlimited.py \
  > FIN_Unlimited_Research_Pipeline.log 2>&1 &
pipeline_pid=$!
printf '%s\n' "$pipeline_pid" > FIN_Unlimited_Research_Pipeline.pid
printf 'FIN unlimited research pipeline started with PID %s\n' "$pipeline_pid"
