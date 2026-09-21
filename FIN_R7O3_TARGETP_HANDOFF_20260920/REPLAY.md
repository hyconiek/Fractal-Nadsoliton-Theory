# Replay instructions

Run from the package root.

Fast final gate (uses completed clean replay ledger):

```bash
python verify.py
python src/global_join_checker.py
python src/hostile_mutations.py
```

Full mathematical replay from certificates (costly, resumable):

```bash
rm -rf checkpoints/clean_math_shards
mkdir -p checkpoints/clean_math_shards
# Run the declared ranges in bounded shards, e.g.
python src/replay_clean_math_multi.py --start 0 --end 200 --workers 8 \
  --output checkpoints/clean_math_shards/shard_00000_00199.json
# Continue ranges until 12425, then:
python src/aggregate_clean_replay.py
python verify.py
```

`certify_fixed` uses no eigensolver or optimizer. It consumes saved `B,c` and outward-rational proof inputs.
