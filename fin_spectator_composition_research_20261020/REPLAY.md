# Exact replay

From repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_spectator_composition_research_20261020/verify_research.py --check --output /tmp/fin-spectator-composition-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The 35 bounded controls use only stdlib Fraction arithmetic. A nine-state
source is independently generated from its scaling flow, not by importing
producer helpers. Actual compression, matrix products, killed-fibre linear
systems and a 27-state uniform deletion are checked. Count deletion uses
the exact hypergeometric kernel. No large count lattice or floating root
is substituted for proof.

General product/autonomy theorems are analytical, not mechanically proved
by those controls. Authoring --record writes only this fresh directory's
results, duplicate search, input list and manifest. Later --check writes
only the requested /tmp result. Old scientific inputs and AGENTS remain
hash-bound and unchanged; do not rewrite historical hashes after future
instruction updates. No network or metadata writer is required for replay.

The --check documentation pass excludes fenced and indented mathematical
code before validating local Markdown links; a formula such as a matrix
factor followed by (P-I) is not treated as a missing file. This corrects a
validation false positive, not a scientific proof or model change.
