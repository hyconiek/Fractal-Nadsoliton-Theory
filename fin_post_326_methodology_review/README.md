# Post-326 methodology reconciliation

[REVIEW_20261014.md](REVIEW_20261014.md) records the scoped methodology findings
and corrections now reflected in AGENTS.md.

Run from the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_post_326_methodology_review/replay_checks.py
```

The checker reads existing artifacts and prints JSON; it does not modify source
campaigns or claim a physical hardware replay. verification.json records the
review-time output.
