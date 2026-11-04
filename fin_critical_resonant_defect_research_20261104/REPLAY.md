# Replay

From the repository root at the bound commit, place this directory beside
the existing source directories, then run:

```sh
python3 fin_critical_resonant_defect_research_20261104/verify_research.py
```

Requires Python 3 standard library only. It writes verification.json in
this NEW directory, never changes inputs or historical manifests, and
exits nonzero on failure. Input hashes bind the actual live AGENTS read
for this study; no historical-prefix substitution is performed.

Controls check finite algebra and input preservation, NOT analytical
estimates, full-PDE time evolution, or independent review. The report is
a stopped partial result with the precise remaining OPEN-CR inequality.

The archive contains this new directory only. Existing repository files
are dependencies, not duplicate deliverables. Subsequent user instructions
authorize repository publication; see REPORT_EN.md.
