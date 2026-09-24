# Audit and provenance gaps

1. **Wave 0–2 raw execution packages are absent from the active runtime.** Their terminal dispositions were consumed by later dependency-closed work, but their individual `REPORT.md`, theorem/counterexample note, `results.json`, replay and manifest cannot be independently checked from this ZIP.
2. **Wave 3 is summary-only in the active runtime.** The summary records the five outcomes, but per-task bytes/manifests are absent here.
3. **Wave 4 raw archive contains provisional stubs for DYN-002 and GATE-001.** The normalized handoff deliberately excludes those stubs; the authoritative versions are the Wave 5 packages. The untouched Wave 4 ZIP is preserved under `03_ORIGINAL_PACKAGES/` for provenance.
4. No external physical experiment or data collection was performed. No GitHub writes were performed during this continuation.

These gaps do not reopen the 42-node execution ledger; they limit the strength of a byte-level global reproducibility claim. `SYN-001` therefore correctly remains fail-closed rather than reporting a blanket PASS.
