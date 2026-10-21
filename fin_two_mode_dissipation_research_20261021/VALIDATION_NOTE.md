# Validation timing note

Date: 2026-10-21.

After a final clarification in VIRIAL_SCOPE_AND_THRESHOLD_GATE.md, the
own metadata writer yielded a live process handle. A read-only check was
then started before that writer finished. It read the old own manifest
and correctly rejected the clarified document's changed hash. The writer
subsequently completed successfully; the check did not modify evidence.

This was an orchestration race, not a failed scientific inequality or a
predecessor hash change. Final metadata authorship and read-only replay
were run in separate completed steps. All source manifests and claims
remain unchanged; no scientific condition was weakened to fix the timing.
