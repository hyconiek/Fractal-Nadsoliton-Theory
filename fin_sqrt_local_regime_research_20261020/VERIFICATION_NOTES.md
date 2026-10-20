# Reproducibility repair retained

Date: 2026-10-20. The first --record run passed all finite controls, input
checks and protected manifests. A subsequent --check rejected only the
DUPLICATE_SEARCH.json comparison: ripgrep's parallel output interleaved
the same matched lines in a different file order.

The diagnostic compared sorted match multisets and searched-file lists;
both were identical. The helper now canonicalizes the ordering of search
lines before recording/comparison. No scientific identity, rate, threshold,
proof, predecessor file or input hash was changed by this repair. The
original failed assertion is recorded here rather than counted as a passed
full-package replay. The final read-only checker must pass after this fix.

This is an output-order reproducibility bug, not a mathematical
counterexample or a reason to modify historical evidence. All new analytical
results remain self-audited conditional claims, not independent review.
