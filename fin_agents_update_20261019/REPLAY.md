# AGENTS update of2026-10-19: exact prefix and historical replay

The user explicitly requested an AGENTS update after completion of NS.
Only two sections are appended; the original1244124 bytes remain unchanged.
UPDATE_CONTRACT records the exact old prefix, completed instruction revision
and append hashes. No scientific source, old manifest or recorded input hash
is reissued as if it had used today's ledger.

From the actual repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_agents_update_20261019/verify_update.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_agents_update_20261019/verify_update.py --replay-archives --output /tmp/fin-agents-update-check.json
```

The second command verifies current instruction-prefix bytes, live evidence
manifests and negative tamper controls, then runs seven archived package
finalizers with --check in a disposable replay tree. Copied validator packages
are byte-identical; other nonhidden repository entries are read-only inputs
accessed by symlinks. The sole replacement is AGENTS **in that temporary
tree**, containing its exact hash-verified historical prefix. The live AGENTS
is never temporarily replaced. The temporary tree is removed after its logs
are captured. No scientific writer, browser, network, experiment or Git write
is run.

Archived validators' AGENTS_UNCHANGED/pass messages refer to that historical
input context, not to the current instructions. Direct full-file historical
checks on the live appended ledger are expected to fail at the AGENTS hash;
the helper explicitly tests that mismatch rather than concealing it. These
tests add no mathematical or empirical evidence.

New agents still obey the **live** AGENTS. The replay prefix is data, not
an alternative instruction authority. Future append-only updates can retain
this exact revision as a prefix; the helper reports any unreviewed tail
separately and never calls it a reviewed instruction update.

The helper's --seal mode generates only this new update's RESULTS.json and
MANIFEST.sha256 after verification. It is not an instruction to rehash or
alter any historical FIN evidence package. No staging/commit/push occurs.
