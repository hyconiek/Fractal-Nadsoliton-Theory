# R7P-001 baseline provenance

Created: 2026-09-13T21:06:03.608629+00:00

## Available mandatory inputs
- Master plan: present and hashed.
- Current accessible AGENTS snapshot: present and hashed, but the exact plan-referenced intake section could not be located in this snapshot.
- Imported handoff and manifest: present and hashed.
- Discord cutoff/source text: present.

## Missing mandatory audit package
The local `fin_handoff_audit/` directory and its required `REPORT.md`, `PROOF.md`, `research.py`, `test_research.py`, `verify.py`, `results.json`, and `verification.json` are not present in the supplied container inputs. R7P-002 is therefore BLOCKED_DEPENDENCY. No baseline expected result is regenerated from handoff seeds.

## Repository state
No Git repository is mounted at `/mnt/data`; commit and dirty-tree provenance cannot be established from this environment. Imported files are copied into `inputs/` and marked read-only by policy.

## Resource revision
The current runtime has ~5.8 GiB RAM and 5 CPUs, materially below the historical 16-GB host assumption. Any later class-M/L execution in this runtime must use smaller caps/checkpoints.
