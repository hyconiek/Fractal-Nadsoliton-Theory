#!/usr/bin/env python3
"""P1517 S4.67: rerun strict coupled theorem with lock and export locked witness set."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1516 = GEN / "p1516_s466_theorem_annex_consistency_lock_for_strict_f_to_lsm_lgr_summary.json"
SUMMARY = GEN / "p1517_s467_locked_witness_set_export_for_strict_coupled_theorem_summary.json"


def main() -> None:
    s1516 = json.loads(P1516.read_text(encoding="utf-8"))
    rows = s1516.get("probe_results", [])

    locked_witness_set = [r for r in rows if bool(r.get("theorem_eval_allowed", False))]
    blocked_set = [r for r in rows if not bool(r.get("theorem_eval_allowed", False))]

    summary = {
        "packet": "P1517",
        "status": "PASS_LOCKED_WITNESS_SET_EXPORTED",
        "scope": "STRICT_ONLY_F_NADSOLITON_TO_LSM_PLUS_LGR",
        "lock_rule": s1516.get("lock_rule"),
        "locked_witness_set": locked_witness_set,
        "blocked_scenarios": blocked_set,
        "counts": {
            "locked_witness_count": len(locked_witness_set),
            "blocked_count": len(blocked_set),
        },
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Implement P1518 strict theorem witness compression: derive the minimal locked witness basis and map it to F->LSM and F->LGR channel obligations.",
        "layman_explanation": "Po włączeniu bezpiecznika dostajemy czystą listę przypadków, którym wolno ufać. To porządkuje dowód: wiemy dokładnie, które scenariusze są legalne, a które odpadają.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1517] status={summary['status']} locked={len(locked_witness_set)} blocked={len(blocked_set)}")


if __name__ == "__main__":
    main()
