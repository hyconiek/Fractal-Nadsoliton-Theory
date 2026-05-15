#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1744 = GEN / "p1744_s694_strict_full_chain_milestone_register_kernel_to_qg_checkpoint.json"
OUT = GEN / "p1745_s695_strict_kernel_to_eom_qg_closure_no_skip_enforcement_checkpoint.json"


def main() -> None:
    p1744 = json.loads(IN1744.read_text(encoding="utf-8"))
    milestones = p1744.get("milestone_register", [])

    milestone_status = {m.get("id", "?"): m.get("status", "UNKNOWN") for m in milestones}

    no_skip_rules = {
        "R1": "M3 must be published before M4 closure claim",
        "R2": "M4 must be published before M5 theorem gate entry",
        "R3": "No strict-core closure claim unless QG theorem witnesses exported",
    }

    payload = {
        "checkpoint": "P1745_S695",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full nonskeleton L_total -> EOM -> milestone register -> no-skip enforcement",
        "full_lagrangian_density_nonskeleton_instantiated": p1744.get(
            "full_lagrangian_density_nonskeleton_instantiated", {}
        ),
        "milestone_status_anchor": milestone_status,
        "no_skip_enforcement_rules": no_skip_rules,
        "current_gate_state": {
            "M3": milestone_status.get("M3", "UNKNOWN"),
            "M4": milestone_status.get("M4", "UNKNOWN"),
            "M5": milestone_status.get("M5", "UNKNOWN"),
            "strict_core_closure_allowed": False,
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Zrealizować M3 i M4 bez przeskoków; dopiero potem otworzyć M5 i dodawać theorem witnesses dla renormalizacji, unitarności i background-independence.",
        "lay_summary": "To blokada przed przeskakiwaniem etapów. Najpierw dwa brakujące testy, potem dopiero wielkie dowody QG — inaczej nie wolno ogłaszać domknięcia.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
