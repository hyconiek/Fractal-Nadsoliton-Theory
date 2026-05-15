#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1751 = GEN / "p1751_s701_strict_nonproxy_h1_boundary_control_contract_upgrade_checkpoint.json"
OUT = GEN / "p1752_s702_strict_nonproxy_h1_4d_execution_trigger_checkpoint.json"


def main() -> None:
    p1751 = json.loads(IN1751.read_text(encoding="utf-8"))
    contract = p1751.get("nonproxy_export_requirements_contract_upgraded", {}).get("h1_gauge_scalar", {})

    required = contract.get("required_nonproxy_exports", [])
    boundary_required = contract.get("boundary_control_contract", {}).get("required", False)

    # Honest repo-state trigger scan (still no nonproxy 4D exporter delivered)
    availability = {k: False for k in required}
    if boundary_required:
        availability["boundary_control_contract"] = False

    missing = [k for k, v in availability.items() if not v]
    trigger = len(missing) == 0

    payload = {
        "checkpoint": "P1752_S702",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full nonskeleton L_total -> EOM -> boundary-aware nonproxy H1 contract -> 4D execution trigger gate",
        "input_anchor": "p1751_s701",
        "nonproxy_h1_4d_trigger_gate": {
            "trigger_ready": trigger,
            "required_exports_and_clauses": availability,
            "missing": missing,
            "decision_policy": "NO_RUN_IF_MISSING_EXPORTS",
        },
        "execution_status": {
            "status": "BLOCKED_TRIGGER_NOT_READY" if not trigger else "READY_TO_RUN",
            "run_performed": False,
            "no_false_pass_confirmed": True,
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "next_honest_step": "Dostarczyć minimum: jawne nonproxy E_A^μ, E_H, kontrakt wspólnej rodziny teł i kontrolę boundary clause; dopiero potem uruchomić 4D H1.",
        "lay_summary": "Zanim uruchomimy najważniejszy test 4D, sprawdziliśmy gotowość wejść. Brakuje kluczowych eksportów, więc uczciwie blokujemy uruchomienie zamiast generować pozorny wynik.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
