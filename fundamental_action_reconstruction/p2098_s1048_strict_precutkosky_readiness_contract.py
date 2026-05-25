#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2096 = GEN / "p2096_s1046_strict_b1_renormalization_closure_contract.json"
IN_2097 = GEN / "p2097_s1047_strict_b1_quotient_closure_stability_minigrid.json"
IN_1953 = GEN / "p1953_s903_strict_dressed_cutkosky_amplitude_availability_audit.json"
OUT = GEN / "p2098_s1048_strict_precutkosky_readiness_contract.json"
MD = GEN / "p2098_s1048_strict_precutkosky_readiness_contract.md"

SCHEMA_VERSION = "p2098_s1048_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2096 = load(IN_2096)
    p2097 = load(IN_2097)
    p1953 = load(IN_1953)

    pre_ready = (
        p2096.get("result_kind") == "PASS_STRICT_B1_RENORMALIZATION_CLOSURE_CONTRACT_WITH_TRACE__QUOTIENT_SCOPE_ONLY"
        and p2097.get("result_kind") == "PASS_STRICT_B1_QUOTIENT_CLOSURE_STABILITY_MINIGRID_WITH_TRACE"
    )

    reqs = [
        {
            "id": "U1",
            "name": "shared_rg_scheme_lock_for_graviton_to_gauge_gauge",
            "status": "OPEN",
            "evidence": "p1953 same-scheme requirement remains open",
        },
        {
            "id": "U2",
            "name": "exact_phase_space_discontinuity_integration_for_optical_theorem",
            "status": "OPEN",
            "evidence": "no exported full discontinuity integral table in current chain",
        },
        {
            "id": "U3",
            "name": "strict_positive_residue_window_with_uncertainty_bounds",
            "status": "OPEN",
            "evidence": "no machine-certified positivity uncertainty table for target channel",
        },
    ]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2098",
        "stage_id": "S1048",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_PRECUTKOSKY_READINESS_CONTRACT_WITH_TRACE__UNITARITY_BLOCKERS_REGISTERED"
            if pre_ready
            else "OPEN_STRICT_PRECUTKOSKY_READINESS_CONTRACT_BLOCKED"
        ),
        "depends_on": {
            "p2096_present": p2096.get("_missing") is None,
            "p2097_present": p2097.get("_missing") is None,
            "p1953_present": p1953.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "precutkosky_contract": {
            "channel": "graviton -> gauge_gauge",
            "renormalization_input_contract": {
                "source": "P2096 quotient-scope closure contract",
                "stability_support": "P2097 mini-grid local stability",
                "scope_limit": "strict B1 quotient scope only",
            },
            "unitarity_blocker_register": reqs,
            "entry_rule": "Task2 execution may start only after at least one of U1/U2/U3 is converted to COMPUTED with machine-checkable witness.",
        },
        "recommended_next_honest_step": {
            "id": "P2099_candidate",
            "goal": "build U1 same-scheme lock witness for graviton->gauge_gauge and export strict RG-scheme compatibility matrix",
        },
        "c3_gate_update": {
            "C3_precutkosky_readiness_contract": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "all_unitarity_blockers_registered_open": all(r["status"] == "OPEN" for r in reqs),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2098 S1048: strict pre-Cutkosky readiness contract",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Registered unitarity blockers: `{len(reqs)}`",
            f"- Preconditions ready: `{pre_ready}`",
            "",
            "This stage prepares entry into Task 2 by exporting a strict pre-Cutkosky readiness contract with explicit blocker register.",
            "No Cutkosky closure or ToE closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
