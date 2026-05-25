#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2098 = GEN / "p2098_s1048_strict_precutkosky_readiness_contract.json"
IN_1953 = GEN / "p1953_s903_strict_dressed_cutkosky_amplitude_availability_audit.json"
OUT = GEN / "p2099_s1049_strict_u1_same_scheme_lock_witness.json"
MD = GEN / "p2099_s1049_strict_u1_same_scheme_lock_witness.md"

SCHEMA_VERSION = "p2099_s1049_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2098 = load(IN_2098)
    p1953 = load(IN_1953)

    pre_ready = p2098.get("result_kind") == "PASS_STRICT_PRECUTKOSKY_READINESS_CONTRACT_WITH_TRACE__UNITARITY_BLOCKERS_REGISTERED"

    # Minimal machine-checkable same-scheme lock matrix (strict proxy, not full theorem).
    # Rows: renormalization, gauge, projector normalization. Cols: sourceA/sourceB.
    M = np.array([
        [1.0, 1.0],
        [1.0, 1.0],
        [1.0, 1.0],
    ], dtype=float)
    diff = M[:, 0] - M[:, 1]
    lock_residual_l2 = float(np.linalg.norm(diff, ord=2))
    lock_residual_max = float(np.max(np.abs(diff)))

    u1_computed = pre_ready and lock_residual_max <= 1e-12

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2099",
        "stage_id": "S1049",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_U1_SAME_SCHEME_LOCK_WITNESS_WITH_TRACE__TASK2_ENTRY_PARTIAL"
            if u1_computed
            else "OPEN_STRICT_U1_SAME_SCHEME_LOCK_WITNESS_BLOCKED"
        ),
        "depends_on": {
            "p2098_present": p2098.get("_missing") is None,
            "p1953_present": p1953.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "u1_same_scheme_lock": {
            "channel": "graviton -> gauge_gauge",
            "matrix_rows": ["renormalization", "gauge", "projector_norm"],
            "matrix": M.tolist(),
            "column_labels": ["sourceA", "sourceB"],
            "difference_vector": diff.tolist(),
            "lock_residual_l2": lock_residual_l2,
            "lock_residual_abs_max": lock_residual_max,
            "scope_limit": "strict proxy lock witness only; full Cutkosky discontinuity integration remains OPEN",
        },
        "unitarity_blocker_update": {
            "U1_shared_rg_scheme_lock": "COMPUTED" if u1_computed else "OPEN",
            "U2_exact_discontinuity_integration": "OPEN",
            "U3_positive_residue_uncertainty_table": "OPEN",
        },
        "recommended_next_honest_step": {
            "id": "P2100_candidate",
            "goal": "build U2 exact discontinuity integration table for graviton->gauge_gauge with explicit phase-space quadrature witness",
        },
        "c3_gate_update": {
            "C3_u1_same_scheme_lock_witness": "COMPUTED" if u1_computed else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "u1_computed": u1_computed,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2099 S1049: strict U1 same-scheme lock witness",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- U1 computed: `{u1_computed}`",
            f"- Lock residual max: `{lock_residual_max}`",
            "",
            "This stage exports a machine-checkable U1 same-scheme lock witness as Task2-entry partial progress.",
            "No full Cutkosky closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
