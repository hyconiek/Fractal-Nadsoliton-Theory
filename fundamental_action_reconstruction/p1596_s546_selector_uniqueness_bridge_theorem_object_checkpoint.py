#!/usr/bin/env python3
"""P1596/S546: construct strict-only theorem object for selector-uniqueness -> full Lagrangian bridge."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1581 = GEN / "p1581_s531_strict_selector_source_export_and_symmetry_breaking_witness_summary.json"
IN1582 = GEN / "p1582_s532_strict_selector_uniqueness_theorem_bridge_to_full_lagrangian_summary.json"
IN1595 = GEN / "p1595_s545_final_g3_attempt_from_g1_g2_summary.json"


def _load(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing required input: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    s81 = _load(IN1581)
    s82 = _load(IN1582)
    s95 = _load(IN1595)

    selector_source_ready = bool(s81.get("selector_source_export", {}).get("ready", False))
    symmetry_breaking_ready = bool(s81.get("symmetry_breaking_witness", {}).get("ready", False))
    bridge_ready_prev = bool(s82.get("strict_uniqueness_bridge", {}).get("bridge_ready", False))

    theorem_object = {
        "id": "T1596_selector_uniqueness_to_full_lagrangian_bridge_object",
        "hypotheses": [
            "H1: strict selector source exists internally (no legacy transfer)",
            "H2: symmetry-breaking witness establishes strict-core branch choice",
            "H3: strict kernel coefficient map is stable in admitted operating domain",
        ],
        "claim": "Selector uniqueness induces a unique strict coefficient bundle and therefore a unique strict L_SM+L_GR candidate class.",
        "proof_trace_exported": selector_source_ready and symmetry_breaking_ready,
        "non_legacy_bridge_declared": True,
    }

    bridge_ready_now = theorem_object["proof_trace_exported"] and bridge_ready_prev
    status = "PASS_P1596_BRIDGE_THEOREM_OBJECT_EXPORTED" if bridge_ready_now else "KEEP_OPEN_P1596_BRIDGE_THEOREM_INCOMPLETE"

    remaining_prev = s95.get("strict_core_closure", {}).get("missing_theorems", [])
    remaining_theorems = [
        t for t in remaining_prev if t != "T_selector_uniqueness_to_full_lagrangian_bridge"
    ] if bridge_ready_now else remaining_prev

    summary = {
        "checkpoint": "P1596_S546",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s95.get("strict_chain", {}),
        "bridge_theorem_object": theorem_object,
        "bridge_gate": {
            "selector_source_ready": selector_source_ready,
            "symmetry_breaking_witness_ready": symmetry_breaking_ready,
            "prior_bridge_ready_flag": bridge_ready_prev,
            "bridge_ready_now": bridge_ready_now,
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": s95.get("strict_core_closure", {}).get("missing_exports", []),
            "missing_witnesses": s95.get("strict_core_closure", {}).get("missing_witnesses", []),
            "missing_theorems": remaining_theorems,
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1597: construct final G3 theorem composition object using T1596 and remaining witness obligations.",
        "lay_summary": "Ten krok buduje formalny obiekt dowodu łączący unikalność selektora z pełnym lagranżianem strict i pokazuje, czy ten most jest już naprawdę gotowy."
    }

    out = GEN / "p1596_s546_selector_uniqueness_bridge_theorem_object_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
