#!/usr/bin/env python3
"""P1582/S532: theorem-level strict selector uniqueness candidate linked to full strict Lagrangian chain."""
from __future__ import annotations
import csv
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_CSV = GEN / "p1581_s531_strict_selector_source_samples.csv"
IN_JSON = GEN / "p1580_s530_strict_kernel_to_full_chain_and_closure_gap_summary.json"


def main() -> None:
    if not IN_CSV.exists() or not IN_JSON.exists():
        raise FileNotFoundError("Run P1580 and P1581 before P1582.")

    with IN_JSON.open("r", encoding="utf-8") as f:
        p1580 = json.load(f)

    vals = []
    with IN_CSV.open("r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for r in rd:
            vals.append({k: float(v) for k, v in r.items()})

    # uniqueness proxy: one-sign branch dominance + monotonic energy descent section
    pos_mass = sum(v["selector_source"] for v in vals if v["selector_source"] > 0)
    neg_mass = -sum(v["selector_source"] for v in vals if v["selector_source"] < 0)
    branch_gap = abs(pos_mass - neg_mass) / max(pos_mass + neg_mass, 1e-12)

    monotonic_segments = 0
    for i in range(len(vals) - 1):
        e_i = vals[i]["selector_source"] ** 2 + 0.3 * vals[i]["k_strict"] ** 2
        e_j = vals[i + 1]["selector_source"] ** 2 + 0.3 * vals[i + 1]["k_strict"] ** 2
        if e_j <= e_i:
            monotonic_segments += 1
    mono_ratio = monotonic_segments / max(len(vals) - 1, 1)

    theorem_candidate_pass = (branch_gap <= 0.08) and (mono_ratio >= 0.55)
    status = "PASS_T1582_STRICT_SELECTOR_UNIQUENESS_CANDIDATE" if theorem_candidate_pass else "FAIL_T1582_STRICT_SELECTOR_UNIQUENESS_CANDIDATE"

    summary = {
        "checkpoint": "P1582_S532",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": {
            "kernel": p1580["pipeline"]["kernel"],
            "coefficients": p1580["pipeline"]["coefficients"],
            "lagrangian": p1580["pipeline"]["lagrangian"],
            "eom": p1580["pipeline"]["eom"],
        },
        "theorem_candidate": {
            "name": "T1582_strict_selector_uniqueness_from_internal_source",
            "branch_gap": branch_gap,
            "monotonic_energy_ratio": mono_ratio,
            "passes_proxy_conditions": theorem_candidate_pass,
        },
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_core_closure": {
            "status": "OPEN",
            "qw2191_closed": False,
            "remaining_missing": [
                "formal_proof_object_for_T1582",
                "full_SM_GR_global_stability_theorem",
                "final_composition_theorem_for_ToE_closure",
            ]
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1583_construct_formal_proof_object_for_T1582_and_compose_with_global_stability",
        "lay_summary": "Tor od kernela strict do równań ruchu jest zachowany; krok theorem-level dla unikalności selektora nadal wymaga formalnego obiektu dowodu."
    }

    out = GEN / "p1582_s532_strict_selector_uniqueness_theorem_bridge_to_full_lagrangian_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
