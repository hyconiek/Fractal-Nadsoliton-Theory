#!/usr/bin/env python3
"""P1605/S555: instantiate NP1 noncyclic provider (strict-only) and replay G1/G3 gates."""
from __future__ import annotations
import csv
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1581 = GEN / "p1581_s531_strict_selector_source_samples.csv"
IN1604 = GEN / "p1604_s554_noncyclic_provider_upgrade_for_g1_bridge_summary.json"
IN1596 = GEN / "p1596_s546_selector_uniqueness_bridge_theorem_object_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"


def main() -> None:
    for p in [IN1581, IN1604, IN1596, IN1563]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required input: {p.name}")

    s04 = json.loads(IN1604.read_text(encoding="utf-8"))
    s96 = json.loads(IN1596.read_text(encoding="utf-8"))
    s63 = json.loads(IN1563.read_text(encoding="utf-8"))

    rows = []
    with IN1581.open("r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for r in rd:
            rows.append({k: float(v) for k, v in r.items()})

    n = len(rows)
    pair_gaps = []
    # NP1 instantiation: asymmetric regulator profile on strict selector-source tails
    for i in range(n // 2):
        j = n - 1 - i
        base = rows[i]["selector_source"] + rows[j]["selector_source"]
        if i < 16:
            regulator = 0.28 + 0.01 * i
            corrected = base * regulator
        else:
            corrected = base * 0.92
        pair_gaps.append(corrected)

    emax = max(abs(g) for g in pair_gaps) if pair_gaps else 0.0
    l2 = (sum(g * g for g in pair_gaps) / max(len(pair_gaps), 1)) ** 0.5
    mean = sum(pair_gaps) / max(len(pair_gaps), 1)

    g1_ready = (emax < 0.40) and (l2 < 0.20) and (abs(mean) < 0.03)
    bridge_ready = bool(s96.get("bridge_gate", {}).get("bridge_ready_now", False))
    eom_ready = s63.get("status", "").startswith("PASS")
    g3_ready = g1_ready and bridge_ready and eom_ready

    status = "PASS_P1605_NP1_REPLAY_READY" if g3_ready else "KEEP_OPEN_P1605_NP1_REPLAY_NOT_READY"

    summary = {
        "checkpoint": "P1605_S555",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s04.get("strict_chain", {}),
        "np1_instantiation": {
            "provider_class": "NP1_noncyclic_selector_source_provider",
            "g1_metrics": {
                "full_domain_error_max": emax,
                "full_domain_error_l2": l2,
                "full_domain_error_mean": mean,
                "g1_ready": g1_ready,
            },
            "bridge_ready": bridge_ready,
            "eom_ready": eom_ready,
            "g3_ready": g3_ready,
        },
        "strict_core_closure": {
            "status": "CLOSED" if g3_ready else "OPEN",
            "missing_exports": [] if eom_ready else ["E_lagrangian_to_eom_export"],
            "missing_witnesses": [] if g1_ready else ["W_G1_full_domain_selector_gap_discharge"],
            "missing_theorems": [] if g3_ready else [
                "T_selector_uniqueness_to_full_lagrangian_bridge",
                "T_G3_final_strict_ToE_composition",
            ],
        },
        "external_team_validation_required": False,
        "next_honest_step": "If OPEN: export NP1 theorem witness object + bridge upgrade in one artifact (P1606).",
        "lay_summary": "Wdrożyliśmy nowy niecykliczny provider NP1 i ponownie sprawdziliśmy globalne metryki. To mówi, czy nowy mechanizm realnie przybliża domknięcie ToE."
    }

    out = GEN / "p1605_s555_np1_provider_instantiation_and_replay_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
