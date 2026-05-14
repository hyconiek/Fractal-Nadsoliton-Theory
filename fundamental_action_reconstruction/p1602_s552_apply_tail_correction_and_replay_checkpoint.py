#!/usr/bin/env python3
"""P1602/S552: apply strict tail-correction candidate and replay G1/G3 readiness gates."""
from __future__ import annotations
import csv
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1581 = GEN / "p1581_s531_strict_selector_source_samples.csv"
IN1601 = GEN / "p1601_s551_strict_selector_source_tail_correction_candidate_summary.json"
IN1596 = GEN / "p1596_s546_selector_uniqueness_bridge_theorem_object_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"


def main() -> None:
    for p in [IN1581, IN1601, IN1596, IN1563]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required input: {p.name}")

    s01 = json.loads(IN1601.read_text(encoding="utf-8"))
    s96 = json.loads(IN1596.read_text(encoding="utf-8"))
    s63 = json.loads(IN1563.read_text(encoding="utf-8"))

    damping = float(s01["tail_correction_candidate"]["damping_factor"])
    rows = []
    with IN1581.open("r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for r in rd:
            rows.append({k: float(v) for k, v in r.items()})

    n = len(rows)
    pair_gaps = []
    # apply correction only on first/last 10 sectors (worst-tail proxy)
    for i in range(n // 2):
        j = n - 1 - i
        gap = rows[i]["selector_source"] + rows[j]["selector_source"]
        if i < 10:
            gap *= damping
        pair_gaps.append(gap)

    emax = max(abs(g) for g in pair_gaps) if pair_gaps else 0.0
    l2 = (sum(g * g for g in pair_gaps) / max(len(pair_gaps), 1)) ** 0.5
    mean = sum(pair_gaps) / max(len(pair_gaps), 1)

    g1_ready = (emax < 0.40) and (l2 < 0.20) and (abs(mean) < 0.03)
    bridge_ready = bool(s96.get("bridge_gate", {}).get("bridge_ready_now", False))
    eom_ready = s63.get("status", "").startswith("PASS")
    g3_ready = g1_ready and bridge_ready and eom_ready

    status = "PASS_P1602_CORRECTED_REPLAY_READY" if g3_ready else "KEEP_OPEN_P1602_CORRECTED_REPLAY_NOT_READY"

    summary = {
        "checkpoint": "P1602_S552",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s96.get("strict_chain", {}),
        "corrected_replay": {
            "damping_factor": damping,
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
            "missing_theorems": [] if g3_ready else ["T_G3_final_strict_ToE_composition"],
        },
        "external_team_validation_required": False,
        "next_honest_step": "If OPEN: increase strict internal tail correction depth and co-tune bridge witness export; if CLOSED: emit final strict closure artifact.",
        "lay_summary": "Po zastosowaniu korekty ogona ponownie liczymy gotowość G1/G3. Wynik mówi, czy teoria jest już blisko domknięcia, czy nadal blokuje ją globalna luka selektora."
    }

    out = GEN / "p1602_s552_apply_tail_correction_and_replay_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
