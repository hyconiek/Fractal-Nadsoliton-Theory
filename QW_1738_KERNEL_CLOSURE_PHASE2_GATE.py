#!/usr/bin/env python3
"""
QW-1738: Phase-2 closure gate for Nadsoliton->Kernel path.

Aggregates 1730..1737 studies into one strict decision.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1738_kernel_closure_phase2_gate.json"
OUT_MD = ROOT / "RAPORT_QW1738_KERNEL_CLOSURE_PHASE2_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1730 = load("report_qw1730_nadsoliton_kernel_chrono_audit.json")
    d1731 = load("report_qw1731_nadsoliton_kernel_node_compatibility.json")
    d1732 = load("report_qw1732_fractal_path_to_hyperbolic_test.json")
    d1734 = load("report_qw1734_micro_beta_tors_derivation.json")
    d1735 = load("report_qw1735_omega_phi_from_lattice_dynamics.json")
    d1736 = load("report_qw1736_kernel_node_bayesian_model_selection.json")
    d1737 = load("report_qw1737_shared_kernel_flavor_gw_cross_constraint.json")

    comp: List[Dict[str, object]] = []

    # 1730
    risk = int(d1730.get("risk_points", 99))
    s1730 = max(0.0, 1.0 - risk / 10.0)
    comp.append(
        {
            "domain": "Chronology consistency (1730)",
            "score": s1730,
            "status": "PASS" if risk <= 3 else "FAIL",
            "note": d1730.get("verdict", ""),
        }
    )

    # 1731
    v1731 = d1731.get("verdict", "")
    s1731 = {"NODE_NARRATIVE_COMPATIBLE": 1.0, "NODE_NARRATIVE_PARTIALLY_INCONSISTENT": 0.4}.get(v1731, 0.0)
    comp.append(
        {
            "domain": "Node compatibility (1731)",
            "score": s1731,
            "status": "PASS" if s1731 >= 0.7 else "FAIL",
            "note": v1731,
        }
    )

    # 1732
    v1732 = d1732.get("verdict", "")
    s1732 = {
        "HYPERBOLIC_REDUCTION_ROBUST": 1.0,
        "HYPERBOLIC_REDUCTION_PLAUSIBLE_BUT_TUNED": 0.45,
    }.get(v1732, 0.1)
    comp.append(
        {
            "domain": "Hyperbolic denominator robustness (1732)",
            "score": s1732,
            "status": "PASS" if s1732 >= 0.7 else "FAIL",
            "note": v1732,
        }
    )

    # 1734
    flags1734 = d1734.get("pass_flags", {})
    s1734 = (
        (0.4 if flags1734.get("fit_quality") else 0.0)
        + (0.3 if flags1734.get("target_beta_match_pm10pct") else 0.0)
        + (0.3 if flags1734.get("oos_stability") else 0.0)
    )
    comp.append(
        {
            "domain": "Microscopic beta derivation (1734)",
            "score": s1734,
            "status": "PASS" if s1734 >= 0.7 else "FAIL",
            "note": d1734.get("verdict", ""),
        }
    )

    # 1735
    flags1735 = d1735.get("pass_flags", {})
    s1735 = (0.6 if flags1735.get("stable_pair") else 0.0) + (0.4 if flags1735.get("near_reference") else 0.0)
    comp.append(
        {
            "domain": "Omega/Phi derivation from dynamics (1735)",
            "score": s1735,
            "status": "PASS" if s1735 >= 0.7 else "FAIL",
            "note": d1735.get("verdict", ""),
        }
    )

    # 1736
    pbest = float(max(d1736.get("posterior_probabilities", {}).values()))
    compat = bool(d1736.get("compatibility_with_characteristic_priors", False))
    s1736 = 0.6 * pbest + (0.4 if compat else 0.0)
    comp.append(
        {
            "domain": "Bayesian node model support (1736)",
            "score": s1736,
            "status": "PASS" if s1736 >= 0.7 else "FAIL",
            "note": d1736.get("verdict", ""),
        }
    )

    # 1737
    gates = d1737.get("empirical_hard_gates", {})
    frac_strict = float(d1737.get("sampling", {}).get("shared_fraction_strict", 0.0))
    ident = float(d1737.get("identifiability", {}).get("relative_entropy_reduction", 0.0))
    s1737 = (
        (0.35 if gates.get("flavor_gate_from_1720_1723") else 0.0)
        + (0.35 if gates.get("gw_gate_from_1725_1726") else 0.0)
        + min(0.2, 4.0 * frac_strict)
        + min(0.1, ident)
    )
    comp.append(
        {
            "domain": "Shared flavor+GW region (1737)",
            "score": s1737,
            "status": "PASS" if s1737 >= 0.7 else "FAIL",
            "note": d1737.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in comp) / len(comp))
    hard_gate = all(c["status"] == "PASS" for c in comp)

    if hard_gate and global_score >= 0.8:
        readiness = "PHASE2_KERNEL_PATH_CLOSED"
    elif global_score >= 0.6:
        readiness = "PHASE2_KERNEL_PATH_PARTIAL"
    else:
        readiness = "PHASE2_KERNEL_PATH_OPEN"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "components": comp,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1738: KERNEL CLOSURE PHASE2 GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{output['hard_gate']}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Components",
    ]
    for c in comp:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1738] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1738] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
