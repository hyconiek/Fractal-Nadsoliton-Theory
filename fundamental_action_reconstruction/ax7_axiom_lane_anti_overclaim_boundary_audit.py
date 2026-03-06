from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
out_dir = root / "generated"

source_packet_path = out_dir / "axiom_lane_closure_packet.json"
source_packet = json.loads(source_packet_path.read_text(encoding="ascii"))

certificate = {
    "lane": "axiom-augmented",
    "source_packet": "generated/axiom_lane_closure_packet.json",
    "source_packet_exists": source_packet_path.exists(),
    "strict_core_status": "outside_strict_core",
    "promotion_policy": {
        "strict_core_promotion_allowed": False,
        "theorem_level_promotion_allowed": False,
        "full_closure_promotion_allowed": False,
    },
    "forbidden_promotions": [
        "strict_core_selector_closure",
        "theorem_level_pass",
        "full_closure_pass",
        "qw_2191_discharge",
        "axiom_lane_equals_strict_core",
        "full_gauge_uniqueness_closure",
    ],
    "retained_axiom_lane_content": {
        "actual_theta": source_packet["actual_theta"],
        "actual_basis_pair": source_packet["actual_basis_pair"],
        "actual_orientation_slice": source_packet["actual_orientation_slice"],
        "bridge_instance": source_packet["bridge_instance"],
        "robustness": source_packet["robustness"],
        "compatibility": source_packet["compatibility"],
    },
    "strict_core_frontier_unchanged": source_packet["residual_frontier"],
    "best_supported_interpretation": "positive external selector lane only",
}

certificate_path = out_dir / "axiom_lane_boundary_certificate.json"
certificate_path.write_text(json.dumps(certificate, indent=2) + "\n", encoding="ascii")

summary = {
    "step": "AX7",
    "status": "AX7_EXECUTED_AXIOM_LANE_ANTI_OVERCLAIM_BOUNDARY_AUDIT_NO_FALSE_PASS",
    "goal": "Certify that AX1..AX6 remain a positive external lane only and may not be promoted into strict core or theorem-level/full-closure status.",
    "created_file": {
        "relative_path": "generated/axiom_lane_boundary_certificate.json",
        "exists_after_step": certificate_path.exists(),
        "content_keys": list(certificate.keys()),
    },
    "result": {
        "strict_core_promotion_allowed": False,
        "theorem_level_promotion_allowed": False,
        "full_closure_promotion_allowed": False,
        "strict_core_frontier_changed": False,
    },
    "residual_frontier": certificate["strict_core_frontier_unchanged"],
    "hard_limits": certificate["forbidden_promotions"],
    "next_step": "AX8",
}

summary_path = out_dir / "ax7_axiom_lane_anti_overclaim_boundary_audit_summary.json"
summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(summary_path)
