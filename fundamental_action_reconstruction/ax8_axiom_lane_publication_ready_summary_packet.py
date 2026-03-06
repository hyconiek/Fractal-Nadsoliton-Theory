from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
out_dir = root / "generated"

closure_packet = json.loads((out_dir / "axiom_lane_closure_packet.json").read_text(encoding="ascii"))
boundary_certificate = json.loads((out_dir / "axiom_lane_boundary_certificate.json").read_text(encoding="ascii"))

packet = {
    "lane": "axiom-augmented",
    "publication_status": "publication_ready_summary_only",
    "selector_axiom": closure_packet["selector_axiom"],
    "selector_family": closure_packet["selector_family"],
    "actual_theta": closure_packet["actual_theta"],
    "actual_basis_pair": closure_packet["actual_basis_pair"],
    "actual_orientation_slice": closure_packet["actual_orientation_slice"],
    "bridge_instance": closure_packet["bridge_instance"],
    "robustness": closure_packet["robustness"],
    "compatibility": closure_packet["compatibility"],
    "boundary": {
        "strict_core_status": boundary_certificate["strict_core_status"],
        "promotion_policy": boundary_certificate["promotion_policy"],
        "forbidden_promotions": boundary_certificate["forbidden_promotions"],
    },
    "strict_core_frontier": closure_packet["residual_frontier"],
}

packet_path = out_dir / "axiom_lane_publication_ready_summary_packet.json"
packet_path.write_text(json.dumps(packet, indent=2) + "\n", encoding="ascii")

summary = {
    "step": "AX8",
    "status": "AX8_EXECUTED_AXIOM_LANE_PUBLICATION_READY_SUMMARY_PACKET_NO_FALSE_PASS",
    "goal": "Assemble AX1..AX7 into one publication-ready summary packet while preserving explicit outside-strict-core boundaries.",
    "created_file": {
        "relative_path": "generated/axiom_lane_publication_ready_summary_packet.json",
        "exists_after_step": packet_path.exists(),
        "content_keys": list(packet.keys()),
    },
    "result": {
        "publication_ready_summary_packet_available": True,
        "strict_core_frontier_changed": False,
        "promotion_into_strict_core_allowed": False,
    },
    "residual_frontier": packet["strict_core_frontier"],
    "hard_limits": packet["boundary"]["forbidden_promotions"],
    "next_step": "AX9",
}

summary_path = out_dir / "ax8_axiom_lane_publication_ready_summary_packet_summary.json"
summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(summary_path)
