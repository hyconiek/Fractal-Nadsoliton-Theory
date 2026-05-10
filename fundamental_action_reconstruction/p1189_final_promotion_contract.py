#!/usr/bin/env python3
"""P1189 final promotion contract combining all key strict gates."""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    p1183 = json.loads((GEN / "p1183_candidate_promotion_gate_summary.json").read_text(encoding="utf-8"))
    p1188 = json.loads((GEN / "p1188_promotion_envelope_check_summary.json").read_text(encoding="utf-8"))
    p1186 = json.loads((GEN / "p1186_post_certification_stability_rate_summary.json").read_text(encoding="utf-8"))

    gate_p1183 = bool(p1183.get("promote_candidate"))
    gate_p1188 = bool(p1188.get("promotion_envelope_pass"))
    stability_rate = float(p1186.get("stability_rate", 0.0))
    gate_p1186 = stability_rate >= 0.95

    contract_pass = gate_p1183 and gate_p1188 and gate_p1186

    out = {
        "packet": "P1189",
        "as_of": "2026-05-10",
        "winner_candidate": p1183.get("winner_candidate"),
        "gates": {
            "p1183_promote_candidate": gate_p1183,
            "p1188_promotion_envelope_pass": gate_p1188,
            "p1186_stability_rate_ge_0_95": gate_p1186,
        },
        "stability_rate": stability_rate,
        "final_promotion_contract_pass": contract_pass,
        "note": "Final contract only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1189_final_promotion_contract_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1189] final_promotion_contract_pass={contract_pass} wrote {out_path}")

if __name__ == "__main__":
    main()
