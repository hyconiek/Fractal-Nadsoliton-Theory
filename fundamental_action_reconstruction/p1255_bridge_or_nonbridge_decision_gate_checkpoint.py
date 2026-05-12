#!/usr/bin/env python3
"""P1255: enforce honest next-step decision gate for legacy->strict bridge question.

This checkpoint explicitly prevents claiming strict-theory closure while the
legacy->strict bridge (or a formal non-bridge theorem) remains unresolved.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1254", type=Path, default=GEN / "p1254_l2_discharge_attempt_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1255_bridge_or_nonbridge_decision_gate_summary.json")
    args = parser.parse_args()

    p1254 = _load_json(args.p1254)

    closure_status = str(p1254.get("theory_closure_status", "OPEN"))
    strict_closure_claim_allowed = bool(p1254.get("strict_closure_claim_allowed", False))
    open_obligations = int(p1254.get("open_obligation_count_after_l2_attempt", 999999))

    bridge_status = "UNRESOLVED"
    nonbridge_theorem_status = "UNRESOLVED"

    next_honest_step = {
        "priority": "HIGHEST",
        "task": "Export either a legacy->strict bridge theorem or a formal non-bridge theorem",
        "classification": "kernel-split-robust",
        "forbidden_shortcut": "No strict-core closure claim from local L1/L2 discharge alone",
    }

    gate_pass = (
        closure_status == "OPEN"
        and (not strict_closure_claim_allowed)
        and open_obligations >= 1
        and bridge_status == "UNRESOLVED"
        and nonbridge_theorem_status == "UNRESOLVED"
    )

    out = {
        "packet": "P1255",
        "as_of": "2026-05-11",
        "input": {"p1254": str(args.p1254)},
        "inherited_state": {
            "theory_closure_status": closure_status,
            "strict_closure_claim_allowed": strict_closure_claim_allowed,
            "open_obligation_count_after_l2_attempt": open_obligations,
        },
        "bridge_status": bridge_status,
        "nonbridge_theorem_status": nonbridge_theorem_status,
        "gate_pass": gate_pass,
        "decision": "PROCEED_TO_BRIDGE_OR_NONBRIDGE_FORMALIZATION_ONLY",
        "next_honest_step": next_honest_step,
        "lay_summary": (
            "Mamy częściowe kroki techniczne, ale kluczowe pytanie pozostaje otwarte: "
            "czy stary i nowy kernel są formalnie równoważne, czy nie. "
            "Dopóki tego nie udowodnimy (albo nie obalimy formalnie), "
            "nie wolno ogłaszać pełnego domknięcia teorii."
        ),
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1255] gate_pass={gate_pass} wrote {args.out}")


if __name__ == "__main__":
    main()
