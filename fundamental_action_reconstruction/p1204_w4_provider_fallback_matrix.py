#!/usr/bin/env python3
"""P1204: provider fallback matrix for W4 symbolic backend selection."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def has_module(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


def main() -> None:
    p1203 = json.loads((GEN / "p1203_w4_no_cas_mode_classifier_summary.json").read_text(encoding="utf-8"))

    providers = [
        {"name": "sympy", "available": has_module("sympy"), "priority": 1},
        {"name": "sage", "available": has_module("sageall") or has_module("sage"), "priority": 2},
        {"name": "external_cas_service", "available": False, "priority": 3},
    ]

    available = [p for p in sorted(providers, key=lambda x: x["priority"]) if p["available"]]
    selected = available[0]["name"] if available else None

    out = {
        "packet": "P1204",
        "as_of": "2026-05-10",
        "input_classification": p1203.get("classification"),
        "providers": providers,
        "selected_provider": selected,
        "fallback_mode_active": selected is None,
        "strict_closure_claim_allowed": False,
        "note": "Provider matrix exported for deterministic backend selection.",
    }

    out_path = GEN / "p1204_w4_provider_fallback_matrix_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1204] selected_provider={selected} fallback_mode_active={out['fallback_mode_active']} wrote {out_path}")


if __name__ == "__main__":
    main()
