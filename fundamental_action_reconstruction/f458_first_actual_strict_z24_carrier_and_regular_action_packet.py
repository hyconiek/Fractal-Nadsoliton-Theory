#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_I24 = GENERATED / "i_24_v1_index_set.json"
OUT_Z24 = GENERATED / "z_24_v1_group.json"
OUT_TAU = GENERATED / "tau_z24_v1_regular_action_on_i_24_v1.json"
OUT_SUMMARY = GENERATED / "f458_first_actual_strict_z24_carrier_and_regular_action_packet_summary.json"


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    as_of = "2026-03-15"

    i24 = {
        "object": "I_24_v1",
        "status": "actual_exported_strict_finite_index_set_object",
        "as_of": as_of,
        "intent": (
            "Provide an explicit strict 24-slot index-set object {0..23} to serve as a typed carrier for cautious "
            "scope-extension probes beyond the physical n=12 scaffold, without implying any Z_12 identification, "
            "phase embedding, theta export, or QW-2191 discharge."
        ),
        "model": {"set": "{0,1,...,23}", "cardinality": 24},
        "occurs_as": {
            "scope_extension_probe_carrier": "Used by probe-level Zn scope-extension work (e.g. P461 scan; P462 Z24 mode-index assignment candidate)."
        },
        "contracts": {
            "observer_free": True,
            "noncyclic_inputs": True,
            "no_theta_inputs": True,
            "no_qw2191_discharge_implied": True,
        },
        "hard_limits": [
            "Does not identify this carrier with the physical n=12 nad12 scaffold.",
            "Does not define a canonical phase embedding into U(1).",
            "Does not export theta_1/theta_2 and does not claim O(2)-cut closure.",
            "Does not claim QW-2191 discharge or ToE closure.",
        ],
    }

    z24 = {
        "object": "Z_24_v1",
        "status": "actual_exported_strict_finite_group_object",
        "as_of": as_of,
        "intent": (
            "Export a typed cyclic group object of order 24 on the carrier I_24_v1, so scope-extension work may reference an explicit "
            "Z_24 action instead of an untyped integer range."
        ),
        "carrier": {"index_set_object": "I_24_v1", "set": "{0,1,...,23}", "cardinality": 24},
        "group_law": {
            "name": "addition_mod_24",
            "formula": "(a + b) mod 24",
            "identity": 0,
            "inverse": "(-a) mod 24",
            "commutative": True,
        },
        "notes": {
            "cyclic": True,
            "order": 24,
            "generator_nonuniqueness": (
                "A choice of generator is not canonical: Aut(Z_24) ≅ (Z/24Z)^× has 8 units {1,5,7,11,13,17,19,23}. "
                "This object does not fix an embedding Z_24 -> U(1) as canonical."
            ),
        },
        "contracts": {"observer_free": True, "noncyclic_inputs": True, "no_selector_closure_implied": True},
        "hard_limits": [
            "Does not export any canonical phase embedding nor any Berry/holonomy primitive.",
            "Does not cut the QW-2191 O(2) family and does not export any theta supply on n=24.",
            "Does not claim any physical identification with the strict n=12 scaffolds.",
            "Does not claim ToE closure.",
        ],
    }

    tau = {
        "object": "tau_Z24_v1_regular_action_on_I_24_v1",
        "status": "actual_exported_strict_group_action_object",
        "as_of": as_of,
        "intent": (
            "Export the regular action of Z_24_v1 on the 24-slot index set I_24_v1, providing a typed group-action primitive "
            "for scope-extension work without hidden conventions."
        ),
        "group": "Z_24_v1",
        "set": "I_24_v1",
        "definition": {"action_name": "tau_Z24_v1", "formula": "tau_Z24_v1(a,k) := (k + a) mod 24"},
        "properties": {"transitive": True, "free": True, "regular": True},
        "contracts": {
            "observer_free": True,
            "noncyclic_inputs": True,
            "no_theta_inputs": True,
            "no_qw2191_discharge_implied": True,
        },
        "hard_limits": [
            "Does not define a canonical embedding of Z_24 into U(1).",
            "Does not provide a strict O(2)-cut selector ingredient (QW-2191 remains open).",
            "Does not export any n=24 theta supply and does not claim ToE closure.",
        ],
    }

    summary = {
        "stage": "F458",
        "status": "F458_EXECUTED_FIRST_ACTUAL_STRICT_Z24_CARRIER_AND_REGULAR_ACTION_PACKET_NO_FALSE_PASS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": {
            "index_set_object": str(OUT_I24.relative_to(REPO)),
            "group_object": str(OUT_Z24.relative_to(REPO)),
            "action_object": str(OUT_TAU.relative_to(REPO)),
        },
        "no_false_pass": True,
    }

    OUT_I24.write_text(json.dumps(i24, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_Z24.write_text(json.dumps(z24, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_TAU.write_text(json.dumps(tau, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

