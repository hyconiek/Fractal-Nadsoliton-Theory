#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_VALUE = GENERATED / "delta_d_sigma_int_positive_window_step_strict_provenance_v1.json"
OUT_SUMMARY = (
    GENERATED
    / "f328_first_actual_strict_sigma_int_positive_window_delta_d_step_strict_provenance_source_upgrade_packet_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


@dataclass(frozen=True)
class KernelTuple:
    omega: float
    phi: float
    beta: float
    eta: float


def corridor_local_window(k: KernelTuple) -> dict[str, float]:
    delta_barrier = (math.pi / 2.0) - abs(k.phi)
    eps_local = 0.5 * delta_barrier
    d_local = eps_local / k.omega
    delta_max = d_local / 11.0
    return {
        "delta_barrier": float(delta_barrier),
        "eps_local": float(eps_local),
        "d_local": float(d_local),
        "delta_max": float(delta_max),
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q2190 = load_json(IN_QW2190)
    kpars = q2190["kernel"]
    k = KernelTuple(
        omega=float(kpars["omega"]),
        phi=float(kpars["phi"]),
        beta=float(kpars["beta"]),
        eta=float(kpars["eta"]),
    )

    corridor = corridor_local_window(k)
    delta_max = float(corridor["delta_max"])

    # Explicit strict-side premise: corridor saturation. This is not strict-derived uniqueness.
    delta_d = float(delta_max)

    obj = {
        "object": "delta_d_sigma_int_positive_window_step_strict_provenance_v1",
        "status": "actual_exported_strict_provenance_value_object__premise_based",
        "as_of": "2026-03-11",
        "intent": (
            "Provide an explicit, observer-free, noncyclic provenance source for the positive-window corridor step "
            "delta_d used by the strict sigma-int → theta candidate lane (T119), removing silent convention laundering."
        ),
        "inputs": {
            "kernel_tuple_source": str(IN_QW2190.relative_to(REPO)),
            "corridor_spec": "T119",
            "target_spec": "T158",
            "carrier_packet": "F328",
        },
        "kernel": {
            "omega": k.omega,
            "phi": k.phi,
            "beta": k.beta,
            "eta": k.eta,
        },
        "corridor": corridor,
        "definition": {
            "delta_max": "d_local/11 from T119 (computed from strict kernel tuple)",
            "delta_d_sigma_int_positive_window_step_strict_provenance_v1": "delta_max",
            "contract": "0 < delta_d <= delta_max",
            "note": "Premise-based corridor saturation; not a strict-derived uniqueness claim.",
        },
        "value": delta_d,
        "provenance": {
            "classification": "strict_source_upgraded",
            "method": "explicit_strict_side_premise_corridor_saturation",
            "premise_id": "P_delta_d_sigma_int_positive_window_corridor_saturation_v1",
            "note": "delta_max is computed from the strict kernel tuple via the T119 corridor; the choice delta_d := delta_max is an explicit premise.",
        },
        "hard_limits": [
            "Does not claim strict derivation or uniqueness of delta_d.",
            "Does not export actual strict-core theta_1/theta_2.",
            "Does not discharge object-support above the export-map object (N395/T130).",
            "Does not imply admissible S_sel_int nor strict-core selector closure.",
            "Does not claim QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F328",
        "object": obj["object"],
        "value": float(obj["value"]),
        "delta_max": float(delta_max),
        "classification": obj["provenance"]["classification"],
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    OUT_VALUE.write_text(
        json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

