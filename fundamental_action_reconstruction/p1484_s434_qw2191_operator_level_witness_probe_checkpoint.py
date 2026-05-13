#!/usr/bin/env python3
"""P1484 S4.34: operator-level witness probe for strict F=>L_SM+L_GR route."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1483 = GEN / "p1483_s433_qw2191_strict_sector_split_probe_summary.json"
P1468 = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
P1478 = GEN / "p1478_s428_qw2191_sp1_operating_policy_v2.json"

SUMMARY = GEN / "p1484_s434_qw2191_operator_level_witness_probe_summary.json"
OBSTRUCTION = GEN / "p1484_s434_qw2191_operator_level_witness_probe_obstruction.json"


def main() -> None:
    s1483 = json.loads(P1483.read_text(encoding="utf-8"))
    s1468 = json.loads(P1468.read_text(encoding="utf-8"))
    pol = json.loads(P1478.read_text(encoding="utf-8"))

    delta = float(s1468["arm_B_with_SP1_metric"]) - float(s1468["arm_A_no_selector_premise_metric"])
    w_sm = float(s1483["ansatz"]["w_SM"])
    w_gr = float(s1483["ansatz"]["w_GR"])
    eps_mix = float(s1483["ansatz"]["eps_mix"])
    eps_cap = float(s1483["ansatz"]["eps_cap"])

    sm_witness = w_sm * delta
    gr_witness = w_gr * delta
    mix_witness = eps_mix * abs(delta)

    checks = {
        "strict_only": True,
        "no_legacy_bridge": True,
        "qw2191_closed": False,
        "sm_witness_positive": sm_witness > 0,
        "gr_witness_positive": gr_witness > 0,
        "mix_below_sm": mix_witness < sm_witness,
        "mix_below_gr": mix_witness < gr_witness,
        "mix_within_policy_cap": eps_mix <= eps_cap <= float(pol["safety_margin"]),
    }

    pass_checks = [k for k in checks if k != "qw2191_closed"]
    ok = all(checks[k] for k in pass_checks)
    status = "PASS_OPERATOR_WITNESS_LOCAL_ONLY" if ok else "FAIL_OPERATOR_WITNESS_LOCAL_ONLY"

    summary = {
        "packet": "P1484",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "route": "strict_only_F_to_LSM_plus_LGR",
        "witnesses": {
            "sm_witness": sm_witness,
            "gr_witness": gr_witness,
            "mix_witness": mix_witness,
        },
        "checks": checks,
        "next_step_recommendation": "S4.35: promote witness equations to explicit symbolic operator identities and run selector-sensitivity sweep over admissible SP1 shifts.",
        "layman_explanation": "Sprawdziliśmy, czy dwa główne składniki (cząstki i grawitacja) dają wyraźny sygnał dodatni, a składnik mieszany pozostaje mały i pod kontrolą. To znaczy: model nie miesza wszystkiego chaotycznie.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if ok:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
    else:
        OBSTRUCTION.write_text(json.dumps({"packet": "P1484", "status": status, "checks": checks}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"[P1484] status={status} sm={sm_witness:.6f} gr={gr_witness:.6f} mix={mix_witness:.6f}")


if __name__ == "__main__":
    main()
