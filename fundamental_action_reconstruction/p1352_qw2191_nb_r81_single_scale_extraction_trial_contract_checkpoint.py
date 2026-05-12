#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


def _require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing required {label}: {path}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "P1352: freeze auditable contract for NB-QW2191-R8.1 single-scale "
            "SM/GR extraction trial under strict-only governance."
        )
    )
    ap.add_argument(
        "--p1347",
        type=Path,
        default=FAR / "P1347_STRICT_HOST_LEVEL_IDENTIFICATION_PACKET_EN_PL.md",
    )
    ap.add_argument(
        "--p1348",
        type=Path,
        default=FAR / "P1348_SINGLE_GLOBAL_CLOSURE_THEOREM_PACKET_EN_PL.md",
    )
    ap.add_argument(
        "--p1349",
        type=Path,
        default=FAR / "P1349_EXTERNAL_BLIND_AUDIT_PACKET_EN_PL.md",
    )
    ap.add_argument("--scale", default="MU_STAR_DECLARED_R81_V1")
    ap.add_argument("--scheme", default="MSBAR_DECLARED_R81_V1")
    ap.add_argument(
        "--out",
        type=Path,
        default=GEN / "p1352_qw2191_nb_r81_single_scale_extraction_trial_contract_summary.json",
    )
    a = ap.parse_args()

    _require_file(a.p1347, "P1347")
    _require_file(a.p1348, "P1348")
    _require_file(a.p1349, "P1349")

    contract = {
        "packet": "P1352",
        "as_of": "2026-05-12",
        "lane": "NB_QW2191_STRICT_ONLY",
        "trial": "NB-QW2191-R8.1-SINGLE-SCALE-EXTRACTION-TRIAL-V1",
        "prerequisites": {
            "p1347_host_identification": "PRESENT",
            "p1348_declared_scope_global_closure": "PRESENT",
            "p1349_external_blind_audit_protocol": "PRESENT",
        },
        "contract": {
            "scale": a.scale,
            "scheme": a.scheme,
            "targets": ["g1", "g2", "g3", "gravity_effective_observable_1"],
            "required_outputs": [
                "residual_table_public_v1",
                "pass_fail_status_v1",
                "incident_log_v1",
            ],
            "governance": {
                "new_axioms_allowed": False,
                "silent_legacy_role_transfer_allowed": False,
                "blind_audit_required": True,
                "rollback_on_fail": True,
            },
        },
        "next_priority": "R8_1_BLIND_AUDIT_EXECUTION_AND_RESIDUAL_PUBLICATION",
    }

    a.out.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1352] wrote {a.out}")


if __name__ == "__main__":
    main()
