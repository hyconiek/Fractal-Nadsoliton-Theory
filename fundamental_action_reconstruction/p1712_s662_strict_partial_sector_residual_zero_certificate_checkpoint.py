#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1710 = GEN / "p1710_s660_strict_nonproxy_gauge_higgs_first_residual_zero_result_checkpoint.json"
IN1711 = GEN / "p1711_s661_strict_nonproxy_fermion_first_residual_zero_result_checkpoint.json"
OUT = GEN / "p1712_s662_strict_partial_sector_residual_zero_certificate_checkpoint.json"


def main() -> None:
    p1710 = json.loads(IN1710.read_text(encoding="utf-8"))
    p1711 = json.loads(IN1711.read_text(encoding="utf-8"))

    gauge_higgs_pass = p1710.get("residual_status") == "PASS_GAUGE_HIGGS_ZERO"
    fermion_pass = p1711.get("residual_status") == "PASS_FERMION_ZERO"

    payload = {
        "checkpoint": "P1712_S662",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full explicit L_total -> sector EL residual checks -> partial-sector certificate",
        "sector_inputs": {
            "gauge_higgs": {
                "source": "P1710_S660",
                "status": p1710.get("residual_status", "UNKNOWN"),
                "flags": p1710.get("residual_zero_flags", {}),
            },
            "fermion": {
                "source": "P1711_S661",
                "status": p1711.get("residual_status", "UNKNOWN"),
                "flags": p1711.get("residual_zero_flags", {}),
            },
        },
        "partial_sector_certificate": {
            "gauge_higgs_zero": gauge_higgs_pass,
            "fermion_zero": fermion_pass,
            "certificate_status": "PASS_PARTIAL_SECTORS_GH_F"
            if gauge_higgs_pass and fermion_pass
            else "FAIL_PARTIAL_SECTOR_CERTIFICATE",
            "coverage": ["gauge", "higgs", "fermion"],
            "missing_for_full_sector_certificate": ["metric", "cross_consistency_bianchi_ward"],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "metric_sector_residual_zero_nonproxy",
            "bianchi_ward_cross_consistency_certificate",
            "global_helmholtz_integrability_nonproxy",
            "brst_nilpotency_and_cohomology_nonproxy",
            "cutkosky_unitarity_full_sector",
            "counterterm_flow_renormalization_closure",
            "background_independence_family_theorem",
            "qw2191_selector_source_or_nonclosure_theorem",
        ],
        "next_honest_step": "Uderzyć w metric_sector_residual_zero_nonproxy jako ostatni brakujący sektor residual-zero, następnie domknąć bianchi_ward_cross_consistency_certificate dla certyfikatu pełnosektorowego.",
        "lay_summary": "Scaliliśmy dwa zaliczone testy (gauge+Higgs i fermiony) w jeden certyfikat częściowy. Do pełnego certyfikatu równań brakuje teraz głównie sektora grawitacyjnego i testów zgodności przekrojowej.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
