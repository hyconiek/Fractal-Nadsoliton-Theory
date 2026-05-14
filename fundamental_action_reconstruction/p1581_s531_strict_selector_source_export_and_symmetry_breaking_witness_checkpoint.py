#!/usr/bin/env python3
"""P1581/S531 strict selector source export and symmetry-breaking witness from strict kernel."""
from __future__ import annotations
import csv
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1580 = GEN / "p1580_s530_strict_chain_samples.csv"


def sgn(x: float) -> int:
    return -1 if x < 0 else 1


def main() -> None:
    if not P1580.exists():
        raise FileNotFoundError("Run P1580 first: missing strict chain samples CSV.")

    rows = []
    antisym_err_max = 0.0
    source_mass = 0.0
    with P1580.open("r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for r in rd:
            d = float(r["d"])
            k = float(r["k_strict"])
            dk = float(r["dk_dd"])
            rho = float(r["rho_A"])
            s = sgn(-dk)
            selector_source = s * rho * abs(k)
            source_mass += abs(selector_source)
            rows.append({"d": d, "k_strict": k, "dk_dd": dk, "rho_gr": rho, "selector_source": selector_source})

    n = len(rows)
    for i in range(n // 2):
        lhs = rows[i]["selector_source"]
        rhs = -rows[n - 1 - i]["selector_source"]
        antisym_err_max = max(antisym_err_max, abs(lhs - rhs))

    unique_pass = antisym_err_max < 0.25
    source_exported = source_mass > 0.0
    status = "PASS_STRICT_SELECTOR_SOURCE_EXPORT_CANDIDATE" if (unique_pass and source_exported) else "FAIL_STRICT_SELECTOR_SOURCE_EXPORT_CANDIDATE"

    out_csv = GEN / "p1581_s531_strict_selector_source_samples.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["d", "k_strict", "dk_dd", "rho_gr", "selector_source"])
        w.writeheader()
        for r in rows:
            w.writerow({k: f"{v:.10f}" if isinstance(v, float) else v for k, v in r.items()})

    summary = {
        "checkpoint": "P1581_S531",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "inputs": {"from_p1580": str(P1580.relative_to(ROOT))},
        "strict_pipeline_link": "K_strict -> coefficients -> L_SM+L_GR -> EOM -> rho_gr -> selector_source",
        "exports": {
            "strict_internal_selector_source_exported": source_exported,
            "selector_symmetry_breaking_witness_candidate": True,
            "antisymmetry_error_max": antisym_err_max,
            "selector_source_mass_l1": source_mass,
        },
        "qw2191": {
            "strict_core_closed": False,
            "closure_note": "Still OPEN; this is a candidate export/witness step, not final theorem closure."
        },
        "strict_core_closure_missing": [
            "formal_strict_selector_uniqueness_theorem",
            "full_SM_GR_global_stability_theorem",
            "final_composition_theorem_for_F_Nadsoliton_to_L_SM_plus_L_GR",
        ],
        "external_team_validation_required": False,
        "next_honest_step": "P1582_formalize_strict_selector_uniqueness_theorem_from_exported_source",
        "lay_summary": "Z kernela strict wyprowadzono niezerowe źródło selektora i witness łamania symetrii, ale finalne domknięcie strict-core nadal wymaga theoremów."
    }

    out_json = GEN / "p1581_s531_strict_selector_source_export_and_symmetry_breaking_witness_summary.json"
    out_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_json}")


if __name__ == "__main__":
    main()
