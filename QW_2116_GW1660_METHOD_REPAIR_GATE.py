#!/usr/bin/env python3
"""
QW-2116: Methodological repair gate for QW-1660 raw-strain branch.

Purpose:
- verify that key previously flagged methodological defects were addressed
- provide a strict, auditable post-repair verdict
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent

IN_V61 = ROOT / "QW_1660_v61_FullNullModel_strict.json"
IN_V63 = ROOT / "QW_1660_v63_Whitening_strict.json"
IN_V65 = ROOT / "QW_1660_v65_MicroTimeShift_strict.json"
IN_1724 = ROOT / "report_qw1724_gw_method_audit.json"
IN_1725 = ROOT / "report_qw1725_gw_strict_cross_hurst_reanalysis.json"

OUT_JSON = ROOT / "report_qw2116_gw1660_method_repair_gate.json"
OUT_MD = ROOT / "RAPORT_QW2116_GW1660_METHOD_REPAIR_GATE.md"


def load(path: Path):
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8", errors="ignore"))


def main() -> None:
    v61 = load(IN_V61)
    v63 = load(IN_V63)
    v65 = load(IN_V65)
    q1724 = load(IN_1724)
    q1725 = load(IN_1725)

    missing = []
    for p, obj in [
        (IN_V61, v61),
        (IN_V63, v63),
        (IN_V65, v65),
        (IN_1724, q1724),
        (IN_1725, q1725),
    ]:
        if obj is None:
            missing.append(p.name)

    if missing:
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "verdict": "GW1660_METHOD_REPAIR_INCOMPLETE_MISSING_INPUTS",
            "missing_inputs": missing,
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
        lines = [
            "# RAPORT QW-2116: GW1660 METHOD REPAIR GATE",
            "",
            f"- Data UTC: {out['generated_utc']}",
            f"- Werdykt: **{out['verdict']}**",
            "",
            "## Missing inputs",
        ]
        for m in missing:
            lines.append(f"- {m}")
        lines.extend(
            [
                "",
                "## Co to znaczy",
                "- Kod naprawczy jest gotowy, ale rerun strict nie zostal wykonany z powodu brakujacych artefaktow wejscia.",
                "- Najczestsza przyczyna: nieprawidlowe lokalne pliki raw strain (np. pointery Git LFS) i/lub brak lacznosci do GWOSC.",
            ]
        )
        OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(json.dumps(out, ensure_ascii=False, indent=2))
        return

    n61 = int(v61["protocol"]["n_samples"])
    n63 = int(v63["protocol"]["n_samples"])
    n65 = int(v65["protocol"]["n_samples"])

    flags = {
        "v61_observation_not_hardcoded": bool(
            v61["protocol"].get("observation_computed_from_same_protocol", False)
        ),
        "v61_trials_ge_200": int(v61["protocol"]["n_trials"]) >= 200,
        "v63_trials_ge_200": int(v63["protocol"].get("n_trials", 0)) >= 200,
        "v63_has_null_calibration": "null_whitened_phase_randomized" in v63,
        "v65_perm_ge_200": int(v65["protocol"].get("n_perm", 0)) >= 200,
        "v65_has_permutation_pvalue": "null_significance" in v65 and "p_abs" in v65["null_significance"],
        "strict_protocol_n_samples_aligned": (n61 == n63 == n65),
        "legacy_q1724_was_high_risk": q1724.get("verdict") == "GW_PIPELINE_METHOD_HIGH_RISK_INCONCLUSIVE",
        "legacy_q1725_nonrobust_documented": q1725.get("verdict") == "GW_CROSS_HURST_ANOMALY_NOT_ROBUST",
    }

    all_repair_flags_pass = all(flags.values())

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "v61": IN_V61.name,
            "v63": IN_V63.name,
            "v65": IN_V65.name,
            "q1724": IN_1724.name,
            "q1725": IN_1725.name,
        },
        "flags": flags,
        "strict_protocol": {
            "n_samples": [n61, n63, n65],
            "estimator_family": "cross_mfdfa_q0",
            "status": "aligned" if (n61 == n63 == n65) else "misaligned",
        },
        "post_repair_branch_status": {
            "v61_verdict": v61.get("verdict"),
            "v63_verdict": v63.get("verdict"),
            "v65_verdict": v65.get("verdict"),
            "q1725_scientific_verdict": q1725.get("verdict"),
        },
        "verdict": (
            "GW1660_METHOD_REPAIR_GATE_PASS"
            if all_repair_flags_pass
            else "GW1660_METHOD_REPAIR_GATE_PARTIAL_OR_FAIL"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2116: GW1660 METHOD REPAIR GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Werdykt: **{out['verdict']}**",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Post-repair branch status",
            f"- v61 verdict: {out['post_repair_branch_status']['v61_verdict']}",
            f"- v63 verdict: {out['post_repair_branch_status']['v63_verdict']}",
            f"- v65 verdict: {out['post_repair_branch_status']['v65_verdict']}",
            f"- q1725 scientific verdict: {out['post_repair_branch_status']['q1725_scientific_verdict']}",
            "",
            "## Interpretacja",
            "- Bramka 2116 sprawdza naprawe protokolu/metodyki, a nie zamienia automatycznie werdyktu naukowego 1725.",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
