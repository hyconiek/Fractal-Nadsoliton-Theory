#!/usr/bin/env python3
"""
QW-1724: Methodological audit of GW phase studies (v47-v65).

Goal:
1) Detect internal inconsistencies across already produced GW results.
2) Quantify methodological risk (cherry-picking, protocol drift, mixed baselines).
3) Produce a strict "evidence reliability" verdict.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1724_gw_method_audit.json"
OUT_MD = ROOT / "RAPORT_QW1724_GW_METHOD_AUDIT.md"


def load(path: Path):
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    except Exception:
        return None


def push_issue(issues: list[dict], severity: str, title: str, detail: str) -> None:
    issues.append({"severity": severity, "title": title, "detail": detail})


def severity_points(sev: str) -> int:
    if sev == "critical":
        return 4
    if sev == "major":
        return 2
    return 1


def main() -> None:
    j48 = load(ROOT / "QW_1660_v48_Pure_Raw_CrossDFA.json")
    j50 = load(ROOT / "QW_1660_v50_MonteCarlo_CrossHurst.json")
    j53 = load(ROOT / "QW_1660_v53_PSDSurrogate.json")
    j55 = load(ROOT / "QW_1660_v55_InjectionTest.json")
    j56 = load(ROOT / "QW_1660_v56_ScaleSeparation.json")
    j59 = load(ROOT / "QW_1660_v59_EnvelopeSwap.json")
    j60 = load(ROOT / "QW_1660_v60_ScaleStability.json")
    j61 = load(ROOT / "QW_1660_v61_FullNullModel.json")
    j62 = load(ROOT / "QW_1660_v62_InterObservatory.json")
    j63 = load(ROOT / "QW_1660_v63_Whitening.json")
    j64 = load(ROOT / "QW_1660_v64_Bandpass.json")
    j65 = load(ROOT / "QW_1660_v65_MicroTimeShift.json")

    issues: list[dict] = []
    missing = []
    for name, obj in [
        ("v48", j48),
        ("v50", j50),
        ("v53", j53),
        ("v55", j55),
        ("v56", j56),
        ("v59", j59),
        ("v60", j60),
        ("v61", j61),
        ("v62", j62),
        ("v63", j63),
        ("v64", j64),
        ("v65", j65),
    ]:
        if obj is None:
            missing.append(name)
    if missing:
        push_issue(
            issues,
            "major",
            "Brak kluczowych artefaktow",
            f"Brak plikow: {missing}. Audit oparty na czesci danych.",
        )

    # 1) Contradiction: v50 says mixed model ~0.51, far from 0.31; v55 claims projection 0.23->0.31 works.
    if j50 and j55:
        cross_mixed = j50.get("mixed_signal", {}).get("Cross_H_mean")
        inj = j55.get("Injected_Cross_H_Output")
        if cross_mixed is not None and inj is not None:
            if abs(cross_mixed - 0.31) > 0.10 and abs(inj - 0.31) < 0.03:
                push_issue(
                    issues,
                    "critical",
                    "Sprzecznosc modelu projekcji",
                    f"v50 daje Cross_H={cross_mixed:.3f} (daleko od 0.31), a v55 twierdzi projekcje do {inj:.3f}~0.31.",
                )

    # 2) Scale invariance claim vs scale drift.
    if j56 and j60:
        h_short = j56.get("H_cross_short_scales_1s_to_25s")
        h_long = j56.get("H_cross_long_scales_25s_to_128s")
        vals = []
        for k, v in j60.items():
            if isinstance(v, dict) and "H_cross" in v:
                vals.append(float(v["H_cross"]))
        if h_short is not None and h_long is not None and vals:
            spread = max(vals) - min(vals)
            if spread > 0.50:
                push_issue(
                    issues,
                    "critical",
                    "Brak inwariancji skali",
                    f"Rozrzut H_cross po skalach = {spread:.3f} (od {min(vals):.3f} do {max(vals):.3f}), mimo claimu o inwariancji.",
                )

    # 3) Whitened high H vs micro time-shift no coherent lag.
    if j63 and j65:
        h_white = j63.get("Whitened_Cross_H")
        best_lag_ms = j65.get("Best_Lag_ms")
        best_corr = j65.get("Best_Correlation")
        if h_white is not None and best_lag_ms is not None and best_corr is not None:
            if h_white > 0.60 and abs(best_corr) < 0.02:
                push_issue(
                    issues,
                    "major",
                    "Niespojnosc fazowa",
                    (
                        f"Po whitening H_cross={h_white:.3f} sugeruje silna strukture, ale "
                        f"test opoznienia ma max korelacje {best_corr:.4f} przy {best_lag_ms:.2f} ms (bardzo slabo)."
                    ),
                )

    # 4) Inter-observatory mismatch.
    if j48 and j62:
        h_h1l1 = j48.get("Cross_H1_L1_Pure_H")
        h_h1v1 = j62.get("H1_V1_Cross_H")
        if h_h1l1 is not None and h_h1v1 is not None:
            if abs(h_h1v1 - h_h1l1) > 0.30:
                push_issue(
                    issues,
                    "major",
                    "Brak reprodukcji miedzy parami detektorow",
                    f"H1-L1={h_h1l1:.3f} vs H1-V1={h_h1v1:.3f}; roznica {abs(h_h1v1-h_h1l1):.3f}.",
                )

    # 5) Full null model may be protocol-drifted (hard-coded observation).
    if j61:
        obs = j61.get("Observation")
        ns = j61.get("Config", {}).get("N_samples")
        if obs == 0.311 and ns is not None and int(ns) == 131072:
            push_issue(
                issues,
                "major",
                "Mozliwy drift protokolu null-testu",
                "W v61 obserwacja jest hard-coded jako 0.311 mimo innej dlugosci sygnalu (N=131072).",
            )

    # 6) Band-specific results crossing stationary/non-stationary regimes.
    if j64:
        band_vals = [float(v) for v in j64.values() if isinstance(v, (int, float))]
        if band_vals and (max(band_vals) - min(band_vals) > 0.25):
            push_issue(
                issues,
                "minor",
                "Silna zaleznosc pasmowa",
                f"H_cross po pasmach: min={min(band_vals):.3f}, max={max(band_vals):.3f}.",
            )

    # 7) Envelope swap indicates strong model sensitivity (possible non-identifiability).
    if j59:
        t1 = j59.get("Test1_Independent_Phases_Real_Envelope", {}).get("mean")
        t2 = j59.get("Test2_Real_Phases_Generic_Envelope", {}).get("value")
        if t1 is not None and t2 is not None and abs(t2 - t1) > 0.30:
            push_issue(
                issues,
                "major",
                "Niska identyfikowalnosc mechanizmu",
                f"Zmiana envelope/fazy przesuwa H_cross o {abs(t2-t1):.3f} (Test1={t1:.3f}, Test2={t2:.3f}).",
            )

    total_risk = sum(severity_points(i["severity"]) for i in issues)
    if total_risk <= 4:
        verdict = "GW_PIPELINE_METHOD_HIGH_CONFIDENCE"
    elif total_risk <= 10:
        verdict = "GW_PIPELINE_METHOD_MODERATE_RISK"
    else:
        verdict = "GW_PIPELINE_METHOD_HIGH_RISK_INCONCLUSIVE"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "issues": issues,
        "total_risk_points": int(total_risk),
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1724: GW METHOD AUDIT",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Risk points: {total_risk}",
        f"- Werdykt: **{verdict}**",
        "",
        "## Zidentyfikowane problemy",
    ]
    if not issues:
        lines.append("- Brak wykrytych niespojnosci metodologicznych.")
    else:
        for i, it in enumerate(issues, start=1):
            lines.append(f"- [{i}] ({it['severity']}) {it['title']}: {it['detail']}")
    lines.extend(
        [
            "",
            "## Artefakty",
            f"- JSON szczegolowy: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1724] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1724] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
