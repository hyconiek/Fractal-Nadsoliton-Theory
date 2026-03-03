#!/usr/bin/env python3
"""
QW-2007: Qualitative benchmark of FIN against major physics programs.

This is a transparent, qualitative scoring framework (not a statistical fit).
Scores are heuristic 0-10 values on explicitly defined axes.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2007_qualitative_theory_benchmark.json"
OUT_MD = ROOT / "RAPORT_QW2007_QUALITATIVE_THEORY_BENCHMARK_EN_PL.md"

AXES = {
    "empirical_coverage": {
        "weight": 0.30,
        "description_en": "How strongly the framework matches existing high-precision data.",
        "description_pl": "Jak silnie teoria zgadza sie z obecnymi danymi wysokiej precyzji.",
    },
    "mathematical_rigor": {
        "weight": 0.25,
        "description_en": "Formal clarity, internal derivational quality, and control of approximations.",
        "description_pl": "Rygor formalny, jakosc wyprowadzen i kontrola przyblizen.",
    },
    "methodological_consistency": {
        "weight": 0.20,
        "description_en": "Stability of protocol across studies (no mixed standards, no silent retunes).",
        "description_pl": "Stabilnosc metodologii miedzy badaniami (bez mieszania standardow i cichych retune).",
    },
    "falsifiability": {
        "weight": 0.15,
        "description_en": "How directly the framework can be tested and potentially falsified.",
        "description_pl": "Jak bezposrednio teoria daje sie testowac i falsyfikowac.",
    },
    "toe_closure_readiness": {
        "weight": 0.10,
        "description_en": "Readiness to claim full unified closure of fundamental sectors.",
        "description_pl": "Gotowosc do claimu pelnego domkniecia sektorow fundamentalnych.",
    },
}

THEORIES = {
    "FIN_Nadsoliton_current_state": {
        "scores": {
            "empirical_coverage": 4.0,
            "mathematical_rigor": 4.0,
            "methodological_consistency": 4.5,
            "falsifiability": 8.0,
            "toe_closure_readiness": 4.0,
        },
        "note_en": "Strong exploratory campaign and falsification intent, but unresolved closure/audit issues remain.",
        "note_pl": "Silna kampania eksploracyjna i falsyfikacyjna, ale pozostaja nierozwiazane luki domkniecia i audytu.",
    },
    "Standard_Model_plus_GR": {
        "scores": {
            "empirical_coverage": 9.5,
            "mathematical_rigor": 9.0,
            "methodological_consistency": 9.0,
            "falsifiability": 9.0,
            "toe_closure_readiness": 6.0,
        },
        "note_en": "Best empirical benchmark; not a full quantum-gravity unification.",
        "note_pl": "Najlepszy benchmark empiryczny; nie jest pelna unifikacja kwantowej grawitacji.",
    },
    "String_Theory_landscape": {
        "scores": {
            "empirical_coverage": 2.0,
            "mathematical_rigor": 8.0,
            "methodological_consistency": 7.0,
            "falsifiability": 3.0,
            "toe_closure_readiness": 5.0,
        },
        "note_en": "Very rich formal machinery, limited direct empirical discrimination so far.",
        "note_pl": "Bardzo bogata formalnie, ale dotad ograniczona bezposrednia rozroznialnosc empiryczna.",
    },
    "Loop_Quantum_Gravity": {
        "scores": {
            "empirical_coverage": 2.0,
            "mathematical_rigor": 6.5,
            "methodological_consistency": 6.0,
            "falsifiability": 4.0,
            "toe_closure_readiness": 4.0,
        },
        "note_en": "Constructive quantum-gravity program with limited precision-contact to data.",
        "note_pl": "Konstruktywny program kwantowej grawitacji, ale z ograniczonym kontaktem precyzyjnym z danymi.",
    },
    "Asymptotic_Safety_gravity": {
        "scores": {
            "empirical_coverage": 2.0,
            "mathematical_rigor": 7.0,
            "methodological_consistency": 6.5,
            "falsifiability": 4.0,
            "toe_closure_readiness": 4.0,
        },
        "note_en": "Strong RG-based logic, but still not a broadly confirmed full ToE framework.",
        "note_pl": "Silna logika RG, ale nadal nie jest szeroko potwierdzona jako pelne ramy ToE.",
    },
    "Causal_Dynamical_Triangulations": {
        "scores": {
            "empirical_coverage": 1.5,
            "mathematical_rigor": 6.0,
            "methodological_consistency": 6.0,
            "falsifiability": 4.0,
            "toe_closure_readiness": 3.0,
        },
        "note_en": "Interesting non-perturbative geometry results, weak direct particle-physics closure.",
        "note_pl": "Ciekawe wyniki geometrii nieperturbacyjnej, slabe domkniecie sektoru czastek.",
    },
    "Causal_Set_program": {
        "scores": {
            "empirical_coverage": 1.0,
            "mathematical_rigor": 5.5,
            "methodological_consistency": 5.5,
            "falsifiability": 3.5,
            "toe_closure_readiness": 3.0,
        },
        "note_en": "Conceptually sharp discreteness idea, still early in full-phenomenology closure.",
        "note_pl": "Koncepcyjnie ostra idea dyskretnosci, ale wczesny etap pelnego domkniecia fenomenologii.",
    },
}


def weighted_total(scores: Dict[str, float]) -> float:
    return sum(scores[k] * AXES[k]["weight"] for k in AXES)


def main() -> None:
    rows: List[dict] = []
    for theory, payload in THEORIES.items():
        total = weighted_total(payload["scores"])
        rows.append(
            {
                "theory": theory,
                "scores": payload["scores"],
                "weighted_total": total,
                "note_en": payload["note_en"],
                "note_pl": payload["note_pl"],
            }
        )

    rows_sorted = sorted(rows, key=lambda x: x["weighted_total"], reverse=True)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "type": "qualitative_benchmark",
        "scale": "0_to_10_per_axis",
        "axes": AXES,
        "rows_ranked": rows_sorted,
        "calibration_anchors": {
            "internal_fin_audits": {
                "QW2006": "PRE1700_METHODOLOGY_MULTI_REGIME_PARTIALLY_RIGOROUS_NOT_FULLY_CLOSED",
                "QW1724": "GW_PIPELINE_METHOD_HIGH_RISK_INCONCLUSIVE",
                "QW1960": "DERIVATION_CONTAINS_MATERIAL_ERRORS_AND_CIRCULAR_STEPS",
            },
            "external_baseline_note": "SM+GR treated as current empirical gold standard in tested domains.",
        },
        "verdict": "FIN_IS_COMPETITIVE_AS_EXPLORATORY_FALSIFIABLE_PROGRAM_BUT_NOT_YET_ABOVE_SM_GR_IN_OVERALL_QUALITY",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines: List[str] = []
    lines.append("# RAPORT QW-2007: QUALITATIVE THEORY BENCHMARK (EN/PL)")
    lines.append("")
    lines.append(f"- Date UTC: {out['generated_utc']}")
    lines.append(f"- Verdict: **{out['verdict']}**")
    lines.append("")
    lines.append("## EN: Method")
    lines.append("- This is a qualitative benchmark with explicit axes and weights.")
    lines.append("- Score scale: 0-10 per axis, then weighted total.")
    lines.append("- Not a proof; a structured comparative snapshot.")
    lines.append("")
    lines.append("## PL: Metoda")
    lines.append("- To benchmark jakosciowy z jawnymi osiami i wagami.")
    lines.append("- Skala: 0-10 na os, potem suma wazona.")
    lines.append("- To nie jest dowod, tylko uporzadkowany obraz porownawczy.")
    lines.append("")
    lines.append("## Axes")
    for ax, desc in AXES.items():
        lines.append(f"- `{ax}` (w={desc['weight']:.2f}): {desc['description_en']}")
    lines.append("")
    lines.append("## Ranking (weighted total)")
    lines.append("| Rank | Theory | Total | Empirical | Math | Method | Falsif. | ToE-ready |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|")

    for idx, r in enumerate(rows_sorted, start=1):
        s = r["scores"]
        lines.append(
            "| "
            + f"{idx} | {r['theory']} | {r['weighted_total']:.2f} | "
            + f"{s['empirical_coverage']:.1f} | {s['mathematical_rigor']:.1f} | "
            + f"{s['methodological_consistency']:.1f} | {s['falsifiability']:.1f} | {s['toe_closure_readiness']:.1f} |"
        )

    lines.append("")
    lines.append("## FIN Position (EN)")
    for r in rows_sorted:
        if r["theory"] == "FIN_Nadsoliton_current_state":
            lines.append(f"- Total: **{r['weighted_total']:.2f}**")
            lines.append(f"- Note: {r['note_en']}")
            break

    lines.append("")
    lines.append("## Pozycja FIN (PL)")
    for r in rows_sorted:
        if r["theory"] == "FIN_Nadsoliton_current_state":
            lines.append(f"- Wynik laczny: **{r['weighted_total']:.2f}**")
            lines.append(f"- Komentarz: {r['note_pl']}")
            break

    lines.append("")
    lines.append("## Internal Anchors Used")
    lines.append("- QW2006: pre-1700 methodology not fully closed.")
    lines.append("- QW1724: GW pipeline high-risk/inconclusive.")
    lines.append("- QW1960: mass derivation contains material errors/circularity.")
    lines.append("")
    lines.append("## Artifacts")
    lines.append(f"- JSON: `{OUT_JSON.name}`")

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2007] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2007] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2007] Verdict: {out['verdict']}")


if __name__ == "__main__":
    main()
