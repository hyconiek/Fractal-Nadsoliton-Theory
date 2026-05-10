#!/usr/bin/env python3
"""Emit a non-strict manifest for the experimental kernel-unfreeze lane."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORTS = {
    "qw2038": ROOT / "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2038_DERIVATION_COMPATIBLE_KERNEL_REFREEZE_SCAN.md",
    "qw2039": ROOT / "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2039_DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE.md",
    "qw2041": ROOT / "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2041_CANONICAL_REFROZEN_REPARAMETERIZATION_AUDIT.md",
    "qw2049": ROOT / "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.md",
}
FITTING_NOTE = ROOT / "ANALIZA_FITTINGU_I_TRIKOW_KOMPENSACYJNYCH.md"
OUTPUT = ROOT / "kernel_unfreeze_sandbox/generated/unfreeze_omega_phi_beta_eta_manifest.json"

TUPLE_RE = re.compile(
    r"(?i)(?:selected kernel |kernel )omega/phi/beta/eta:\s*"
    r"([0-9.]+)\s*/\s*([0-9.]+)\s*/\s*([0-9.]+)\s*/\s*([0-9.]+)"
)
CANONICAL_RE = re.compile(
    r"Canonical TeX omega/phi/beta/eta:\s*"
    r"([0-9.]+)\s*/\s*([0-9.]+)\s*/\s*([0-9.]+)\s*/\s*([0-9.]+)"
)
REFROZEN_2041_RE = re.compile(
    r"Refrozen QW-2039 omega/phi/beta/eta:\s*"
    r"([0-9.]+)\s*/\s*([0-9.]+)\s*/\s*([0-9.]+)\s*/\s*([0-9.]+)"
)
CI_RE = re.compile(
    r"beta CI95:\s*\[([0-9.]+),\s*([0-9.]+)\].*?eta CI95:\s*\[([0-9.]+),\s*([0-9.]+)\]",
    re.S,
)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def parse_match(match: re.Match[str] | None) -> dict[str, float]:
    if not match:
        raise ValueError("Could not parse omega/phi/beta/eta tuple")
    omega, phi, beta, eta = (float(group) for group in match.groups())
    return {"omega": omega, "phi": phi, "beta": beta, "eta": eta}


def parse_tuple(text: str) -> dict[str, float]:
    return parse_match(TUPLE_RE.search(text))


def midpoint(lhs: dict[str, float], rhs: dict[str, float]) -> dict[str, float]:
    return {
        key: round((lhs[key] + rhs[key]) / 2.0, 6)
        for key in ("omega", "phi", "beta", "eta")
    }


def distance(lhs: dict[str, float], rhs: dict[str, float]) -> float:
    return math.sqrt(
        sum((lhs[key] - rhs[key]) ** 2 for key in ("omega", "phi", "beta", "eta"))
    )


def parse_ci(text: str) -> dict[str, list[float]]:
    match = CI_RE.search(text)
    if not match:
        return {}
    beta_lo, beta_hi, eta_lo, eta_hi = (float(group) for group in match.groups())
    return {"beta": [beta_lo, beta_hi], "eta": [eta_lo, eta_hi]}


def main() -> None:
    qw2038_text = read_text(REPORTS["qw2038"])
    qw2039_text = read_text(REPORTS["qw2039"])
    qw2041_text = read_text(REPORTS["qw2041"])
    qw2049_text = read_text(REPORTS["qw2049"])
    fitting_note_text = read_text(FITTING_NOTE)

    qw2038 = parse_tuple(qw2038_text)
    qw2039 = parse_tuple(qw2039_text)
    canonical = parse_match(CANONICAL_RE.search(qw2041_text))
    qw2041_refrozen = parse_match(REFROZEN_2041_RE.search(qw2041_text))
    qw2049 = parse_tuple(qw2049_text)
    ci = parse_ci(qw2049_text)

    manifest = {
        "lane_status": "NONSTRICT_FIT_DERIVED_PARAMETER_UNFREEZE_SANDBOX",
        "branch": "experimental/unfreeze-omega-phi-beta-eta",
        "root": str(ROOT),
        "provenance": {
            "qw2038_best_row": qw2038,
            "qw2039_refrozen": qw2039,
            "qw2041_canonical": canonical,
            "qw2041_refrozen_side": qw2041_refrozen,
            "qw2049_working_tuple": qw2049,
            "qw2049_micro_ci95": ci,
            "fitting_note_present": "FITTING" in fitting_note_text.upper(),
        },
        "interpretation": {
            "admitted": "Experimental unfreeze is allowed only because the tuple has scan/gate provenance and this lane stays explicitly non-strict.",
            "forbidden": [
                "strict discharge",
                "T176 export",
                "QW-2191 solved",
                "legacy/strict bridge theorem",
            ],
        },
        "search_windows": {
            "local_continuity": {
                "omega": [0.14575, 0.22575],
                "phi": [0.12250, 0.20250],
                "beta": [0.85, 1.15],
                "eta": [1.6, 2.0],
            },
            "canonical_tension": {
                "omega": [0.18575, round(math.pi / 4.0, 6)],
                "phi": [0.16250, round(math.pi / 6.0, 6)],
                "beta": [0.01, 1.0],
                "eta": [1.0, 1.8],
            },
        },
        "seed_candidates": {
            "q2049_hold": qw2049,
            "q2039_refrozen": qw2039,
            "local_beta_midpoint": midpoint(qw2039, qw2049),
            "canonical_bridge_midpoint": midpoint(canonical, qw2049),
            "canonical_refrozen_midpoint": midpoint(qw2041_refrozen, qw2049),
        },
        "distances": {
            "q2049_to_q2039": round(distance(qw2049, qw2039), 6),
            "q2049_to_canonical": round(distance(qw2049, canonical), 6),
            "q2049_to_q2041_refrozen": round(distance(qw2049, qw2041_refrozen), 6),
            "q2039_to_canonical": round(distance(qw2039, canonical), 6),
        },
        "score_axes": [
            "operational_continuity",
            "canonical_tension",
            "selector_relevance",
            "anti_fitting_discipline",
        ],
        "stop_rules": [
            "better fit without selector signal",
            "compensatory freedom explosion",
            "strict-language drift",
            "silent legacy-role transfer",
        ],
    }

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    print("Wrote:", OUTPUT)
    print("QW-2049 working tuple:", qw2049)
    print("QW-2039 refrozen tuple:", qw2039)
    print("QW-2041 refrozen side:", qw2041_refrozen)
    print("Canonical tuple:", canonical)
    print("Distance q2049->canonical:", manifest["distances"]["q2049_to_canonical"])


if __name__ == "__main__":
    main()
