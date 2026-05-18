#!/usr/bin/env python3
"""P1971 S921 strict Higgs branch-lock / invertibility lemma probe.

P1970 identified the remaining PO2 leftover as delta_H.  This probe asks the
next honest question: can the current strict/PO3 exports force delta_H=0?

It proves a conditional algebraic branch-lock lemma: if the Higgs-row map is
strictly monotone/invertible on the admissible interval (encoded by a positive
lower bound for the effective derivative), then equal source/curvature/kinetic
Higgs EOM rows imply delta_H=0.  It then audits current PO3 artifacts and finds
that this positivity/invertibility premise is not exported.  Therefore P1971
sharpens the exact missing object but does not close PO2.
"""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1971_s921_strict_higgs_branch_lock_invertibility_lemma_probe.json"

GREP_FILES = [
    "P1970_STRICT_HIGGS_YUKAWA_CURVATURE_VARIATIONAL_NORMAL_FORM_ATTEMPT_PACKET.md",
    "P1965_STRICT_DELTA_BG_YF_EOM_NORMAL_FORM_EXTRACTION_NONAVAILABILITY_PACKET.md",
    "P1963_STRICT_PO3_DOUBLE_RUN_MACHINE_CHECKER_PACKET.md",
    "P1964_STRICT_PO2_CONDITIONAL_SUFFICIENCY_AND_EOM_GAP_CERTIFICATE_PACKET.md",
]
GREP_TERMS = ["delta_H", "branch", "invert", "PO3", "Higgs", "Omega_unexported", "EOM", "normal-form"]


def load_generated(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode("utf-8")).hexdigest()


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else "MISSING"


def grep_evidence() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name in GREP_FILES:
        path = ROOT / name
        if not path.exists():
            rows.append({"path": name, "present": False, "matches": []})
            continue
        matches = []
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            if any(term.lower() in line.lower() for term in GREP_TERMS):
                matches.append({"line": lineno, "text": line.strip()[:240]})
            if len(matches) >= 10:
                break
        rows.append({"path": name, "present": True, "sha256": file_sha(path), "matches": matches})
    return rows


def conditional_branch_lock_lemma() -> dict[str, Any]:
    H_A, H_B, t = sp.symbols("H_A H_B t", real=True)
    mu2, lam, xi, R = sp.symbols("mu_H_sq lambda_H xi_H R", real=True)
    m_min = sp.symbols("m_eff_min", positive=True)

    H_t = H_B + t * (H_A - H_B)
    fprime = sp.expand(mu2 + xi * R + 3 * lam * H_t**2)
    integral = sp.integrate(fprime, (t, 0, 1))
    f_A_minus_f_B = sp.expand((mu2 + xi * R) * (H_A - H_B) + lam * (H_A**3 - H_B**3))
    mean_value_identity_residual = sp.simplify(f_A_minus_f_B - (H_A - H_B) * integral)

    return {
        "higgs_force_map": "F(H;R)=mu_H_sq*H + lambda_H*H**3 + xi_H*R*H",
        "derivative_along_segment": sp.sstr(fprime),
        "integral_derivative_0_1": sp.sstr(sp.factor(integral)),
        "difference_identity": "F(H_A;R)-F(H_B;R)=(H_A-H_B)*Integral_0^1 F'(H_B+t(H_A-H_B);R) dt",
        "identity_residual_sympy_zero": bool(mean_value_identity_residual == 0),
        "conditional_premise_needed": "F'(H_B+t(H_A-H_B);R) >= m_eff_min > 0 for all t in [0,1] on the admissible branch interval",
        "conditional_conclusion": "If EOM_A-EOM_B=0 with equal source/curvature/kinetic terms and the positive derivative premise holds, then H_A-H_B=0.",
        "proof_sketch_machine_form": "0=(H_A-H_B)*I and I>=m_eff_min>0, hence H_A-H_B=0 over reals.",
        "strictly_positive_symbol": str(m_min),
    }


def audit_current_exports(p1963: dict[str, Any], p1940: dict[str, Any]) -> dict[str, Any]:
    joined = json.dumps({"p1963": p1963, "p1940": p1940}, sort_keys=True).lower()
    required_terms = ["mu_h", "lambda_h", "xi_h", "m_eff", "invert", "monotone", "higgs branch"]
    found_terms = {term: (term in joined) for term in required_terms}
    positivity_exported = any(found_terms[term] for term in ["m_eff", "invert", "monotone", "higgs branch"])
    po3_formal_nonempty = p1963.get("po3_restamp", {}).get("after") == "PASS_MACHINE_CHECKED_FORMAL_DOMAIN_NONEMPTY"
    return {
        "po3_formal_nonempty_pass": po3_formal_nonempty,
        "searched_terms_in_current_po3_exports": found_terms,
        "higgs_invertibility_or_positive_m_eff_exported": bool(positivity_exported),
        "verdict": "MISSING_HIGGS_BRANCH_LOCK_POSITIVITY_PREMISE" if not positivity_exported else "POTENTIAL_PREMISE_FOUND_REVIEW_REQUIRED",
        "meaning": "PO3 proves the formal domain is inhabited, but current exported summaries do not bind that domain to a positive Higgs-row derivative/invertibility condition.",
    }


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    p1970 = load_generated("p1970_s920_strict_higgs_yukawa_curvature_variational_normal_form_attempt.json")
    p1963 = load_generated("p1963_s913_strict_po3_double_run_machine_checker.json")
    p1965 = load_generated("p1965_s915_strict_delta_bg_yf_eom_normal_form_extraction_nonavailability.json")
    p1940 = load_generated("p1940_s890_strict_po3_coeff_inequality_and_machinecheck_transcript_probe.json")

    lemma = conditional_branch_lock_lemma()
    audit = audit_current_exports(p1963, p1940)
    premise_exported = audit["higgs_invertibility_or_positive_m_eff_exported"]
    conditional_lemma_ok = lemma["identity_residual_sympy_zero"]

    out = {
        "packet_id": "P1971",
        "stage_id": "S921",
        "status": "CONDITIONAL_HIGGS_BRANCH_LOCK_LEMMA_PROVED__POSITIVITY_PREMISE_NOT_EXPORTED__PO2_STILL_OPEN",
        "route": "strict_only_higgs_branch_lock_probe",
        "legacy_bridge_used": False,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "repo_grep_evidence": grep_evidence(),
        "input_hashes": {
            "p1970_sha256": digest(p1970),
            "p1963_sha256": digest(p1963),
            "p1965_sha256": digest(p1965),
            "p1940_sha256": digest(p1940),
        },
        "conditional_branch_lock_lemma": lemma,
        "current_export_audit": audit,
        "machine_checks": {
            "conditional_branch_lock_identity_proved": conditional_lemma_ok,
            "po3_formal_domain_nonempty": audit["po3_formal_nonempty_pass"],
            "higgs_positive_derivative_premise_exported": premise_exported,
            "delta_H_forced_zero_in_current_exports": bool(conditional_lemma_ok and premise_exported),
            "full_po2_closed": False,
        },
        "p1970_leftover_restamp": {
            "before": p1970.get("machine_checks", {}).get("sharp_new_blocker_identified"),
            "after": "Exact missing premise: positive Higgs-row derivative / branch-lock invertibility on the admissible PO3 branch interval.",
        },
        "theorem_export": {
            "conditional_positive_statement": "If the admissible branch interval exports F'_H >= m_eff_min > 0 for the Higgs force map, then equal Higgs EOM rows force delta_H=0 and remove the P1970 Yukawa leftover.",
            "current_state_negative_statement": "The current PO3 exports prove formal nonemptiness but do not export the required Higgs invertibility/positive-derivative premise, so delta_H is not forced zero in the current repo state.",
            "not_promoted_to": [
                "full PO2 sufficiency",
                "global background-independence closure",
                "strict-core ToE closure",
                "QW-2191 discharge",
            ],
        },
        "false_pass_guard": "A conditional branch-lock lemma is not a current strict PO2 closure unless its positivity premise is exported on the actual admissible domain.",
        "next_honest_step": "Build P1972: attempt to export a concrete m_eff_min>0 witness from the strict coefficient/domain data; if it fails, freeze delta_H as the exact PO2 normal-form obstruction term.",
        "lay_explanation": "Udowodniliśmy, jaki warunek wystarczyłby do zablokowania różnicy tła Higgsa: równanie Higgsa musi być jednoznaczne/rosnące na dopuszczalnym zakresie. Repo ma dowód, że dopuszczalny zakres nie jest pusty, ale nie ma jeszcze dowodu tej jednoznaczności. Dlatego wiemy dokładniej, czego brakuje, ale nadal nie wolno ogłosić pełnego domknięcia.",
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
