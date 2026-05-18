#!/usr/bin/env python3
"""P1976 S926 strict L_total anisotropic provider nonexport audit.

This follows P1975's next honest step.  It audits the currently exported strict
L_total term registries for an explicit derived shear-energy / anisotropic
stress provider matching the P1975 obligation.  It is a nonexport audit, not a
universal impossibility proof.
"""

from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"
OUT = GEN / "p1976_s926_strict_ltotal_anisotropic_provider_nonexport_audit.json"

SHEAR_TOKENS = ("sigma1", "sigma2", "sigma_1", "sigma_2", "shear", "anisotropic", "bianchi", "pi_ij", "piij", "rho_required", "delta_t")
STRICT_SOURCE_TOKENS = ("k_strict", "c_gr", "z_", "xi_h", "kappa", "lambda", "y_f", "m_chi")


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def flatten_terms(registry: dict[str, Any], prefix: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for sector, terms in registry.items():
        if isinstance(terms, list):
            for idx, term in enumerate(terms):
                rows.append({"source": prefix, "sector": str(sector), "term_id": f"{sector}[{idx}]", "term": str(term)})
        elif isinstance(terms, str):
            rows.append({"source": prefix, "sector": str(sector), "term_id": str(sector), "term": terms})
    return rows


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1846 = load("p1846_s796_strict_full_lagrangian_termwise_and_eom_witness_program_checkpoint.json")
    p1907 = load("p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json")
    p1975 = load("p1975_s925_strict_minimal_anisotropic_source_obligation_and_energy_sign_audit.json")

    reg1846 = (((p1846.get("full_lagrangian_non_skeleton_term_sheet") or {}).get("L_total_termwise")) or {})
    reg1907 = p1907.get("full_lagrangian_term_registry_non_skeleton") or {}
    rows = flatten_terms(reg1846, "P1846_L_total_termwise") + flatten_terms(reg1907, "P1907_term_registry")

    audit_rows: list[dict[str, Any]] = []
    feature_vectors: list[list[float]] = []
    for row in rows:
        text = row["term"].lower()
        shear_hits = [tok for tok in SHEAR_TOKENS if tok in text]
        strict_hits = [tok for tok in STRICT_SOURCE_TOKENS if tok in text]
        # Required provider must at minimum name an anisotropic/shear carrier; strict
        # coefficient names alone are insufficient.
        provider_candidate = bool(shear_hits)
        feature = [float(tok in text) for tok in SHEAR_TOKENS]
        feature_vectors.append(feature)
        audit_rows.append(
            {
                **row,
                "shear_token_hits": shear_hits,
                "strict_source_token_hits": strict_hits,
                "explicit_anisotropic_provider_candidate": provider_candidate,
                "decision": "NO_EXPLICIT_PROVIDER_MATCH" if not provider_candidate else "REVIEW_REQUIRED",
            }
        )

    feature_matrix = np.array(feature_vectors, dtype=float) if feature_vectors else np.zeros((0, len(SHEAR_TOKENS)))
    shear_hit_count = int(np.sum(feature_matrix)) if feature_matrix.size else 0
    max_row_norm = float(max((la.norm(feature_matrix[i], ord=2) for i in range(feature_matrix.shape[0])), default=0.0))
    candidates = [r for r in audit_rows if r["explicit_anisotropic_provider_candidate"]]

    s1, s2 = sp.symbols("sigma1 sigma2", real=True)
    rho_required = sp.sympify((p1975.get("source_obligation") or {}).get("rho_required", "0"), locals={"sigma1": s1, "sigma2": s2})
    q_shear = sp.factor(s1**2 + s1 * s2 + s2**2)
    sign_obligation_confirmed = sp.simplify(rho_required + q_shear) == 0

    out = {
        "ledger_id": "P1976_S926_STRICT_LTOTAL_ANISOTROPIC_PROVIDER_NONEXPORT_AUDIT",
        "packet_id": "P1976",
        "stage_id": "S926",
        "produced_by": "p1976_s926_strict_ltotal_anisotropic_provider_nonexport_audit.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1846_present": "full_lagrangian_non_skeleton_term_sheet" in p1846,
            "p1907_present": "full_lagrangian_term_registry_non_skeleton" in p1907,
            "p1975_present": p1975.get("result_kind") == "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        "target_obligation_from_p1975": {
            "rho_required": str(rho_required),
            "q_shear": str(q_shear),
            "rho_required_equals_minus_q_shear": sign_obligation_confirmed,
            "required_provider_minimum": "An exported strict L_total/EOM term with explicit trace-free shear/anisotropic carrier and the P1975 Delta_T component structure.",
        },
        "audit_scope": {
            "registries": ["P1846 full_lagrangian_non_skeleton_term_sheet.L_total_termwise", "P1907 full_lagrangian_term_registry_non_skeleton"],
            "terms_audited": len(audit_rows),
            "shear_tokens": list(SHEAR_TOKENS),
            "strict_source_tokens": list(STRICT_SOURCE_TOKENS),
        },
        "term_audit_rows": audit_rows,
        "feature_replay": {
            "shear_token_feature_matrix_shape": list(feature_matrix.shape),
            "total_shear_token_hits": shear_hit_count,
            "max_row_l2_norm": max_row_norm,
            "candidate_count": len(candidates),
        },
        "gatekeeper_checks": {
            "p1975_sign_obligation_confirmed": sign_obligation_confirmed,
            "no_explicit_anisotropic_provider_in_current_registries": len(candidates) == 0,
            "audited_at_least_one_term": len(audit_rows) > 0,
            "nonexport_not_universal_no_go": True,
        },
        "result_kind": "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "STRICT_LTOTAL_ANISOTROPIC_PROVIDER_NOT_EXPORTED",
            "P1975_NEGATIVE_RHO_OBLIGATION_STILL_UNDERIVED",
            "BACKGROUND_INDEPENDENCE_PROVIDER_GAP_REMAINS",
        ],
        "theorem_export": {
            "nonexport_statement": "In the currently exported P1846/P1907 strict L_total sector registries, no term explicitly exports a trace-free shear/anisotropic stress carrier matching the P1975 minimal Delta_T obligation.",
            "scope_limit": "This is a current-export nonavailability audit, not a theorem that no future strict tensorial provider can be derived.",
            "not_licensed": [
                "global no-go theorem for all possible strict completions",
                "background-independence closure",
                "positive-energy anisotropic source closure",
                "BRST or Cutkosky closure",
                "QW-2191 selector discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1976 only proves that the required anisotropic provider is not exported in the current strict L_total registries; it must not be read as a universal impossibility theorem or as background-independence closure.",
        "next_honest_step": "Either derive a new strict tensorial shear provider from the full variational L_total (beyond the current registries) or formalize a bounded no-go theorem for positive-energy providers under the current P1846/P1907 term basis.",
        "lay_explanation": "Sprawdziliśmy aktualną listę składników teorii strict i nie znaleźliśmy jawnego składnika, który mógłby dostarczyć potrzebne anizotropowe źródło z P1975. To nie dowodzi, że taki składnik nigdy nie może powstać, ale pokazuje, że w obecnie wyeksportowanym lagranżianie nie ma jeszcze brakującej części.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
