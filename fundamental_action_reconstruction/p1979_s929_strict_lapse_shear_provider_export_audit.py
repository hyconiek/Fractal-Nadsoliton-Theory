#!/usr/bin/env python3
"""P1979 S929 strict lapse/shear provider export audit.

After P1978, the remaining honest route is non-energy-neutral: an actual strict
variational shear/lapse provider with U00=-Q_shear must be exported, or the
current operator bundle must be shown not to export one.  This script performs
the bounded current-export audit only.  It is not a universal variational no-go.
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
OUT = GEN / "p1979_s929_strict_lapse_shear_provider_export_audit.json"

LAPSE_SHEAR_TOKENS = (
    "lapse",
    "n(t)",
    "adm",
    "k_ij",
    "extrinsic",
    "sigma1",
    "sigma2",
    "sigma_1",
    "sigma_2",
    "shear",
    "anisotropic",
    "bianchi",
    "u00",
    "q_shear",
    "delta_n",
)
COVARIANT_CURVATURE_TOKENS = ("r", "r^2", "ricci", "riemann", "g_{gb}", "g_gb", "lambda")


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
    p1978 = load("p1978_s928_strict_energy_neutral_tensor_transport_obstruction.json")

    reg1846 = (((p1846.get("full_lagrangian_non_skeleton_term_sheet") or {}).get("L_total_termwise")) or {})
    reg1907 = p1907.get("full_lagrangian_term_registry_non_skeleton") or {}
    rows = flatten_terms(reg1846, "P1846_L_total_termwise") + flatten_terms(reg1907, "P1907_term_registry")

    sigma1, sigma2, H, dsigma1, dsigma2 = sp.symbols("sigma1 sigma2 H dsigma1 dsigma2", real=True)
    q_shear = sp.factor(sigma1**2 + sigma1 * sigma2 + sigma2**2)
    required_u = [
        sp.factor(-q_shear),
        sp.factor(3 * H * sigma1 + dsigma1 - q_shear),
        sp.factor(3 * H * sigma2 + dsigma2 - q_shear),
        sp.factor(-3 * H * sigma1 - 3 * H * sigma2 - dsigma1 - dsigma2 - q_shear),
    ]
    required_spatial_trace = sp.factor(sum(required_u[1:]))
    q_matrix = sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    scipy_eigs = la.eigvalsh(np.array(q_matrix.tolist(), dtype=float))

    audit_rows: list[dict[str, Any]] = []
    feature_vectors: list[list[float]] = []
    for row in rows:
        text = row["term"].lower()
        lapse_shear_hits = [tok for tok in LAPSE_SHEAR_TOKENS if tok in text]
        curvature_hits = [tok for tok in COVARIANT_CURVATURE_TOKENS if tok in text]
        exact_signature_hits = [tok for tok in ("u00=-q_shear", "u00 = -q_shear", "-q_shear", "sigma1**2 + sigma1*sigma2 + sigma2**2") if tok in text]
        feature_vectors.append([float(tok in text) for tok in LAPSE_SHEAR_TOKENS])
        audit_rows.append(
            {
                **row,
                "lapse_shear_token_hits": lapse_shear_hits,
                "covariant_curvature_token_hits": curvature_hits,
                "exact_required_signature_hits": exact_signature_hits,
                "exports_required_lapse_shear_provider": bool(lapse_shear_hits and exact_signature_hits),
                "decision": "NO_EXPORTED_LAPSE_SHEAR_PROVIDER_CERTIFICATE" if not (lapse_shear_hits and exact_signature_hits) else "REVIEW_REQUIRED",
            }
        )

    feature_matrix = np.array(feature_vectors, dtype=float) if feature_vectors else np.zeros((0, len(LAPSE_SHEAR_TOKENS)))
    direct_provider_rows = [r for r in audit_rows if r["exports_required_lapse_shear_provider"]]
    lapse_shear_hit_count = int(np.sum(feature_matrix)) if feature_matrix.size else 0
    max_row_l2_norm = float(max((la.norm(feature_matrix[i], ord=2) for i in range(feature_matrix.shape[0])), default=0.0))

    p1978_pass = p1978.get("result_kind") == "PASS_BOUNDED_OBSTRUCTION"
    requires_energy_slot = sp.simplify(required_u[0] + q_shear) == 0
    required_provider_not_tracefree = sp.simplify(required_spatial_trace + 3 * q_shear) == 0
    no_exported_provider = len(direct_provider_rows) == 0
    audited_terms = len(audit_rows) > 0
    audit_passed = bool(p1978_pass and requires_energy_slot and required_provider_not_tracefree and no_exported_provider and audited_terms)

    out = {
        "ledger_id": "P1979_S929_STRICT_LAPSE_SHEAR_PROVIDER_EXPORT_AUDIT",
        "packet_id": "P1979",
        "stage_id": "S929",
        "produced_by": "p1979_s929_strict_lapse_shear_provider_export_audit.py",
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1846_present": "full_lagrangian_non_skeleton_term_sheet" in p1846,
            "p1907_present": "full_lagrangian_term_registry_non_skeleton" in p1907,
            "p1978_bounded_obstruction_present": p1978_pass,
        },
        "required_non_energy_neutral_provider_signature": {
            "Q_shear": str(q_shear),
            "required_U_components": [str(expr) for expr in required_u],
            "required_U00": str(required_u[0]),
            "required_spatial_trace": str(required_spatial_trace),
            "requires_lapse_or_energy_slot": bool(requires_energy_slot),
            "required_provider_is_not_tracefree_spatial": bool(required_provider_not_tracefree),
            "Q_eigenvalues_exact": ["1/2", "3/2"],
            "Q_eigenvalues_scipy": [float(ev) for ev in scipy_eigs],
        },
        "audit_scope": {
            "registries": ["P1846 full_lagrangian_non_skeleton_term_sheet.L_total_termwise", "P1907 full_lagrangian_term_registry_non_skeleton"],
            "terms_audited": len(audit_rows),
            "lapse_shear_tokens": list(LAPSE_SHEAR_TOKENS),
            "covariant_curvature_tokens": list(COVARIANT_CURVATURE_TOKENS),
            "scope_warning": "A covariant curvature term may require a future explicit ADM/Bianchi-I variation before it can be accepted or rejected as a provider; absence of an exported certificate is not a universal no-go.",
        },
        "term_audit_rows": audit_rows,
        "feature_replay": {
            "lapse_shear_feature_matrix_shape": list(feature_matrix.shape),
            "total_lapse_shear_token_hits": lapse_shear_hit_count,
            "max_row_l2_norm": max_row_l2_norm,
            "direct_provider_candidate_count": len(direct_provider_rows),
        },
        "gatekeeper_checks": {
            "p1978_energy_neutral_route_closed_boundedly": p1978_pass,
            "required_provider_has_U00_minus_Q": bool(requires_energy_slot),
            "required_provider_has_nonzero_spatial_trace_minus_3Q": bool(required_provider_not_tracefree),
            "no_exported_lapse_shear_provider_certificate_in_current_registries": no_exported_provider,
            "audited_at_least_one_term": audited_terms,
            "bounded_current_export_audit_passed": audit_passed,
        },
        "result_kind": "PASS_CURRENT_EXPORT_NONAVAILABILITY_AUDIT" if audit_passed else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "NON_ENERGY_NEUTRAL_LAPSE_SHEAR_PROVIDER_NOT_EXPORTED",
            "CURRENT_REGISTRY_HAS_COVARIANT_TERMS_BUT_NO_ADM_BIANCHI_PROVIDER_CERTIFICATE",
            "BACKGROUND_INDEPENDENCE_REMAINS_OPEN",
            "GLOBAL_TOE_CLOSURE_STILL_OPEN",
        ],
        "theorem_export": {
            "current_export_statement": "The current P1846/P1907 strict term registries do not export a lapse/shear provider certificate with U00=-Q_shear and spatial trace -3*Q_shear for the P1974 Bianchi-I residual.",
            "not_a_no_go_clause": "This is not a theorem that the covariant curvature terms cannot yield such a contribution after a future explicit ADM/Bianchi-I variation; it only says that no such provider is currently exported.",
            "remaining_required_work": "Perform a real ADM/Bianchi-I variation of the strict covariant gravity sector, including lapse N(t), for R, R^2, Ricci^2, Riemann^2, and G_GB; then test whether the generated E_N equation supplies the required U00=-Q_shear without ad hoc negative-energy matter.",
            "not_licensed": [
                "global no-go theorem for all strict variational completions",
                "background-independence closure",
                "PO2/PO3 closure",
                "BRST or Cutkosky closure",
                "QW-2191 selector discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1979 is only a current-export audit.  It must not be cited as proving that curvature-squared terms cannot generate the needed lapse/shear equation after a full ADM/Bianchi-I variational calculation.",
        "next_honest_step": "Run the explicit ADM/Bianchi-I lapse variation for the strict gravity sector terms R, R^2, Ricci^2, Riemann^2, and G_GB, and compare the resulting E_N shear terms against U00=-Q_shear.",
        "lay_explanation": "Po P1978 wiemy, że sama korekta przestrzenna nie wystarczy.  Teraz sprawdziliśmy, czy obecny wyeksportowany lagranżian strict ma już jawny certyfikat składnika energii/lapse potrzebnego do naprawy.  Nie ma go.  To nie dowodzi, że taki składnik nie wyjdzie z cięższego rachunku wariacyjnego; mówi tylko, że trzeba go naprawdę policzyć, zamiast go zakładać.",
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
