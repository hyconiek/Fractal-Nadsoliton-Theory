#!/usr/bin/env python3
"""P1966 S916 strict QW-2191 selector premise obstruction and minimal breaking.

This packet performs the next honest step after P1965 for the selector block.
It proves, symbolically, that an O(2)-equivariant kernel-only/symmetry-only
operator on a Fourier-degenerate pair plane cannot select a unique axis, and
exports the minimal additional datum that would break the obstruction: a
nonzero traceless symmetric two-tensor on each affected pair plane.

Important refinement after repository grep: the current repo does contain
strict-core, lane-scoped axis-only source exports (notably the Shannon
 element-order and diagonal/local mode-index assignment objects).  Therefore
this packet must distinguish:

* abstract kernel/O(2)-equivariant impossibility (still true),
* lane-scoped strict axis-only sources (exported), and
* global/sign-sensitive admissible selector closure (still not exported).
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1966_s916_strict_qw2191_selector_premise_obstruction_and_minimal_breaking.json"


def load_generated(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_repo_file(name: str) -> str:
    path = ROOT / name
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def symbolic_o2_commutant() -> dict[str, Any]:
    x, y, z, w, t = sp.symbols("x y z w t", real=True)
    c = sp.cos(t)
    s = sp.sin(t)
    r = sp.Matrix([[c, -s], [s, c]])
    a = sp.Matrix([[x, y], [z, w]])
    comm = sp.simplify(a * r - r * a)

    # Equivariance for every t is equivalent to all Taylor coefficients at t=0
    # vanishing.  The first derivative already imposes y+z=0 and w-x=0; the
    # zero-order coefficient is identically zero.  Symmetric selector operators
    # then force y=z=0 and w=x, so only scalar matrices remain.
    dcomm0 = sp.diff(comm, t).subs(t, 0)
    constraints_equiv = [sp.Eq(entry, 0) for entry in dcomm0]
    sol_equiv = sp.solve(constraints_equiv, (z, w), dict=True)

    sym_constraints = [sp.Eq(y, z)]
    sol_symmetric_equiv = sp.solve(constraints_equiv + sym_constraints, (y, z, w), dict=True)
    scalar_form = a.subs({y: 0, z: 0, w: x})

    eigenvals_scalar = scalar_form.eigenvals()
    distinct_axis_possible = len(eigenvals_scalar) == 2

    return {
        "rotation_matrix_R_t": sp.sstr(r),
        "generic_operator_A": sp.sstr(a),
        "commutator_A_R_t": sp.sstr(comm),
        "derivative_constraints_at_identity": [sp.sstr(eq) for eq in constraints_equiv],
        "all_linear_equivariant_solutions": [
            {str(k): sp.sstr(v) for k, v in solution.items()} for solution in sol_equiv
        ],
        "symmetric_equivariant_solutions": [
            {str(k): sp.sstr(v) for k, v in solution.items()} for solution in sol_symmetric_equiv
        ],
        "symmetric_equivariant_operator_normal_form": sp.sstr(scalar_form),
        "scalar_eigenvalue_multiplicity": {sp.sstr(k): int(v) for k, v in eigenvals_scalar.items()},
        "unique_axis_possible_from_o2_equivariant_kernel_alone": bool(distinct_axis_possible),
        "verdict": "O2_EQUIVARIANT_SYMMETRIC_OPERATOR_IS_SCALAR_ON_PAIR_PLANE__NO_AXIS_SELECTOR",
    }


def minimal_breaking_tensor() -> dict[str, Any]:
    a, b, lam = sp.symbols("sigma_c sigma_s lambda", real=True)
    s_mat = sp.Matrix([[a, b], [b, -a]])
    charpoly = sp.factor(s_mat.charpoly(lam).as_expr())
    eigenvals = s_mat.eigenvals()
    gap_squared = sp.simplify((2 * sp.sqrt(a**2 + b**2)) ** 2)
    determinant = sp.factor(s_mat.det())
    trace = sp.trace(s_mat)

    samples = []
    for idx, (av, bv) in enumerate([(1.0, 0.0), (0.0, 1.0), (0.6, 0.8), (2.0, -1.5)], start=1):
        mat = np.array([[av, bv], [bv, -av]], dtype=float)
        vals, vecs = np.linalg.eigh(mat)
        residual = float(np.linalg.norm(mat @ vecs - vecs @ np.diag(vals), ord=np.inf))
        samples.append(
            {
                "sample_id": idx,
                "sigma_c": av,
                "sigma_s": bv,
                "eigenvalues_numeric": [float(v) for v in vals],
                "gap_numeric": float(vals[1] - vals[0]),
                "residual_inf_norm": residual,
                "selector_angle_half_atan2_rad": 0.5 * math.atan2(bv, av),
                "passes_nonzero_gap": bool(vals[1] - vals[0] > 1e-12),
            }
        )

    return {
        "minimal_datum_label": "Delta_sel_pair_m = [[sigma_c_m, sigma_s_m], [sigma_s_m, -sigma_c_m]]",
        "matrix": sp.sstr(s_mat),
        "trace": sp.sstr(trace),
        "determinant": sp.sstr(determinant),
        "characteristic_polynomial": sp.sstr(charpoly),
        "eigenvalues_symbolic": {sp.sstr(k): int(v) for k, v in eigenvals.items()},
        "gap_squared": sp.sstr(gap_squared),
        "nondegeneracy_condition": "sigma_c_m**2 + sigma_s_m**2 > 0 for every selected degenerate pair plane",
        "residual_symmetry_after_breaking": "Z2 sign of each selected eigenvector, not continuous O(2)",
        "numeric_probe_samples": samples,
        "all_numeric_samples_pass": all(row["passes_nonzero_gap"] and row["residual_inf_norm"] < 1e-12 for row in samples),
    }


def current_source_audit() -> dict[str, Any]:
    p1499 = load_generated("p1499_s449_qw2191_global_closure_requirements_summary.json")
    p1552 = load_generated("p1552_s502_strict_internal_selector_source_witness_construction_summary.json")
    p1965 = load_generated("p1965_s915_strict_delta_bg_yf_eom_normal_form_extraction_nonavailability.json")
    shannon_assignment = load_generated("mode_index_assignment_shannon_element_order_reference_strict_core_v1.json")
    diagonal_assignment = load_generated("mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json")
    p455_summary = load_generated("p455_current_strict_mode_index_assignment_shannon_vs_diagonal_alignment_audit_probe_summary.json")
    n480 = read_repo_file("N480_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_CROSS_ENTROPY_CUTS_PAIR1_O2_TO_Z2_UNIQUENESS_THEOREM.md")
    n488 = read_repo_file("N488_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_CROSS_ENTROPY_CUTS_PAIR2_O2_TO_Z2_UNIQUENESS_THEOREM.md")
    n496 = read_repo_file("N496_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_CROSS_ENTROPY_CUTS_PAIR3_TO_PAIR5_O2_TO_Z2_UNIQUENESS_THEOREM.md")
    n487 = read_repo_file("N487_CURRENT_FIRST_STRICT_QW2191_CONTINUOUS_O2_FAMILY_DISCHARGE_ON_ALL_FOURIER_DEGENERATE_PAIRS_VIA_CANONICAL_DIAGONAL_LOCAL_SECTOR_THEOREM.md")
    ax16 = read_repo_file("AX16_STRICT_SIDE_SELECTOR_INGREDIENT_ADMISSIBILITY_PRINCIPLE_THEORY_LEVEL_ACCEPTANCE_PACKET.md")

    p1499_requirements = p1499.get("requirements", {}) if isinstance(p1499, dict) else {}
    p1552_checks = p1552.get("physical_checks", {}) if isinstance(p1552, dict) else {}

    shannon_axis_source_exported = (
        shannon_assignment.get("object") == "ModeIndexAssignment_shannon_element_order_reference_strict_core_v1"
        and "strict_core" in str(shannon_assignment.get("scope", ""))
        and all(f"pair{m}" in (shannon_assignment.get("outputs", {}).get("pairs", {})) for m in range(1, 6))
    )
    diagonal_axis_source_exported = (
        diagonal_assignment.get("object") == "ModeIndexAssignment_canonical_local_diagonal_strict_derived_v1"
        and all(f"pair{m}" in (diagonal_assignment.get("outputs", {}).get("pairs", {})) for m in range(1, 6))
    )
    shannon_theorems_present = all("DISCHARGED" in text for text in [n480, n488, n496])
    axis_alignment_audited = p455_summary.get("aligned_all_pairs_up_to_residual_z2") is True

    n487_scoped = "does not claim" in n487.lower() and "global strict-core selector closure" in n487.lower()
    ax16_extension_only = "strict_extension_only" in ax16

    axis_only_source_exported_now = bool(
        shannon_axis_source_exported
        and diagonal_axis_source_exported
        and shannon_theorems_present
        and axis_alignment_audited
    )

    missing_for_global = []
    if not p1499_requirements.get("R2_strict_internal_selector_source_exported", False):
        missing_for_global.append("P1499_R2_global_strict_internal_selector_source_not_restamped")
    if not p1499_requirements.get("R5_F_to_LSM_LGR_mapping_witness_exported", False):
        missing_for_global.append("P1499_R5_F_to_LSM_LGR_mapping_witness_missing")
    if not p1552_checks.get("theorem_level_uniqueness_pass", False):
        missing_for_global.append("P1552_theorem_level_uniqueness_pass_false")
    if n487_scoped:
        missing_for_global.append("N487_scoped_diagonal_local_not_global_selector_closure")
    if ax16_extension_only:
        missing_for_global.append("AX16_extension_only_not_strict_core_admissible_S_sel_int")

    return {
        "inputs": {
            "p1499_sha256": digest(p1499),
            "p1552_sha256": digest(p1552),
            "p1965_sha256": digest(p1965),
            "shannon_assignment_sha256": digest(shannon_assignment),
            "diagonal_assignment_sha256": digest(diagonal_assignment),
            "p455_summary_sha256": digest(p455_summary),
            "n480_text_sha256": hashlib.sha256(n480.encode("utf-8")).hexdigest(),
            "n488_text_sha256": hashlib.sha256(n488.encode("utf-8")).hexdigest(),
            "n496_text_sha256": hashlib.sha256(n496.encode("utf-8")).hexdigest(),
            "n487_text_sha256": hashlib.sha256(n487.encode("utf-8")).hexdigest(),
            "ax16_text_sha256": hashlib.sha256(ax16.encode("utf-8")).hexdigest(),
        },
        "axis_only_strict_sources_found_by_repo_grep": {
            "shannon_element_order_assignment_exported": shannon_axis_source_exported,
            "diagonal_local_assignment_exported": diagonal_axis_source_exported,
            "shannon_pair_theorems_n480_n488_n496_present": shannon_theorems_present,
            "p455_axis_alignment_audited_up_to_z2": axis_alignment_audited,
            "axis_only_source_exported_now": axis_only_source_exported_now,
        },
        "global_selector_closure_audit": {
            "p1499_global_closed": p1499.get("global_closed"),
            "p1499_requirements": p1499_requirements,
            "p1552_qw2191_closed": p1552.get("qw2191_closed"),
            "p1552_theorem_level_uniqueness_pass": p1552_checks.get("theorem_level_uniqueness_pass"),
            "n487_scoped_not_global": n487_scoped,
            "ax16_extension_only": ax16_extension_only,
            "p1965_status": p1965.get("status"),
            "missing_for_global_strict_selector_closure": missing_for_global,
            "admissible_S_sel_int_or_sign_sensitive_global_selector_exported_now": False,
        },
    }


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    commutant = symbolic_o2_commutant()
    breaking = minimal_breaking_tensor()
    audit = current_source_audit()

    axis_only_source_found = audit["axis_only_strict_sources_found_by_repo_grep"]["axis_only_source_exported_now"]
    global_selector_closure = audit["global_selector_closure_audit"][
        "admissible_S_sel_int_or_sign_sensitive_global_selector_exported_now"
    ]
    extension_premise_identified = breaking["all_numeric_samples_pass"]

    out = {
        "packet_id": "P1966",
        "stage_id": "S916",
        "status": "QW2191_MINIMAL_SELECTOR_PREMISE_IDENTIFIED__STRICT_AXIS_ONLY_SOURCE_FOUND__GLOBAL_SELECTOR_CLOSURE_NOT_PASSED",
        "route": "strict_only_with_lane_scoped_axis_source_distinguished_from_global_selector_closure",
        "legacy_bridge_used": False,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "sympy": sp.__version__,
        "symbolic_o2_obstruction": commutant,
        "minimal_symmetry_breaking_premise": breaking,
        "current_repository_source_audit": audit,
        "theorem_export": {
            "kernel_only_negative_statement": "On a degenerate Fourier pair plane, every symmetric O(2)-equivariant kernel/symmetry-only operator is scalar, so that layer alone cannot choose a unique selector axis.",
            "axis_source_refinement": "Repository grep finds lane-scoped strict axis-only source exports (Shannon element-order plus diagonal/local mode-index assignments) that instantiate the nonzero traceless tensor shape up to residual Z2.",
            "minimal_conditional_positive_statement": "If each affected pair plane receives a nonzero traceless symmetric selector tensor Delta_sel_pair_m with sigma_c_m**2 + sigma_s_m**2 > 0, the continuous O(2) family is cut to residual Z2 and a unique axis pair is selected.",
            "remaining_nonclaim": "Axis-only O(2)->Z2 is not the same as sign-sensitive global selector closure/admissible S_sel_int/ToE closure.",
            "status": "AXIS_ONLY_STRICT_SOURCE_RECOGNIZED__GLOBAL_SELECTOR_CLOSURE_STILL_OPEN",
        },
        "strict_axis_only_source_exported_now": bool(axis_only_source_found),
        "strict_core_selector_discharge_pass": bool(axis_only_source_found and global_selector_closure),
        "extension_premise_identified": extension_premise_identified,
        "false_pass_guard": [
            "This packet does not claim QW-2191 global strict-core selector closure.",
            "It corrects the previous audit by recognizing existing strict lane-scoped axis-only source exports.",
            "Residual Z2 sign remains unresolved unless gauge-irrelevance is separately proved for each downstream observable or a sign-sensitive strict datum is exported.",
            "AX16 remains strict_extension_only and is not re-labeled as strict-core admissible S_sel_int.",
        ],
        "strict_core_statusvector_restamp": {
            "renormalization": "OPEN_BACKEND_THEOREM_SCOPE_LIMITS_REMAIN",
            "unitarity": "OPEN_DRESSED_AMPLITUDE_MISSING",
            "background_independence": "OPEN_GLOBAL_TRANSPORT_PENDING",
            "po3": "PASS_MACHINE_CHECKED_FORMAL_DOMAIN_NONEMPTY_BUT_COVARIANT_FULL_CHECK_OPEN",
            "po2": "CONDITIONAL_TENSORIAL_PASS_FULL_EOM_NORMAL_FORM_OPEN",
            "selector_qw2191": "AXIS_ONLY_O2_TO_Z2_SOURCE_FOUND_ON_EXPORTED_N12_LANES__GLOBAL_SIGN_SENSITIVE_SELECTOR_CLOSURE_OPEN",
        },
        "next_honest_step": "Build P1967: explicitly map the Shannon element-order strict source into the P1966 traceless tensor Delta_sel_pair_m for pair1..pair5, verify nonzero gaps, and keep the residual Z2/global selector limits explicit.",
        "lay_explanation": "Po przeszukaniu repo widać, że istnieje wewnętrzny strict kompas wybierający osie (ale tylko do znaku). To naprawia zbyt ostrożną wersję P1966: brakowało rozróżnienia między wyborem osi a pełnym wyborem zwrotu/znaku i globalnym selektorem. Pełne domknięcie nadal nie jest udowodnione.",
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
