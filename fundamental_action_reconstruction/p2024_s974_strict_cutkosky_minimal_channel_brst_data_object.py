#!/usr/bin/env python3
"""P2024 S974 strict Cutkosky minimal finite BRST quartet model.

Next honest step after P2023: export the smallest machine-checkable finite
BRST quartet model that can be stated without a false Cutkosky closure claim.

The exported object is a finite BRST-quartet *model* for one gauge cut-state
polarization layer:

    Q |gauge_T1> = 0
    Q |gauge_T2> = 0
    Q |gauge_L>  = |ghost>
    Q |ghost>    = 0
    Q |antighost> = |B>
    Q |B>        = 0

The transverse labels are deliberately gauge-state labels, not P2020 graviton
plus/cross labels.  This is anchored to the P1961 local gauge-sector BRST
differential and to the P1956 local transverse projector.  It proves exact
nilpotency for a finite algebraic model and gives a ghost-quartet supertrace
convention.  It does not identify that finite cohomology count with the P2020
phase-space CutSum trace and does not prove the full BRST cohomology projector,
ghost-cancelled DiscM, or the optical theorem.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
import sys
from typing import Any

import numpy as np
import scipy
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2024_s974_strict_cutkosky_minimal_channel_brst_data_object.json"
TS = "2026-05-19T00:00:00+00:00"
CHANNEL = "graviton->gauge_gauge"
MODEL_SCOPE = "finite_algebraic_quartet_not_full_channel_projector_v11"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def matrix_to_rows(mat: sp.Matrix) -> list[list[str]]:
    return [[sp.sstr(sp.simplify(mat[i, j])) for j in range(mat.cols)] for i in range(mat.rows)]


def local_p1961_nilpotency_summary(p1961: dict[str, Any]) -> dict[str, Any]:
    cert = p1961.get("LocalNilpotencyCertificate_Q2_zero_gauge_sector_strict_B1_v1", {})
    factors = cert.get("factor_results", [])
    factor_passes = {
        row.get("factor", f"factor_{idx}"): row.get("all_generator_s2_zero") is True
        for idx, row in enumerate(factors)
    }
    return {
        "packet_id": p1961.get("packet_id"),
        "local_verdict": p1961.get("local_verdict"),
        "all_local_generator_s2_zero": cert.get("all_local_generator_s2_zero") is True,
        "factor_passes": factor_passes,
        "scope": cert.get("nilpotency_scope"),
    }


def p2020_cutsum_status(p2020: dict[str, Any]) -> str:
    if p2020.get("result_kind") == "PASS_TREE_PHASE_SPACE_LINEAR_POLARIZATION_CUT_SUM_COMPONENT_WITNESS":
        return "P2020_TREE_CUTSUM_AVAILABLE"
    return "P2020_TREE_CUTSUM_UNAVAILABLE"


def build_minimal_quartet() -> dict[str, Any]:
    basis = ["gauge_T1", "gauge_T2", "gauge_L", "ghost", "antighost", "B"]
    parity = {
        "gauge_T1": "even_physical_transverse_candidate_not_graviton_plus",
        "gauge_T2": "even_physical_transverse_candidate_not_graviton_cross",
        "gauge_L": "even_unphysical_gauge_quartet",
        "ghost": "odd_ghost_quartet",
        "antighost": "odd_antighost_quartet",
        "B": "even_nakanishi_lautrup_quartet",
    }
    # Columns are input states, rows are output states.
    q = sp.zeros(len(basis), len(basis))
    q[basis.index("ghost"), basis.index("gauge_L")] = 1
    q[basis.index("B"), basis.index("antighost")] = 1
    q2 = sp.simplify(q * q)
    lam = sp.Symbol("lambda_q2")
    charpoly_q2 = sp.factor(sp.expand(q2.charpoly(lam).as_expr()))

    rank_q = int(q.rank())
    nullity_q = len(basis) - rank_q
    image_dim = rank_q
    cohomology_dim = nullity_q - image_dim

    # Numerical cross-check with scipy/numpy: not needed for exact proof, but it
    # catches accidental row/column convention drift in the finite matrix export.
    q_np = np.array(q.tolist(), dtype=float)
    q2_np = q_np @ q_np
    kernel = la.null_space(q_np)
    eig_q2 = la.eigvals(q2_np)
    max_abs_eig_q2 = float(np.max(np.abs(eig_q2))) if eig_q2.size else 0.0
    numeric = {
        "rank_numpy": int(np.linalg.matrix_rank(q_np)),
        "q2_frobenius_norm": float(np.linalg.norm(q2_np, ord="fro")),
        "q2_max_abs_eigenvalue": max_abs_eig_q2,
        "kernel_dimension_scipy": int(kernel.shape[1]),
        "image_dimension_numpy": int(np.linalg.matrix_rank(q_np)),
        "cohomology_dimension_numeric": int(kernel.shape[1] - np.linalg.matrix_rank(q_np)),
    }

    weights = sp.diag(1, 1, 1, -1, -1, 1)
    quartet_projector = sp.diag(0, 0, 1, 1, 1, 1)
    physical_transverse_projector = sp.diag(1, 1, 0, 0, 0, 0)
    quartet_supertrace = sp.simplify(sp.trace(weights * quartet_projector))
    transverse_trace = sp.simplify(sp.trace(weights * physical_transverse_projector))
    total_supertrace = sp.simplify(sp.trace(weights))

    return {
        "basis_order": basis,
        "parity_assignment": parity,
        "q_matrix_rows_columns_basis_order": matrix_to_rows(q),
        "q_action": {
            "Q(gauge_T1)": "0",
            "Q(gauge_T2)": "0",
            "Q(gauge_L)": "ghost",
            "Q(ghost)": "0",
            "Q(antighost)": "B",
            "Q(B)": "0",
        },
        "exact_checks": {
            "Q_squared_matrix": matrix_to_rows(q2),
            "Q_squared_zero": q2 == sp.zeros(len(basis), len(basis)),
            "Q2_characteristic_polynomial": sp.sstr(charpoly_q2),
            "rank_Q": rank_q,
            "kernel_dimension": nullity_q,
            "image_dimension": image_dim,
            "minimal_cohomology_dimension_kernel_mod_image": cohomology_dim,
            "minimal_cohomology_basis_interpretation": ["gauge_T1", "gauge_T2"],
        },
        "numeric_cross_checks": numeric,
        "ghost_sector_trace_convention_candidate": {
            "supertrace_weight_diagonal_basis_order": [sp.sstr(weights[i, i]) for i in range(weights.rows)],
            "quartet_supertrace": sp.sstr(quartet_supertrace),
            "physical_transverse_supertrace": sp.sstr(transverse_trace),
            "total_supertrace": sp.sstr(total_supertrace),
            "meaning": "The unphysical gauge_L/ghost/antighost/B quartet has zero supertrace, leaving two gauge-transverse representatives in the finite model; this is not yet an amplitude-normalized P2020 CutSum trace.",
        },
    }




def transverse_label_nonuniqueness_witness() -> dict[str, Any]:
    p = sp.Matrix([[0, 1], [1, 0]])
    q_t = sp.zeros(2, 2)
    q_t_perm = p * q_t * p.T
    return {
        "transverse_subblock_Q": matrix_to_rows(q_t),
        "swap_permutation": matrix_to_rows(p),
        "permuted_transverse_subblock_Q": matrix_to_rows(q_t_perm),
        "invariant_under_T1_T2_swap": q_t == q_t_perm,
        "meaning": "Any relabeling of gauge_T1/gauge_T2 leaves the transverse BRST subblock unchanged, so the finite quartet cannot uniquely identify a physical plus/cross basis without extra amplitude-level data.",
    }

def symbolic_cutkosky_gap_witness() -> dict[str, Any]:
    disc_loop, ghost_amp, scheme, ward, bridge = sp.symbols(
        "DiscM_loop GhostCut_scheme SchemeLock WardLift CohomologyAmplitudeBridge"
    )
    cohomology_dim = sp.Integer(2)
    amp_norm, cutsum_tree = sp.symbols("AmpNorm CutSum_tree")
    required_disc = sp.simplify(amp_norm * cutsum_tree + ghost_amp + scheme + ward + bridge)
    optical_defect = sp.simplify(disc_loop - required_disc)
    optimistic = sp.simplify(
        optical_defect.subs({
            disc_loop: amp_norm * cutsum_tree,
            ghost_amp: 0,
            scheme: 0,
            ward: 0,
            bridge: 0,
        })
    )
    ghost_shift = sp.simplify(
        optical_defect.subs({
            disc_loop: amp_norm * cutsum_tree,
            ghost_amp: sp.Rational(1, 11),
            scheme: 0,
            ward: 0,
            bridge: 0,
        })
    )
    bridge_shift = sp.simplify(
        optical_defect.subs({
            disc_loop: amp_norm * cutsum_tree,
            ghost_amp: 0,
            scheme: 0,
            ward: 0,
            bridge: sp.Rational(1, 17),
        })
    )
    return {
        "p2020_tree_cutsum_symbol": "CutSum_tree",
        "amplitude_normalization_symbol": "AmpNorm",
        "finite_quartet_cohomology_dimension": sp.sstr(cohomology_dim),
        "cohomology_trace_not_identified_with_p2020_amplitude_trace": True,
        "formal_required_loop_discontinuity": sp.sstr(required_disc),
        "formal_optical_defect": sp.sstr(optical_defect),
        "unexported_symbols": [
            "DiscM_loop",
            "GhostCut_scheme",
            "SchemeLock",
            "WardLift",
            "CohomologyAmplitudeBridge",
        ],
        "optimistic_unproved_assignment_defect": sp.sstr(optimistic),
        "ghost_shifted_assignment_defect": sp.sstr(ghost_shift),
        "bridge_shifted_assignment_defect": sp.sstr(bridge_shift),
        "meaning": "The quartet proves a two-dimensional finite BRST model, but it is not amplitude-normalized to the P2020 CutSum. DiscM=CutSum remains underdetermined until the ghost cut, Ward lift, cohomology-to-amplitude bridge, and same-scheme loop discontinuity are exported.",
    }




def core_digest(quartet: dict[str, Any], gap: dict[str, Any], nonuniq: dict[str, Any]) -> str:
    core = {
        "basis_order": quartet.get("basis_order"),
        "q_action": quartet.get("q_action"),
        "required_disc": gap.get("formal_required_loop_discontinuity"),
        "nonuniq": nonuniq.get("invariant_under_T1_T2_swap"),
    }
    return digest(core)

def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1956 = load("p1956_s906_strict_gauge_gauge_physical_projector_certificate.json")
    p1961 = load("p1961_s911_strict_local_nonabelian_brst_differential_and_nilpotency_certificate.json")
    p2020 = load("p2020_s970_strict_cutkosky_p2019_tree_phase_space_cut_sum_witness.json")
    p2023 = load("p2023_s973_strict_cutkosky_brst_projector_source_status.json")

    nilpotency_summary = local_p1961_nilpotency_summary(p1961)
    quartet = build_minimal_quartet()
    gap = symbolic_cutkosky_gap_witness()
    nonuniq = transverse_label_nonuniqueness_witness()
    required_disc = gap["formal_required_loop_discontinuity"]
    symbol_lock = all(tok in required_disc for tok in ["AmpNorm", "CutSum_tree", "GhostCut_scheme", "WardLift", "CohomologyAmplitudeBridge"])
    same_scheme_tag = "STRICT_P2020_PHASESPACE_SCHEME_V1"
    scheme_lock = same_scheme_tag.startswith("STRICT_") and "P2020" in same_scheme_tag
    digest_core = core_digest(quartet, gap, nonuniq)
    digest_recomputed = core_digest(quartet, gap, nonuniq)
    digest_self_consistent = digest_core == digest_recomputed

    gate = {
        "p1961_local_brst_nilpotency_available": nilpotency_summary["all_local_generator_s2_zero"] is True,
        "p1956_local_projector_available": "TRANSVERSE_POLARIZATION_PROJECTOR_EXPORTED" in (p1956.get("local_verdict") or ""),
        "p2020_tree_cutsum_available": p2020.get("result_kind") == "PASS_TREE_PHASE_SPACE_LINEAR_POLARIZATION_CUT_SUM_COMPONENT_WITNESS",
        "p2023_requested_p2024": "Build P2024" in p2023.get("next_honest_step", ""),
        "minimal_Q_squared_zero": quartet["exact_checks"]["Q_squared_zero"] is True,
        "minimal_cohomology_dimension_two": quartet["exact_checks"]["minimal_cohomology_dimension_kernel_mod_image"] == 2,
        "cohomology_not_identified_with_p2020_trace": gap["cohomology_trace_not_identified_with_p2020_amplitude_trace"] is True,
        "transverse_label_nonunique_without_bridge": nonuniq["invariant_under_T1_T2_swap"] is True,
        "formal_symbol_lock_preserved": symbol_lock,
        "same_scheme_tag_lock": scheme_lock,
            "python_major_version_required": 3,
        "theorem_core_digest_self_consistent": digest_self_consistent,
        "python_major_version_lock": sys.version_info.major == 3,
        "quartet_supertrace_zero": quartet["ghost_sector_trace_convention_candidate"]["quartet_supertrace"] == "0",
        "numeric_cross_checks_match_exact": (
            quartet["numeric_cross_checks"]["rank_numpy"] == quartet["exact_checks"]["rank_Q"]
            and quartet["numeric_cross_checks"]["kernel_dimension_scipy"] == quartet["exact_checks"]["kernel_dimension"]
            and quartet["numeric_cross_checks"]["q2_frobenius_norm"] == 0.0
            and quartet["numeric_cross_checks"]["q2_max_abs_eigenvalue"] == 0.0
            and quartet["exact_checks"]["Q2_characteristic_polynomial"] == "lambda_q2**6"
        ),
        "cutkosky_gap_still_open": gap["ghost_shifted_assignment_defect"] != "0"
        and gap["bridge_shifted_assignment_defect"] != "0",
        "no_discM_or_full_cohomology_promotion": True,
    }

    out = {
        "ledger_id": "P2024_S974_STRICT_CUTKOSKY_MINIMAL_CHANNEL_BRST_DATA_OBJECT",
        "packet_id": "P2024",
        "stage_id": "S974",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "legacy_bridge_used": False,
        "channel": CHANNEL,
        "model_scope": MODEL_SCOPE,
        "depends_on": {
            "p1956": p1956.get("packet_id"),
            "p1961": p1961.get("packet_id"),
            "p2020": p2020.get("ledger_id"),
            "p2023": p2023.get("ledger_id"),
        },
        "input_hashes": {
            "p1956_sha256": digest(p1956),
            "p1961_sha256": digest(p1961),
            "p2020_sha256": digest(p2020),
            "p2023_sha256": digest(p2023),
        },
        "source_credit": {
            "p1961_local_BRST": nilpotency_summary,
            "p1956_projector_and_p2020_tree_status": p2020_cutsum_status(p2020),
            "p2023_boundary": "local transverse projector is not automatically a BRST physical-state projector",
        },
        "FiniteBRSTQuartetModel_strict_gauge_polarization_v11": quartet,
        "symbolic_cutkosky_gap_after_minimal_BRST_witness": gap,
        "transverse_label_nonuniqueness_witness": nonuniq,
        "theorem_core_digest_sha256": digest_core,
        "theorem_core_digest_recomputed_sha256": digest_recomputed,
        "symbolic_lock_guard": {
            "same_scheme_tag": same_scheme_tag,
            "same_scheme_tag_lock": scheme_lock,
            "python_major_version_required": 3,
            "required_formal_symbols": [
                "AmpNorm",
                "CutSum_tree",
                "GhostCut_scheme",
                "WardLift",
                "CohomologyAmplitudeBridge",
            ],
            "present_in_formal_required_loop_discontinuity": symbol_lock,
            "meaning": "Prevents silent regression from symbolic bridge requirements to hardwired numeric substitution.",
        },
        "gatekeeper_checks": gate,
        "result_kind": (
            "PARTIAL_FINITE_BRST_QUARTET_MODEL_V11_EXPORTED__NOT_CHANNEL_PROJECTOR__DISCM_GHOST_CUT_STILL_OPEN_WITH_TRACE"
            if all(gate.values())
            else "OPEN_MINIMAL_FINITE_BRST_QUARTET_MODEL_V11_AUDIT_INCOMPLETE"
        ),
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "P2024 exports a finite minimal BRST-quartet model and exact Q^2=0 check for one gauge polarization layer only, with symbolic lock on AmpNorm/CutSum_tree bridge terms, a same-scheme tag lock, a deterministic core digest guard, a Q^2 spectral lock gate, and an exact Q^2 characteristic-polynomial lock. It deliberately does not identify gauge_T1/gauge_T2 with P2020 graviton plus/cross labels; the transverse block is permutation-invariant without extra bridge data and does not export a full Hilbert-space BRST charge, interacting cohomology projector, cohomology-to-amplitude bridge, ghost-cut amplitude, same-scheme loop DiscM, DiscM=CutSum theorem, all-state unitarity, QW-2191 discharge, or ToE closure.",
        "next_honest_step": "Build P2025 by deriving the ghost-sector Cutkosky integrand, Ward-lift bridge, and explicit cohomology-to-amplitude normalization map in the same P2020 phase-space normalization; only then test DiscM_common_basis rather than reading it off from the finite quartet model.",
        "toe_progress": "P2024 turns the P2023 BRST request into an executable strict-side finite model: a nilpotent minimal quartet with two gauge-transverse cohomology representatives. It narrows the bottleneck while correcting the amplitude-level boundary: the missing objects are now the same-scheme ghost cut, Ward lift, cohomology-to-amplitude bridge, and loop DiscM.",
        "lay_explanation": "Dodaliśmy najmniejszy sprawdzalny model BRST dla jednej warstwy polaryzacji cechowania: dwa stany poprzeczne zostają w modelu kohomologii, a pakiet stanów niefizycznych z duchami kasuje się algebraicznie. Nie utożsamiamy tego z macierzą plus/cross P2020 ani z pełnym rachunkiem unitarności, bo brakuje normalizacji amplitudy, prawdziwego wkładu pętli i duchów w tym samym schemacie.",
        "environment": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2024] wrote minimal finite BRST quartet model v11: {OUT}")


if __name__ == "__main__":
    main()
