#!/usr/bin/env python3
"""P1954 S904 strict dressed amplitude nonavailability theorem.

This packet follows P1953.  It does not try to invent a dressed
graviton->gauge_gauge amplitude from sector-level L_total scaffolds.  It
proves the narrower repository-state statement that the current strict exports
under-determine that amplitude: the required 4D hAA vertex, BRST physical
projector, dressed propagator residues, and same-scheme DiscM/CutSum values are
not exported.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1954_s904_strict_dressed_amplitude_nonavailability_theorem.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def has_key_deep(obj: Any, key_name: str) -> bool:
    if isinstance(obj, dict):
        return any(k == key_name or has_key_deep(v, key_name) for k, v in obj.items())
    if isinstance(obj, list):
        return any(has_key_deep(v, key_name) for v in obj)
    return False


def text_contains(obj: Any, needle: str) -> bool:
    return needle.lower() in json.dumps(obj, sort_keys=True, ensure_ascii=False).lower()


def req_row(
    req_id: str,
    required_object: str,
    observed_sources: list[str],
    verdict: str,
    exact_missing_data: list[str],
) -> dict[str, object]:
    return {
        "req_id": req_id,
        "required_object": required_object,
        "observed_sources": observed_sources,
        "verdict": verdict,
        "exact_missing_data": exact_missing_data,
    }


def main() -> None:
    p1852 = load("p1852_s802_strict_b1_brst_anomaly_and_cutkosky_seed_witness_checkpoint.json")
    p1862 = load("p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json")
    p1866 = load("p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json")
    p1906 = load("p1906_s856_strict_c1_gr_diagram_inventory_stub_probe.json")
    p1907 = load("p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json")
    p1908 = load("p1908_s858_strict_c1_gr_divergent_coefficients_table_v1_probe.json")
    p1909 = load("p1909_s859_strict_c1_gr_counterterm_st_cutkosky_binding_probe.json")
    p1910 = load("p1910_s860_strict_c1_gr_coefficient_table_v1_probe.json")
    p1913 = load("p1913_s863_strict_c1_gr_common_basis_unitarity_background_probe.json")
    p1953 = load("p1953_s903_strict_dressed_cutkosky_amplitude_availability_audit.json")

    interface = p1953.get("dressed_amplitude_interface_contract") or {}
    required_fields = interface.get("required_fields") or []

    ltotal_has_sector_registry = "full_lagrangian_term_registry_non_skeleton" in p1907
    ltotal_has_1d_proxy = "eom_phi_proxy_1d" in (p1866.get("eom_export") or {})
    has_hAA_vertex = (
        has_key_deep(p1866, "hAA_vertex_tensor")
        or has_key_deep(p1907, "hAA_vertex_tensor")
        or text_contains(p1907, "graviton_gauge_gauge_vertex")
        or text_contains(p1866, "graviton_gauge_gauge_vertex")
    )
    has_brt_projector = (
        has_key_deep(p1852, "brst_physical_state_projector")
        or text_contains(p1852, "physical-state projector")
    )
    p1862_note = ((p1862.get("dressed_pole_residue_seed_table") or {}).get("note") or "")
    has_computed_dressed_residues = (
        "seed" not in p1862_note.lower()
        and has_key_deep(p1862, "residue_values_per_pole")
    )
    p1913_scheme = ((p1913.get("common_basis_definition") or {}).get("scheme") or "MISSING")
    scheme_locked = p1913_scheme == "MSbar_B1_seed"

    grmix_rows = [
        row for row in (p1913.get("unitarity_common_basis_rows_v1") or [])
        if row.get("channel") == "grmix"
    ]
    grmix_evaluated = bool(
        grmix_rows
        and all(row.get("status") in {"PASS_ZERO", "PASS_SYMBOLIC_CHANNEL_GLOBAL"} for row in grmix_rows)
    )

    inventory_has_gravity_mixed = text_contains(p1906, "d_gravity_mixed_counterterm_support")
    coeffs_are_placeholders = text_contains(p1910, "A_gr/epsilon") or text_contains(p1908, "OPEN_NO_NUMERIC")
    binding_is_symbolic_only = text_contains(p1909, "OPEN_SYMBOLIC_ONLY")

    missing_matrix = [
        req_row(
            "M1",
            "4D graviton-gauge-gauge Feynman vertex V_hAA^{mu nu;rho sigma;ab}(p,q)",
            [
                f"P1907 sector registry present={ltotal_has_sector_registry}",
                f"P1866 1D proxy present={ltotal_has_1d_proxy}",
                f"hAA vertex exported={has_hAA_vertex}",
            ],
            "MISSING",
            [
                "metric perturbation convention g=eta+kappa*h",
                "delta^3 S / delta h delta A delta A tensor",
                "momentum routing and gauge-index normalization",
                "contact terms from R^2/Ricci^2/Riemann^2 sector if included",
            ],
        ),
        req_row(
            "M2",
            "BRST physical-state projector for the gauge_gauge final state",
            [
                f"P1852 seed contract present={'cutkosky_seed_contract' in p1852}",
                f"BRST physical projector exported={has_brt_projector}",
            ],
            "MISSING",
            [
                "physical polarization sum",
                "ghost subtraction/cancellation trace",
                "proof that the cut sum is restricted to BRST cohomology",
            ],
        ),
        req_row(
            "M3",
            "computed dressed graviton propagator pole list and residues",
            [
                f"P1862 note={p1862_note}",
                f"computed residues exported={has_computed_dressed_residues}",
            ],
            "SEED_ONLY",
            [
                "inverse quadratic operator in fixed gauge",
                "pole polynomial after one-loop dressing",
                "residue derivative at each physical pole",
                "ghost-pole exclusion theorem",
            ],
        ),
        req_row(
            "M4",
            "same-scheme DiscM_grmix and CutSum_grmix coefficients",
            [
                f"P1913 scheme={p1913_scheme}",
                f"grmix row evaluated={grmix_evaluated}",
                f"scheme_locked_to_MSbar_B1_seed={scheme_locked}",
            ],
            "OPEN_PLACEHOLDER",
            [
                "alpha_gr, beta_gr, gamma_gr values",
                "alpha_gr_cut, beta_gr_cut, gamma_gr_cut values",
                "DiscM_minus_CutSum simplified expression",
                "scheme tag MSbar_B1_seed across both sides",
            ],
        ),
        req_row(
            "M5",
            "diagram coefficient and integral reductions for gravity/gauge channel",
            [
                f"gravity mixed diagram inventory present={inventory_has_gravity_mixed}",
                f"coefficients are placeholders={coeffs_are_placeholders}",
                f"binding symbolic only={binding_is_symbolic_only}",
            ],
            "OPEN_PLACEHOLDER",
            [
                "master integral reductions for d_gravity_mixed_counterterm_support",
                "finite and pole parts in the same basis",
                "channel normalization linking coefficient table to amplitude",
            ],
        ),
    ]

    all_acceptance_fields_absent_or_failed = not (
        has_hAA_vertex
        and has_brt_projector
        and has_computed_dressed_residues
        and grmix_evaluated
        and scheme_locked
        and not coeffs_are_placeholders
    )

    underdetermination_witness = {
        "free_symbols_if_one_attempts_derivation_now": [
            "V_hAA_tensor_free",
            "Pi_hh_dressed_free",
            "P_BRST_phys_free",
            "alpha_gr_minus_alpha_gr_cut",
            "beta_gr_minus_beta_gr_cut",
            "gamma_gr_minus_gamma_gr_cut",
        ],
        "meaning": (
            "The current sector-level L_total and placeholder common-basis rows admit multiple "
            "incompatible dressed amplitude completions while preserving all currently exported "
            "seed/proxy artifacts. Therefore a unique theorem-grade Cutkosky amplitude is not "
            "derivable from current exports."
        ),
    }

    out = {
        "packet_id": "P1954",
        "stage_id": "S904",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "local_verdict": "FORMAL_NONAVAILABILITY_THEOREM_EXPORTED__DRESSED_AMPLITUDE_UNDERDETERMINED",
        "route": "strict_only",
        "legacy_bridge_used": False,
        "depends_on": {
            "p1852_present": "cutkosky_seed_contract" in p1852,
            "p1862_present": "dressed_pole_residue_seed_table" in p1862,
            "p1866_present": "full_lagrangian_non_skeleton" in p1866,
            "p1907_present": "full_lagrangian_term_registry_non_skeleton" in p1907,
            "p1913_present": "unitarity_common_basis_rows_v1" in p1913,
            "p1953_present": "dressed_amplitude_interface_contract" in p1953,
        },
        "input_hashes": {
            "p1852_sha256": digest(p1852),
            "p1862_sha256": digest(p1862),
            "p1866_sha256": digest(p1866),
            "p1907_sha256": digest(p1907),
            "p1908_sha256": digest(p1908),
            "p1909_sha256": digest(p1909),
            "p1910_sha256": digest(p1910),
            "p1913_sha256": digest(p1913),
            "p1953_sha256": digest(p1953),
        },
        "p1953_required_interface_fields": required_fields,
        "minimal_missing_data_matrix": missing_matrix,
        "acceptance_condition_currently_failed": bool(all_acceptance_fields_absent_or_failed),
        "formal_nonavailability_theorem": {
            "statement": (
                "On the current repository state, the strict dressed "
                "graviton->gauge_gauge Cutkosky amplitude object is underdetermined and "
                "not derivable as a theorem-grade export."
            ),
            "proof_trace": [
                "P1907 exports a sector-level non-skeleton registry, but not a 4D hAA Feynman vertex tensor.",
                "P1866 exports a 1D symbolic proxy chain, not a 4D perturbative expansion around g=eta+kappa*h.",
                "P1852 exports BRST/Cutkosky seed contracts, not a channel BRST physical-state projector.",
                "P1862 explicitly marks dressed residues as seed-inherited and still requiring full propagator computation.",
                "P1913 leaves grmix DiscM/CutSum in open symbolic placeholders under MSbar_candidate rather than MSbar_B1_seed.",
                "Therefore the P1953 acceptance contract cannot be satisfied from current strict exports.",
            ],
            "underdetermination_witness": underdetermination_witness,
        },
        "safe_consequence": {
            "may_continue_using": [
                "P1950 declared B1 counterterm cancellation",
                "P1951 seed phase-space positivity",
                "P1952 QW-2049 seed-local positivity rectangle",
            ],
            "must_not_claim": [
                "full dressed graviton->gauge_gauge Cutkosky equality",
                "BRST-projected optical theorem closure",
                "ghost-free dressed propagator theorem",
                "global UR_link theorem",
            ],
        },
        "next_solver_input_contract": {
            "recommended_packet": "P1955",
            "minimum_new_exports_required": [
                "V_hAA_tensor_strict_B1_v1",
                "BRSTPhysicalProjector_gauge_gauge_strict_B1_v1",
                "DressedGravitonPropagatorPoleResidueTable_strict_B1_v1",
                "DiscM_CutSum_grmix_MSbar_B1_seed_common_basis_v1",
                "SchemeLock_MSbar_B1_seed_for_P1913_grmix_v1",
            ],
        },
        "false_pass_guard": "This is a nonavailability theorem for the current export state, not a no-go theorem about the theory in principle.",
        "higher_reasoning_required_for_next_step": True,
        "next_honest_step": "Build P1955 with high reasoning: derive V_hAA_tensor_strict_B1_v1 from a 4D metric perturbation expansion of L_total, or prove that the current L_total export lacks the metric-density detail needed for that derivation.",
        "lay_explanation": "Nie da sie uczciwie policzyc pelnej amplitudy z obecnych danych, bo brakuje reguly wierzcholka, projektora stanow fizycznych i policzonych residuow propagatora. To nie obala teorii; mowi dokladnie, jaki kawalek matematycznej maszynerii trzeba teraz zbudowac.",
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
