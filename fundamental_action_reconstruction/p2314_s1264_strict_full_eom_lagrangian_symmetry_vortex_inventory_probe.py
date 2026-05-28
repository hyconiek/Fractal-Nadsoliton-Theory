#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

IN_ALPHA = GEN / "alpha_geo_strict_derived_v1.json"
IN_2030 = GEN / "p2030_s980_strict_tensor_resolved_projection_source_audit.json"
IN_2031 = GEN / "p2031_s981_strict_b1_tensor_component_profile_table_scaffold.json"
IN_2032 = GEN / "p2032_s982_strict_b1_metric_gauge_component_projection_rule_audit.json"
IN_2033 = GEN / "p2033_s983_strict_curved_b1_metric_ansatz_nonavailability_theorem.json"
IN_2300 = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json"
IN_2313 = GEN / "p2313_s1263_strict_provider_margin_route_stop_and_source_class_contract_probe.json"
IN_RELEASE7 = REPO / "RELEASE_7_STRICT_PROJECTIVE_OPERATIONAL_MODEL_BRIEF.md"
IN_LEGACY_VORTEX = REPO / "28 FULL LEPTON MASS HIERARCHY VIA NONLINEAR HYDRODYNAMIC COUPLING.py"
OUT = GEN / "p2314_s1264_strict_full_eom_lagrangian_symmetry_vortex_inventory_probe.json"
MD = GEN / "p2314_s1264_strict_full_eom_lagrangian_symmetry_vortex_inventory_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def count_terms(text: str, terms: list[str]) -> dict[str, int]:
    return {term: len(re.findall(re.escape(term), text, flags=re.IGNORECASE)) for term in terms}


def main() -> None:
    GEN.mkdir(exist_ok=True)
    alpha = load(IN_ALPHA)
    p2030 = load(IN_2030)
    p2031 = load(IN_2031)
    p2032 = load(IN_2032)
    p2033 = load(IN_2033)
    p2300 = load(IN_2300)
    p2313 = load(IN_2313)
    release7_text = read_text(IN_RELEASE7)
    legacy_vortex_text = read_text(IN_LEGACY_VORTEX)

    schematic_lagrangian_present = "\\mathcal{L}_{ZTP}" in release7_text and "schematic representative form" in release7_text
    release7_candidate_layer_present = "strict core candidate Lagrangian / EOM layer" in release7_text
    release7_not_finished = "not yet" in release7_text and "finished particle-level emergence theorem" in release7_text

    p2030_checks = p2030.get("gatekeeper_checks", {}) or {}
    p2031_checks = p2031.get("gatekeeper_checks", {}) or {}
    p2032_checks = p2032.get("gatekeeper_checks", {}) or {}
    p2033_checks = p2033.get("gatekeeper_checks", {}) or {}
    p2300_checks = p2300.get("gatekeeper_checks", {}) or {}
    p2313_checks = p2313.get("gatekeeper_checks", {}) or {}

    eom_lagrangian_inventory = {
        "schematic_candidate_lagrangian_present": schematic_lagrangian_present,
        "candidate_lagrangian_eom_layer_present_in_release7": release7_candidate_layer_present,
        "release7_explicitly_not_finished_particle_emergence_theorem": release7_not_finished,
        "partial_adm_bianchi_lapse_chain_present": bool(p2030_checks.get("adm_lapse_chain_present") is True),
        "partial_spatial_eom_provider_matrix_present": bool(p2300_checks.get("provider_matrix_consistent") is True and p2300_checks.get("canonical_solution_exported") is True),
        "full_tensor_component_profile_table_missing": bool(p2030_checks.get("tensor_component_profile_table_missing") is True and p2031_checks.get("all_entries_marked_missing") is True),
        "curved_metric_projection_rule_missing": bool(p2032_checks.get("rule_ready") is False and p2033_checks.get("minimal_curved_b1_ansatz_not_exported") is True),
        "full_strict_eom_exported": False,
        "full_theorem_grade_lagrangian_for_task3_exported": False,
        "answer": "No: the repo has a schematic/candidate strict core Lagrangian/EOM layer and several partial ADM/Bianchi-I EOM witnesses, but not a full tensor-resolved EOM or theorem-grade full Lagrangian sufficient for Task-3 G1/G3 closure.",
    }

    n = 16
    p = 1.0 / n
    entropy_nats = -n * p * math.log(p)
    alpha_value_ok = alpha.get("value") == "4 ln(2)" and abs(entropy_nats - 4.0 * math.log(2.0)) < 1e-15
    sample_permutation = [3, 0, 15, 1, 14, 2, 13, 4, 12, 5, 11, 6, 10, 7, 9, 8]
    uniform = [p for _ in range(n)]
    permuted = [uniform[i] for i in sample_permutation]
    entropy_permuted = -sum(x * math.log(x) for x in permuted)
    alpha_symmetry_audit = {
        "alpha_geo_strict_source_loaded": alpha.get("status") == "actual_exported_strict_derived_source_upgrade_value",
        "definition": alpha.get("definition", {}),
        "computed_entropy_nats": entropy_nats,
        "computed_entropy_equals_4_ln2": alpha_value_ok,
        "sample_permutation_entropy": entropy_permuted,
        "permutation_invariant_for_uniform_measure": abs(entropy_permuted - entropy_nats) < 1e-15,
        "breaks_symmetry_by_itself": False,
        "can_weight_candidate_branch": True,
        "strict_selector_source_exported": False,
        "answer": "No: 4 ln 2 is a strict-side scalar entropy/source value.  For the uniform 16-state witness it is permutation-invariant, so by itself it does not choose an orientation, branch, chirality, or selector.  It can weight a future candidate only if an additional typed selector/orientation theorem is exported.",
    }

    vortex_terms = ["vortex", "twist", "torsion", "K_torsion", "beta_tors", "skręc", "wir"]
    legacy_vortex_counts = count_terms(legacy_vortex_text, vortex_terms)
    strict_generated_vortex_hits = []
    for path in GEN.glob("*.json"):
        if path == OUT:
            continue
        text = read_text(path)
        if re.search(r"vortex|twist|torsion|skręc|\bwir\b", text, re.IGNORECASE):
            strict_generated_vortex_hits.append(path.name)
    vortex_torsion_audit = {
        "legacy_or_numeric_vortex_terms_detected": any(count > 0 for count in legacy_vortex_counts.values()),
        "legacy_or_numeric_vortex_term_counts": legacy_vortex_counts,
        "generated_artifacts_with_vortex_or_torsion_terms": strict_generated_vortex_hits[:40],
        "strict_twisted_vortex_source_exported": False,
        "strict_torsion_to_task3_provider_margin_bridge_exported": False,
        "legacy_beta_tors_role_transfer_allowed": False,
        "answer": "Not as a strict bridge: the repository contains legacy/numerical vortex/torsion material, but P2314 finds no current strict exported twisted-vortex/torsion source that supplies the P2281/P2302 margin orientation.  Any relation to a twisted vortex remains a future source-class hypothesis, not a current strict result.",
    }

    computational_answers = {
        "is_full_eom_exported": eom_lagrangian_inventory["full_strict_eom_exported"],
        "is_full_lagrangian_exported": eom_lagrangian_inventory["full_theorem_grade_lagrangian_for_task3_exported"],
        "does_4ln2_break_symmetry": alpha_symmetry_audit["breaks_symmetry_by_itself"],
        "does_current_strict_route_have_twisted_vortex_bridge": vortex_torsion_audit["strict_twisted_vortex_source_exported"],
        "current_provider_margin_route_stopped": bool(p2313_checks.get("all_current_route_blockers_active") is True),
    }

    theorem_export = {
        "statement_id": "P2314_FULL_EOM_LAGRANGIAN_SYMMETRY_VORTEX_INVENTORY_THEOREM",
        "formal_statement": (
            "In the current strict export set, a schematic candidate nadsoliton Lagrangian/EOM layer and partial ADM/Bianchi-I lapse/spatial-EOM witnesses exist, "
            "but a full tensor-resolved EOM and theorem-grade full Lagrangian sufficient for Task-3 G1/G3 are not exported.  The strict value 4 ln 2 is "
            "a permutation-invariant Shannon entropy scalar for the uniform 16-state witness and does not by itself break selector/orientation symmetry.  Legacy/numerical vortex/torsion material is present, but no current strict twisted-vortex/torsion bridge to the P2281/P2302 margin is exported."
        ),
        "proof_bits": {
            "eom_lagrangian_inventory": eom_lagrangian_inventory,
            "alpha_symmetry_audit": alpha_symmetry_audit,
            "vortex_torsion_audit": vortex_torsion_audit,
            "computational_answers": computational_answers,
        },
        "not_claimed": [
            "full strict EOM",
            "full theorem-grade Lagrangian for Task-3",
            "4 ln 2 selector/symmetry-breaking theorem",
            "strict twisted-vortex/torsion bridge",
            "strict G1 closure",
            "strict G3 closure",
            "QW-2191 discharge",
            "legacy-kernel role transfer",
            "ToE closure",
        ],
    }
    theorem_fingerprint = sha256_json(theorem_export)

    gatekeeper_checks = {
        "schematic_candidate_lagrangian_present": schematic_lagrangian_present,
        "full_tensor_eom_not_exported": eom_lagrangian_inventory["full_strict_eom_exported"] is False,
        "full_task3_lagrangian_not_exported": eom_lagrangian_inventory["full_theorem_grade_lagrangian_for_task3_exported"] is False,
        "p2030_p2033_tensor_gaps_loaded": eom_lagrangian_inventory["full_tensor_component_profile_table_missing"] and eom_lagrangian_inventory["curved_metric_projection_rule_missing"],
        "alpha_geo_strict_loaded": alpha_symmetry_audit["alpha_geo_strict_source_loaded"],
        "alpha_entropy_computation_matches_4ln2": alpha_symmetry_audit["computed_entropy_equals_4_ln2"],
        "alpha_permutation_invariant_not_selector": alpha_symmetry_audit["permutation_invariant_for_uniform_measure"] and not alpha_symmetry_audit["breaks_symmetry_by_itself"],
        "strict_twisted_vortex_bridge_not_exported": vortex_torsion_audit["strict_twisted_vortex_source_exported"] is False,
        "legacy_torsion_role_not_transferred": vortex_torsion_audit["legacy_beta_tors_role_transfer_allowed"] is False,
        "p2313_route_stop_loaded": computational_answers["current_provider_margin_route_stopped"],
        "strict_g1_g3_not_updated": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2314_s1264_v1",
        "packet_id": "P2314",
        "stage_id": "S1264",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_INVENTORY_FULL_EOM_LAGRANGIAN_SYMMETRY_VORTEX_NOT_CLOSED",
        "result_kind": "STRICT_COMPUTATIONAL_INVENTORY_ANSWERS_FULL_EOM_LAGRANGIAN_4LN2_VORTEX_QUESTIONS",
        "strict_full_eom_lagrangian_symmetry_vortex_inventory_probe": {
            "probe_id": "P2314_S1264_STRICT_FULL_EOM_LAGRANGIAN_SYMMETRY_VORTEX_INVENTORY",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p2030": "generated/p2030_s980_strict_tensor_resolved_projection_source_audit.json",
                "p2031": "generated/p2031_s981_strict_b1_tensor_component_profile_table_scaffold.json",
                "p2032": "generated/p2032_s982_strict_b1_metric_gauge_component_projection_rule_audit.json",
                "p2033": "generated/p2033_s983_strict_curved_b1_metric_ansatz_nonavailability_theorem.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
                "p2313": "generated/p2313_s1263_strict_provider_margin_route_stop_and_source_class_contract_probe.json",
                "release7": "RELEASE_7_STRICT_PROJECTIVE_OPERATIONAL_MODEL_BRIEF.md",
                "legacy_vortex_scan": "28 FULL LEPTON MASS HIERARCHY VIA NONLINEAR HYDRODYNAMIC COUPLING.py",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p2030_sha256": sha256_file(IN_2030),
                "p2031_sha256": sha256_file(IN_2031),
                "p2032_sha256": sha256_file(IN_2032),
                "p2033_sha256": sha256_file(IN_2033),
                "p2300_sha256": sha256_file(IN_2300),
                "p2313_sha256": sha256_file(IN_2313),
                "release7_sha256": sha256_file(IN_RELEASE7),
                "legacy_vortex_scan_sha256": sha256_file(IN_LEGACY_VORTEX),
            },
            "eom_lagrangian_inventory": eom_lagrangian_inventory,
            "alpha_4ln2_symmetry_audit": alpha_symmetry_audit,
            "vortex_torsion_audit": vortex_torsion_audit,
            "computational_answers": computational_answers,
            "task3_g1_g3_update": {
                "G1_reduction_certainty": "OPEN_UNCHANGED",
                "G2_nonlinear_trajectory_realism": "CLOSED_FROM_P2301_UNCHANGED",
                "G3_operational_policy_rule": "OPEN_UNCHANGED",
                "reason": "P2314 is an inventory/proof audit; it exports no full EOM, no full Task-3 Lagrangian, no 4ln2 selector theorem, and no strict twisted-vortex bridge.",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2315_candidate",
            "goal": "If continuing strictly, choose one new source class and prove it: either a full tensor-resolved EOM/Lagrangian completion, an internal 4ln2 selector/orientation theorem, or a strict twisted-vortex/torsion bridge; then replay P2302/P2281 and recompute P2282 rows.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_WITH_COMPUTATIONAL_INVENTORY_NO_FALSE_PASS",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2314/S1264 — strict full-EOM/Lagrangian, 4ln2 symmetry, and vortex inventory",
            "",
            f"- Status: `{payload['status']}`",
            f"- Full EOM exported: `{computational_answers['is_full_eom_exported']}`",
            f"- Full theorem-grade Lagrangian for Task-3 exported: `{computational_answers['is_full_lagrangian_exported']}`",
            f"- `4 ln 2` breaks symmetry by itself: `{computational_answers['does_4ln2_break_symmetry']}`",
            f"- Current strict twisted-vortex bridge exported: `{computational_answers['does_current_strict_route_have_twisted_vortex_bridge']}`",
            "- G1/G3 update: `OPEN_UNCHANGED`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Concrete answers",
            "1. Full EOM: not yet; only partial ADM/Bianchi-I and tensor-gap witnesses are exported.",
            "2. Full Lagrangian: not theorem-grade for Task-3; the repo has a schematic/candidate strict core Lagrangian layer.",
            "3. `4 ln 2`: strict-derived entropy source, but permutation-invariant and not a selector by itself.",
            "4. Twisted vortex/torsion: present as legacy/numerical material, not as a strict exported bridge to provider_lift_per_step.",
            "",
            "## Guardrail statement",
            "P2314 does not close G1/G3, does not discharge QW-2191, does not transfer legacy torsion/vortex roles, and does not claim ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
