#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2325_s1275_strict_signed_source_to_axis_susceptibility_bridge_audit_probe.json"
MD = GEN / "p2325_s1275_strict_signed_source_to_axis_susceptibility_bridge_audit_probe.md"

SOURCE_FILES = {
    "P2324_AXIS_BRANCH_SUSCEPTIBILITY": GEN / "p2324_s1274_strict_axis_branch_susceptibility_nonpromotion_probe.json",
    "P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION": GEN / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.json",
    "P2308_CURRENT_INTERFACE_NONEXISTENCE": GEN / "p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json",
    "P1967_SHANNON_AXIS_TO_DELTA_SEL": GEN / "p1967_s917_strict_shannon_axis_source_to_delta_sel_tensor_map.json",
    "P729_PAIR12_ORBIT_DIRECTION_BRIDGE_NONEXPORT": GEN / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json",
    "P731_W_BREAK_PAIR12_PROMOTION_NONEXPORT": GEN / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json",
    "P734_SCOPE_TOPOLOGY_DESCENT_NONEXPORT": GEN / "p734_current_strict_t188_declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json",
    "P735_LOCAL_SCALAR_BIND_DESCENT_NONEXPORT": GEN / "p735_current_strict_t189_local_source_side_scalar_bind_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json",
}

GREP_PATTERNS = (
    "delta_sel",
    "WBreak",
    "w_break",
    "orbit direction",
    "signed source",
    "selector theorem",
    "descent bridge",
    "nonexport",
    "provider_lift_per_step",
    "QW-2191",
)

REQUIRED_SIGNED_SOURCE_FIELDS = (
    "strict_internal_signed_mu_source",
    "theta_ref_or_branch_reference",
    "pairwise_sign_descent_to_axis_branch",
    "task3_margin_response_functional",
    "p2281_replay_semantics",
    "selector_free_admissibility_theorem",
)


def load_json(path: Path) -> dict[str, Any]:
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
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def corpus_hits() -> list[dict[str, Any]]:
    paths = sorted(
        set(SOURCE_FILES.values())
        | set(ROOT.glob("P7*_CURRENT_STRICT*SELECTOR*.md"))
        | set(ROOT.glob("N4*_CURRENT_FIRST_STRICT*SHANNON*.md"))
        | set(GEN.glob("p7*_current_strict*selector*json"))
        | set(GEN.glob("p7*_current_strict*orbit*json"))
        | set(GEN.glob("p1967*json"))
        | set(GEN.glob("p2324*json"))
    )
    self_paths = {OUT.resolve(), MD.resolve(), Path(__file__).resolve()}
    rows: list[dict[str, Any]] = []
    for path in paths:
        if not path.exists() or path.is_dir() or path.resolve() in self_paths:
            continue
        text = read_text(path)
        lowered = text.lower()
        count = sum(lowered.count(pattern.lower()) for pattern in GREP_PATTERNS)
        if count == 0:
            continue
        first_line = 0
        first_excerpt = ""
        for idx, line in enumerate(text.splitlines(), start=1):
            if any(pattern.lower() in line.lower() for pattern in GREP_PATTERNS):
                first_line = idx
                first_excerpt = line.strip()[:240]
                break
        rows.append({
            "path": path.relative_to(REPO).as_posix(),
            "pattern_hit_count": count,
            "first_hit_line": first_line,
            "first_hit_excerpt": first_excerpt,
        })
    rows.sort(key=lambda row: (-int(row["pattern_hit_count"]), row["path"]))
    return rows


def candidate_source_rows(artifacts: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    p1967 = artifacts["P1967_SHANNON_AXIS_TO_DELTA_SEL"]
    p729 = artifacts["P729_PAIR12_ORBIT_DIRECTION_BRIDGE_NONEXPORT"]
    p731 = artifacts["P731_W_BREAK_PAIR12_PROMOTION_NONEXPORT"]
    p734 = artifacts["P734_SCOPE_TOPOLOGY_DESCENT_NONEXPORT"]
    p735 = artifacts["P735_LOCAL_SCALAR_BIND_DESCENT_NONEXPORT"]
    return [
        {
            "candidate_id": "P1967_SHANNON_DELTA_SEL_AXIS_SOURCE",
            "status": p1967.get("status"),
            "positive_evidence": [
                "strict Shannon axis source is exported",
                "legacy_bridge_used is false" if p1967.get("legacy_bridge_used") is False else "legacy bridge status unclear",
                "delta_sel tensor map exists at packet level",
            ],
            "strict_internal_signed_mu_source": False,
            "theta_ref_or_branch_reference": False,
            "pairwise_sign_descent_to_axis_branch": False,
            "task3_margin_response_functional": False,
            "p2281_replay_semantics": False,
            "selector_free_admissibility_theorem": False,
            "blocking_reason": "P1967 provides an axis/delta_sel typed source but still leaves global selector closure open and does not export the signed mu/theta_ref response needed by P2324.",
        },
        {
            "candidate_id": "P729_PAIR12_ORBIT_DIRECTION_LOCALIZATION",
            "status": p729.get("status"),
            "positive_evidence": [
                "pair1/pair2 ambiguity localized as opposite orbit directions",
                f"target_exported={p729.get('t183_target_exported_on_current_repo_state')}",
            ],
            "strict_internal_signed_mu_source": False,
            "theta_ref_or_branch_reference": True,
            "pairwise_sign_descent_to_axis_branch": False,
            "task3_margin_response_functional": False,
            "p2281_replay_semantics": False,
            "selector_free_admissibility_theorem": False,
            "blocking_reason": "The pair12 direction split is localized but the bridge target is explicitly not exported, so it cannot supply P2324 signed perturbation semantics.",
        },
        {
            "candidate_id": "P731_W_BREAK_PAIR12_ANTISYMMETRIC_SCORES",
            "status": p731.get("status"),
            "positive_evidence": [
                "w_break separates pair12 orbit-direction branches",
                f"pair1_sign={p731.get('pair1_w_break_branch_score_sign')}",
                f"pair2_sign={p731.get('pair2_w_break_branch_score_sign')}",
                f"antisymmetric={p731.get('w_break_pair12_branch_scores_are_antisymmetric')}",
            ],
            "strict_internal_signed_mu_source": False,
            "theta_ref_or_branch_reference": True,
            "pairwise_sign_descent_to_axis_branch": False,
            "task3_margin_response_functional": False,
            "p2281_replay_semantics": False,
            "selector_free_admissibility_theorem": False,
            "blocking_reason": "P731 is the closest signed-looking witness, but its promotion bridge is explicitly non-exported; sign scores do not descend into a strict mu/theta_ref response theorem.",
        },
        {
            "candidate_id": "P734_SCOPE_TOPOLOGY_SELECTOR_DESCENT",
            "status": p734.get("status"),
            "positive_evidence": [
                "declared scope source topology selector theorem exists",
                f"quotient_class_only={p734.get('current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only')}",
            ],
            "strict_internal_signed_mu_source": False,
            "theta_ref_or_branch_reference": False,
            "pairwise_sign_descent_to_axis_branch": False,
            "task3_margin_response_functional": False,
            "p2281_replay_semantics": False,
            "selector_free_admissibility_theorem": False,
            "blocking_reason": "The theorem remains quotient-class only and P731 does not descend to it, so it is not a signed branch source for P2324.",
        },
        {
            "candidate_id": "P735_LOCAL_SOURCE_SIDE_SCALAR_BIND",
            "status": p735.get("status"),
            "positive_evidence": [
                "local source-side scalar bind audited",
                f"branch_blind={p735.get('current_local_source_side_scalar_bind_is_pair12_branch_blind')}",
            ],
            "strict_internal_signed_mu_source": False,
            "theta_ref_or_branch_reference": False,
            "pairwise_sign_descent_to_axis_branch": False,
            "task3_margin_response_functional": False,
            "p2281_replay_semantics": False,
            "selector_free_admissibility_theorem": False,
            "blocking_reason": "The local scalar bind factors through shared cos-phi data and is branch-blind, so it cannot create the signed P2324 tilt.",
        },
    ]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    p2324_probe = artifacts["P2324_AXIS_BRANCH_SUSCEPTIBILITY"].get("strict_axis_branch_susceptibility_nonpromotion_probe", {}) or {}
    p2318_probe = artifacts["P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION"].get("strict_selector_lane_to_task3_margin_bridge_obligation_probe", {}) or {}
    required_lift = float((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("required_lift_per_step", 0.0))
    missing_fields = list((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("missing_required_bridge_fields", []))
    mu_by_pair = (p2324_probe.get("susceptibility_certificate", {}) or {}).get("mu_needed_for_required_lift_spread_by_pair", {})

    candidates = candidate_source_rows(artifacts)
    predicate_matrix = []
    for candidate in candidates:
        failed = [field for field in REQUIRED_SIGNED_SOURCE_FIELDS if not candidate[field]]
        predicate_matrix.append({
            "candidate_id": candidate["candidate_id"],
            "passed_predicates": [field for field in REQUIRED_SIGNED_SOURCE_FIELDS if candidate[field]],
            "failed_predicates": failed,
            "admissible_for_p2324_signed_mu_bridge": len(failed) == 0,
            "blocking_reason": candidate["blocking_reason"],
        })

    admissible = [row["candidate_id"] for row in predicate_matrix if row["admissible_for_p2324_signed_mu_bridge"]]
    bridge_audit_certificate = {
        "candidate_count": len(candidates),
        "admissible_candidate_count": len(admissible),
        "admissible_candidates": admissible,
        "all_candidates_fail_at_least_one_required_predicate": all(row["failed_predicates"] for row in predicate_matrix),
        "closest_candidate": "P731_W_BREAK_PAIR12_ANTISYMMETRIC_SCORES",
        "closest_candidate_still_nonexported": True,
        "p2324_required_mu_by_pair_loaded": mu_by_pair,
        "required_lift_per_step": required_lift,
        "formal_status": "NO_CURRENT_EXPORTED_SIGNED_SOURCE_TO_P2324_MU_RESPONSE_BRIDGE",
        "lay_symptom": "The repo contains several signs that look like a possible compass needle, but every current candidate either stays quotient-level, branch-blind, nonexported, or lacks replay/response semantics.",
    }

    bridge_obligation_update = {
        "required_lift_per_step": required_lift,
        "p2318_missing_fields_loaded": missing_fields,
        "p2324_mu_values_loaded": mu_by_pair,
        "signed_source_to_mu_bridge_exported": False,
        "fields_still_unfilled_after_signed_source_audit": missing_fields,
        "g1_g3_update_allowed": False,
        "reason": "No current signed-looking source satisfies the joint predicates needed to turn P2324 susceptibility into a strict Task-3 provider_lift_per_step bridge.",
    }

    theorem_export = {
        "theorem_name": "P2325 current signed-source to P2324 susceptibility bridge nonexistence audit",
        "formal_statement": (
            "Within the current exported source class P1967/P729/P731/P734/P735, there is no admissible signed-source bridge to the P2324 perturbation parameter mu.  "
            "P731 supplies the closest antisymmetric sign-looking witness, and P2324 shows finite mu values would be sufficient numerically, but the promotion/descent targets are explicitly nonexported or branch-blind and no candidate exports a Task-3 margin response functional, P2281 replay semantics, or selector-free admissibility theorem.  "
            "Therefore the P2324 susceptibility trace remains non-promoted and G1/G3 remain open."
        ),
        "proof_bits": {
            "candidate_count": bridge_audit_certificate["candidate_count"],
            "admissible_candidate_count": bridge_audit_certificate["admissible_candidate_count"],
            "all_candidates_fail_at_least_one_required_predicate": bridge_audit_certificate["all_candidates_fail_at_least_one_required_predicate"],
            "closest_candidate": bridge_audit_certificate["closest_candidate"],
            "p2324_required_mu_by_pair_loaded": mu_by_pair,
            "p2318_missing_field_count": len(missing_fields),
        },
        "scope_limits": [
            "finite current-export source-class audit, not a universal future nonexistence theorem",
            "does not deny that a future strict signed source could be exported",
            "does not promote P731 sign scores into a branch selector",
            "does not discharge QW-2191, update G1/G3, or close ToE",
        ],
        "nonpromotion_rule": "Do not use P1967/P729/P731/P734/P735 as the P2324 signed perturbation source unless a new export supplies mu/theta_ref descent, response functional, replay semantics, and selector-free admissibility.",
    }
    theorem_fingerprint = sha256_json(theorem_export)

    probe = {
        "probe_id": "P2325_S1275_STRICT_SIGNED_SOURCE_TO_AXIS_SUSCEPTIBILITY_BRIDGE_AUDIT",
        "source_packets": {key: path.relative_to(REPO).as_posix() for key, path in SOURCE_FILES.items()},
        "source_hashes": {f"{key}_sha256": sha256_file(path) for key, path in SOURCE_FILES.items()},
        "far_signed_source_grep_audit": {
            "patterns": list(GREP_PATTERNS),
            "hit_count": len(corpus_hits()),
            "top_hits": corpus_hits()[:30],
        },
        "required_signed_source_fields": list(REQUIRED_SIGNED_SOURCE_FIELDS),
        "candidate_source_rows": candidates,
        "predicate_matrix": predicate_matrix,
        "bridge_audit_certificate": bridge_audit_certificate,
        "bridge_obligation_update": bridge_obligation_update,
        "existing_blocker_context": {
            "p2324_result_kind": artifacts["P2324_AXIS_BRANCH_SUSCEPTIBILITY"].get("result_kind"),
            "p2318_result_kind": artifacts["P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION"].get("result_kind"),
            "p2308_result_kind": artifacts["P2308_CURRENT_INTERFACE_NONEXISTENCE"].get("result_kind"),
        },
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": theorem_fingerprint,
    }

    gatekeeper_checks = {
        "far_signed_source_grep_hits_found": len(probe["far_signed_source_grep_audit"]["top_hits"]) >= 5,
        "p2324_loaded": artifacts["P2324_AXIS_BRANCH_SUSCEPTIBILITY"].get("packet_id") == "P2324",
        "p2318_loaded": artifacts["P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION"].get("packet_id") == "P2318",
        "p1967_loaded": artifacts["P1967_SHANNON_AXIS_TO_DELTA_SEL"].get("packet_id") == "P1967",
        "candidate_class_nonempty": len(candidates) == 5,
        "p2324_mu_values_loaded": len(mu_by_pair) == 5,
        "no_admissible_signed_source_bridge": len(admissible) == 0,
        "all_candidates_fail_at_least_one_required_predicate": bridge_audit_certificate["all_candidates_fail_at_least_one_required_predicate"],
        "p731_closest_candidate_not_promoted": bridge_audit_certificate["closest_candidate_still_nonexported"],
        "p2318_bridge_fields_still_unfilled": len(bridge_obligation_update["fields_still_unfilled_after_signed_source_audit"]) == 6,
        "strict_g1_g3_not_updated": not bridge_obligation_update["g1_g3_update_allowed"],
        "no_selector_closure_claimed": True,
        "no_qw2191_discharge_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2325_s1275_v1",
        "packet_id": "P2325",
        "stage_id": "S1275",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_NO_CURRENT_SIGNED_SOURCE_TO_P2324_MU_BRIDGE",
        "result_kind": "STRICT_SIGNED_SOURCE_TO_AXIS_SUSCEPTIBILITY_BRIDGE_AUDIT_NO_G1_G3_UPDATE",
        "strict_signed_source_to_axis_susceptibility_bridge_audit_probe": probe,
        "recommended_next_honest_step": "Either export a real strict signed source with mu/theta_ref descent plus response/replay semantics, or prove a broader nonexistence theorem over larger signed-source classes.",
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_SIGNED_SOURCE_BRIDGE_NONEXISTENCE_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2325/S1275 — signed-source to P2324 susceptibility bridge audit\n\n"
        f"- Result kind: `{payload['result_kind']}`.\n"
        f"- Candidate source count: `{bridge_audit_certificate['candidate_count']}`.\n"
        f"- Admissible candidate count: `{bridge_audit_certificate['admissible_candidate_count']}`.\n"
        f"- Closest candidate: `{bridge_audit_certificate['closest_candidate']}` (still nonexported).\n"
        f"- P2324 mu values loaded: `{mu_by_pair}`.\n"
        "- Lay symptom: current exports contain compass-like signs, but none is wired into the stable nadsoliton axis landscape as a strict signed tilt with response/replay semantics.\n"
        "- Verdict: no current signed-source bridge to P2324 mu; no G1/G3 update, QW-2191 discharge, selector closure, or ToE closure.\n",
        encoding="utf-8",
    )
    print(json.dumps({"wrote": str(OUT.relative_to(REPO)), "theorem_fingerprint_sha256": theorem_fingerprint}, indent=2))


if __name__ == "__main__":
    main()
