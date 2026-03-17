#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = GENERATED / "p672_current_strict_t172_readiness_probe.json"
OUT_SUMMARY = GENERATED / "p672_current_strict_t172_readiness_probe_summary.json"


def read_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads(read_text(repo_relative_path))


def exists(repo_relative_path: str) -> bool:
    return (REPO / repo_relative_path).exists()


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq_summaries = {
        "F469": "fundamental_action_reconstruction/generated/f469_current_strict_global_selector_atlas_and_transition_object_export_packet_summary.json",
        "F470": "fundamental_action_reconstruction/generated/f470_current_strict_global_projective_selector_state_object_export_packet_summary.json",
        "F473": "fundamental_action_reconstruction/generated/f473_current_strict_t164_z12_generator_orientation_canonical_fixing_datum_export_packet_summary.json",
        "F474": "fundamental_action_reconstruction/generated/f474_current_strict_t171_global_directed_selector_state_datum_export_packet_summary.json",
        "N550": "fundamental_action_reconstruction/generated/n550_current_strict_global_selector_bridge_operator_promotion_from_seed_v1_chain_on_c_v1_discharge_theorem_summary.json",
        "N551": "fundamental_action_reconstruction/generated/n551_current_strict_global_selector_reduction_operator_promotion_from_seed_v1_chain_on_c_v1_discharge_theorem_summary.json",
        "N552": "fundamental_action_reconstruction/generated/n552_current_strict_global_selector_output_operator_promotion_from_seed_v1_chain_on_c_v1_discharge_theorem_summary.json",
        "N553": "fundamental_action_reconstruction/generated/n553_current_strict_global_downstream_completion_branch_discharge_for_promoted_seed_v1_chain_on_c_v1_discharge_theorem_summary.json",
    }

    loaded: dict[str, dict[str, Any] | None] = {}
    load_errors: dict[str, str] = {}
    for key, path in prereq_summaries.items():
        if not exists(path):
            loaded[key] = None
            load_errors[key] = f"missing: {path}"
            continue
        try:
            loaded[key] = load_json(path)
        except Exception as exc:  # pragma: no cover - probe robustness
            loaded[key] = None
            load_errors[key] = f"failed_to_parse_json: {path}: {exc}"

    t172_spec_path = "fundamental_action_reconstruction/T172_CURRENT_STRICT_GLOBAL_QW2191_DISCHARGE_AND_SELECTOR_CLOSURE_TARGET_SPEC.md"
    t172_spec_text = read_text(t172_spec_path) if exists(t172_spec_path) else ""

    def get_bool(obj: dict[str, Any] | None, dotted_key: str) -> bool | None:
        if obj is None:
            return None
        cur: Any = obj
        for part in dotted_key.split("."):
            if not isinstance(cur, dict) or part not in cur:
                return None
            cur = cur[part]
        if isinstance(cur, bool):
            return cur
        if isinstance(cur, (int, float)):
            return bool(cur)
        return None

    f469 = loaded["F469"]
    f470 = loaded["F470"]
    f473 = loaded["F473"]
    f474 = loaded["F474"]
    n550 = loaded["N550"]
    n551 = loaded["N551"]
    n552 = loaded["N552"]
    n553 = loaded["N553"]

    checks_spec = [
        {
            "id": "t172_target_spec_present",
            "actual": exists(t172_spec_path),
            "expected": True,
            "meaning": "T172 target spec exists (this probe audits readiness; it does not discharge T172).",
        },
        {
            "id": "t172_target_spec_no_false_pass_marked",
            "actual": "NO_FALSE_PASS" in t172_spec_text,
            "expected": True,
            "meaning": "T172 spec explicitly carries the no-false-pass discipline tag.",
        },
        {
            "id": "t170_discharged",
            "actual": get_bool(f469, "t170_discharged"),
            "expected": True,
            "meaning": "Global atlas + transition objects on C_v1 are exported (T170 discharged via F469/N515).",
        },
        {
            "id": "global_projective_state_exported",
            "actual": bool(get_bool(f470, "no_false_pass")) if f470 is not None else None,
            "expected": True,
            "meaning": "Global projective selector state object on C_v1 is exported (F470/N516).",
        },
        {
            "id": "t164_fixing_datum_exported_premise_based",
            "actual": (get_bool(f473, "t164_discharged") is True) and (get_bool(f473, "premise_based") is True),
            "expected": True,
            "meaning": "T164 fixing datum is exported and premise-based (explicit sign/orientation convention; no Aut(Z_12)-invariant overclaim).",
        },
        {
            "id": "t171_discharged_directed_state_exported",
            "actual": get_bool(f474, "t171_discharged"),
            "expected": True,
            "meaning": "Global directed/sign-sensitive selector state object is exported in the declared premise-tracked scope (T171 via F474/N524).",
        },
        {
            "id": "seed_v1_global_promotions_present",
            "actual": all(
                get_bool(obj, "theorem_result.discharged") is True
                for obj in (n550, n551, n552, n553)
                if obj is not None
            )
            and all(obj is not None for obj in (n550, n551, n552, n553)),
            "expected": True,
            "meaning": "Seed-v1 strict-core internal selector-source lane is promoted to global objects (N550–N553).",
        },
        {
            "id": "seed_v1_global_promotions_projector_section_level_only",
            "actual": all(
                get_bool(obj, "theorem_result.projector_section_level_only") is True
                for obj in (n550, n551, n552, n553)
                if obj is not None
            )
            and all(obj is not None for obj in (n550, n551, n552, n553)),
            "expected": True,
            "meaning": "Global promotions are explicitly projector/section-level only (respect N512 boundary; no operator-level groupoid promotion).",
        },
        {
            "id": "n553_residual_sign_gauge_explicit",
            "actual": get_bool(n553, "theorem_result.residual_sign_gauge_explicit"),
            "expected": True,
            "meaning": "Residual Z2 sign gauge is explicitly tracked at global level for promoted seed-v1 chain (N553).",
        },
    ]

    checks = []
    blocking_mismatches: list[str] = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )
        if not ok:
            blocking_mismatches.append(item["id"])

    # T172 missing target objects (expected on current repo state).
    closure_candidate_paths = [
        "fundamental_action_reconstruction/generated/selector_closure_global_c_v1_projective_strict_v1.json",
        "fundamental_action_reconstruction/generated/selector_closure_global_c_v1_directed_strict_v1.json",
    ]
    closure_objects_found = [p for p in closure_candidate_paths if exists(p)]

    # Current exported theorem summaries explicitly keep closure/QW2191 open (info-only).
    closure_claims_info = {
        "strict_core_selector_closure_claimed": bool(get_bool(n553, "theorem_result.strict_core_selector_closure")),
        "QW2191_discharge_claimed": bool(get_bool(n553, "theorem_result.QW2191_discharge")),
    }

    if blocking_mismatches:
        status = "P672_NOT_READY_MISSING_OR_FAILED_PREREQS_FOR_T172"
        artifact: dict[str, Any] = {
            "stage": "P672",
            "lane": "t172_readiness_probe_only",
            "status": status,
            "checks": checks,
            "blocking_mismatches": blocking_mismatches,
            "load_errors": load_errors,
            "closure_objects_found": closure_objects_found,
            "closure_claims_info": closure_claims_info,
            "no_false_pass": True,
        }
        summary: dict[str, Any] = {
            "stage": "P672",
            "lane": artifact["lane"],
            "status": status,
            "blocking_mismatches": blocking_mismatches,
            "no_false_pass": True,
        }
    else:
        if closure_objects_found:
            status = "P672_REVIEW_REQUIRED_CLOSURE_OBJECT_PRESENT_T172_PARTIALLY_STARTED"
        else:
            status = "P672_PREREQS_PRESENT_T172_OPEN_NO_GLOBAL_SELECTOR_CLOSURE_OBJECT_EXPORTED_YET"

        missing_target_objects = []
        if not closure_objects_found:
            missing_target_objects.append("SelectorClosure_global_C_v1_(projective|directed)_strict_v1 (choose one declared scope)")
        missing_target_objects.append("Theorem-level uniqueness/boundary statement resolving the strict meaning of global QW-2191 discharge in the chosen closure scope")

        artifact = {
            "stage": "P672",
            "lane": "t172_readiness_probe_only",
            "status": status,
            "checks": checks,
            "blocking_mismatches": [],
            "closure_objects_found": closure_objects_found,
            "closure_claims_info": closure_claims_info,
            "t172_missing_target_objects": missing_target_objects,
            "recommended_next_actions": [
                "Export one explicit global closure object on C_v1 (projective OR directed; keep sign discipline explicit).",
                "Export one theorem-level uniqueness/boundary statement for global QW-2191 in the chosen closure scope.",
            ],
            "hard_limits": [
                "this_is_a_readiness_probe_only",
                "no_strict_core_selector_closure_claim",
                "no_global_QW2191_discharge_claim",
                "no_ToE_closure_claim",
            ],
            "no_false_pass": True,
        }
        summary = {
            "stage": "P672",
            "lane": artifact["lane"],
            "status": status,
            "t172_missing_target_objects": missing_target_objects,
            "closure_objects_found": closure_objects_found,
            "closure_claims_info": closure_claims_info,
            "no_false_pass": True,
        }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

