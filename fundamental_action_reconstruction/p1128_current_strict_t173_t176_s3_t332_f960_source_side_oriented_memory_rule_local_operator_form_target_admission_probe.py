#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-04-25"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_S3 = ROOT / "S3_CURRENT_INFORMATION_COMPRESSION_DECOMPRESSION_ROUTE_CLASSIFICATION_NOTE.md"
IN_T332 = ROOT / "T332_CURRENT_STRICT_T173_T176_EXISTING_S3_F966_F972_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_TARGET_SPEC.md"
IN_N966 = ROOT / "N966_CURRENT_STRICT_T173_T176_EXISTING_S3_F966_F972_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_TARGET_ADMISSION_THEOREM.md"
IN_F969 = ROOT / "F969_CURRENT_STRICT_T173_T176_MINIMAL_ORIENTED_NONRECIPROCAL_DEPHASING_NEW_IMPORT_BOUNDARY_PACKET.md"
IN_N717 = ROOT / "N717_CURRENT_STRICT_T176_SOURCE_TOPOLOGY_BASIS_FREE_QW2191_SAFE_PROVIDER_NONUPGRADE_BOUNDARY_THEOREM.md"
IN_N718 = ROOT / "N718_CURRENT_STRICT_T177_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_NONEXPORT_BOUNDARY_THEOREM.md"
IN_N719 = ROOT / "N719_CURRENT_STRICT_T178_SOURCE_TOPOLOGY_TO_ATLAS_CHART_SEED_SELECTION_BRIDGE_NONEXPORT_BOUNDARY_THEOREM.md"
IN_C6 = ROOT / "C6_PROJECTED_SECOND_VARIATION_SOURCE_AUDIT.md"
IN_F456 = ROOT / "F456_CURRENT_STRICT_SIGMA_INT_ORIENTATION_SLICE_TO_A1_PAIR1_OPERATOR_BRIDGE_PACKET.md"

OUT_JSON = GENERATED / "p1128_current_strict_t173_t176_s3_t332_f960_source_side_oriented_memory_rule_local_operator_form_target_admission_probe.json"
OUT_SUMMARY = GENERATED / "p1128_current_strict_t173_t176_s3_t332_f960_source_side_oriented_memory_rule_local_operator_form_target_admission_probe_summary.json"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
PARENT_TARGET = "SourceSideOrientedMemoryRuleTarget_against_ResidualDatumPair12OrbitDirectionSelectionBridge_v1"
TARGET_NAME = "LocalSourceSideOrientedMemoryRuleOperatorFormTarget_against_ResidualDatumPair12OrbitDirectionSelectionBridge_v1"


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_existing_target_hits() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded = {
        Path(__file__).name,
        "P1128_CURRENT_STRICT_T173_T176_EXISTING_S3_T332_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ADMISSION_PROBE.md",
        "T333_CURRENT_STRICT_T173_T176_EXISTING_S3_T332_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_SPEC.md",
        "N967_CURRENT_STRICT_T173_T176_EXISTING_S3_T332_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ADMISSION_THEOREM.md",
        OUT_JSON.name,
        OUT_SUMMARY.name,
    }
    hits: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if TARGET_NAME in text:
                hits.append(rel(path))
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_S3, IN_T332, IN_N966, IN_F969, IN_N717, IN_N718, IN_N719, IN_C6, IN_F456]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1128",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    s3_text = load_text(IN_S3)
    t332_text = load_text(IN_T332)
    n966_text = load_text(IN_N966)
    f969_text = load_text(IN_F969)
    n717_text = load_text(IN_N717)
    n718_text = load_text(IN_N718)
    n719_text = load_text(IN_N719)
    c6_text = load_text(IN_C6)
    f456_text = load_text(IN_F456)
    existing_hits = scan_existing_target_hits()

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    s3_requires_local_operator_form_and_hooks = all(
        needle in s3_text
        for needle in [
            "explicit **local** operator form",
            "admissible well-posedness control on the RG side",
            "admissible positivity control on the QFT side",
            "Without those locality and positivity hooks, the idea remains heuristic only.",
        ]
    )

    t332_n966_parent_target_frozen = all(
        needle in t332_text + n966_text
        for needle in [
            PARENT_TARGET,
            ACTIVE_BRIDGE,
            "future_rg_locality_audit_hook_required := yes",
            "future_qft_positivity_audit_hook_required := yes",
            "candidate_provider_class_seed_only",
        ]
    )

    f969_oriented_memory_rule_component_explicit = all(
        needle in f969_text
        for needle in [
            "oriented_memory_rule_component := oriented_memory_rule",
            ACTIVE_BRIDGE,
            "admissible_as_candidate_provider_class_seed_only := yes",
        ]
    )

    chart_sensitive_geometry_still_required = all(
        needle in n717_text + n718_text + n719_text
        for needle in [
            "chart-sensitive transported signed flux/current-like section",
            "chart-sensitive transported flux/current-like section",
            "physics-facing but still chart-blind",
            "source-topology to atlas chart-seed selection bridge",
        ]
    )

    repo_already_admits_hook_only_narrowing = all(
        needle in c6_text
        for needle in [
            "positivity only in declared strict scope.",
            "packet-ready source tuple",
            "no_explicit_positivity_certified_second_variation_block_on_that_exported_subspace",
        ]
    )

    repo_already_admits_operator_bridge_object = all(
        needle in f456_text
        for needle in [
            "operator-level bridge object",
            "A_1(pair1)",
            "strict-core, slot-free operator object",
        ]
    )

    no_existing_local_operator_form_target_export = len(existing_hits) == 0

    strongest_honest_next_move_is_freeze_local_operator_form_target = (
        s3_requires_local_operator_form_and_hooks
        and t332_n966_parent_target_frozen
        and f969_oriented_memory_rule_component_explicit
        and chart_sensitive_geometry_still_required
        and repo_already_admits_hook_only_narrowing
        and repo_already_admits_operator_bridge_object
        and no_existing_local_operator_form_target_export
    )

    add_check(
        "s3_requires_local_operator_form_and_hooks",
        s3_requires_local_operator_form_and_hooks,
        True,
        "S3 already requires explicit local operator form plus locality/positivity hooks to avoid heuristic floating.",
    )
    add_check(
        "t332_n966_parent_target_frozen",
        t332_n966_parent_target_frozen,
        True,
        "T332/N966 already freeze the parent oriented-memory-rule target below bridge realization.",
    )
    add_check(
        "f969_oriented_memory_rule_component_explicit",
        f969_oriented_memory_rule_component_explicit,
        True,
        "F969 already keeps oriented_memory_rule explicit as one irreducible novelty component.",
    )
    add_check(
        "chart_sensitive_geometry_still_required",
        chart_sensitive_geometry_still_required,
        True,
        "N717/N718/N719 keep the required chart-sensitive transported-section geometry explicit.",
    )
    add_check(
        "repo_already_admits_hook_only_narrowing",
        repo_already_admits_hook_only_narrowing,
        True,
        "C6 already shows the repo admits one narrower move that freezes source tuples and missing positivity hooks without faking the certificate.",
    )
    add_check(
        "repo_already_admits_operator_bridge_object",
        repo_already_admits_operator_bridge_object,
        True,
        "F456 already shows the repo admits one operator-level bridge object below closure.",
    )
    add_check(
        "no_existing_local_operator_form_target_export",
        no_existing_local_operator_form_target_export,
        True,
        "No current exact target already isolates the local operator-form target beneath T332.",
    )
    add_check(
        "strongest_honest_next_move_is_freeze_local_operator_form_target",
        strongest_honest_next_move_is_freeze_local_operator_form_target,
        True,
        "Therefore the strongest honest next export is one exact local operator-form target beneath T332 with locality/positivity hooks and chart-sensitive compatibility.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_S3_T332_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ADMISSION_AUDITED"
        if not blocking and strongest_honest_next_move_is_freeze_local_operator_form_target
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_S3_T332_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ADMISSION_AUDIT"
    )

    artifact = {
        "stage": "P1128",
        "status": status,
        "as_of": AS_OF,
        "parent_target_name": PARENT_TARGET,
        "target_name": TARGET_NAME,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "s3_requires_local_operator_form_and_hooks": s3_requires_local_operator_form_and_hooks,
        "t332_n966_parent_target_frozen": t332_n966_parent_target_frozen,
        "f969_oriented_memory_rule_component_explicit": f969_oriented_memory_rule_component_explicit,
        "chart_sensitive_geometry_still_required": chart_sensitive_geometry_still_required,
        "repo_already_admits_hook_only_narrowing": repo_already_admits_hook_only_narrowing,
        "repo_already_admits_operator_bridge_object": repo_already_admits_operator_bridge_object,
        "no_existing_local_operator_form_target_export": no_existing_local_operator_form_target_export,
        "current_repo_already_exports_same_target_hits": existing_hits,
        "strongest_honest_next_move_is_freeze_local_operator_form_target": strongest_honest_next_move_is_freeze_local_operator_form_target,
        "counts_as_exact_bridge_reduction": False,
        "counts_as_lawful_supplier": False,
        "counts_as_strict_physical_orientation_datum": False,
        "counts_as_t183_discharge": False,
        "counts_as_t176_discharge": False,
        "counts_as_qw2191_discharge": False,
        "counts_as_qw2580_resolution": False,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "parent_target_name": artifact["parent_target_name"],
        "target_name": artifact["target_name"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "strongest_honest_next_move_is_freeze_local_operator_form_target": artifact[
            "strongest_honest_next_move_is_freeze_local_operator_form_target"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
