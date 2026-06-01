#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.json"
MD = GEN / "p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.md"

SOURCE_FILES = {
    "P1563_KERNEL_TO_EOM_EXPORT": GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json",
    "P1622_FULL_STRICT_LAGRANGIAN_EOM_EXPORT": GEN / "p1622_s572_full_strict_lagrangian_density_and_eom_summary.json",
    "P1653_NONSKELETON_FULL_LAGRANGIAN": GEN / "p1653_s603_strict_full_lagrangian_nonskeleton_and_bidirectional_obligation_summary.json",
    "P1688_SECTOR_EOM_EXPORT": GEN / "p1688_s638_strict_full_lagrangian_to_sector_eom_export_summary.json",
    "P1693_MULTISECTOR_SYMPY_BRIDGE": GEN / "p1693_s643_strict_full_lagrangian_multisector_sympy_eom_bridge.json",
    "P1866_SYMBOLIC_FULL_LAGRANGIAN_EXPORT": GEN / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json",
    "P2030_TENSOR_PROJECTION_AUDIT": GEN / "p2030_s980_strict_tensor_resolved_projection_source_audit.json",
    "P2031_B1_TENSOR_COMPONENT_SCAFFOLD": GEN / "p2031_s981_strict_b1_tensor_component_profile_table_scaffold.json",
    "P2032_METRIC_GAUGE_PROJECTION_RULE_AUDIT": GEN / "p2032_s982_strict_b1_metric_gauge_component_projection_rule_audit.json",
    "P2033_CURVED_B1_ANSATZ_NONAVAILABILITY": GEN / "p2033_s983_strict_curved_b1_metric_ansatz_nonavailability_theorem.json",
    "P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT": GEN / "p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.json",
    "P2314_FULL_EOM_LAGRANGIAN_INVENTORY": GEN / "p2314_s1264_strict_full_eom_lagrangian_symmetry_vortex_inventory_probe.json",
    "P2315_SCHEMATIC_LAGRANGIAN_SPECTRUM": GEN / "p2315_s1265_strict_schematic_lagrangian_eom_kernel_spectrum_probe.json",
    "RELEASE7_SCHEMATIC_CANDIDATE_LAYER": REPO / "RELEASE_7_STRICT_PROJECTIVE_OPERATIONAL_MODEL_BRIEF.md",
}

GREP_PATTERNS = (
    "Lagrang",
    "L_total",
    "Euler-Lagrange",
    "full Lagrangian",
    "strict full",
    "termwise",
    "sector EOM",
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


def first_present(data: dict[str, Any], keys: list[str], default: Any = None) -> Any:
    for key in keys:
        if key in data:
            return data[key]
    return default


def short_string(value: Any, limit: int = 420) -> str:
    if isinstance(value, (dict, list)):
        text = json.dumps(value, sort_keys=True, ensure_ascii=False)
    else:
        text = str(value)
    return text if len(text) <= limit else text[: limit - 3] + "..."


def collect_repo_search_hits() -> list[dict[str, Any]]:
    candidate_paths = [p for p in SOURCE_FILES.values() if p.exists()]
    extra_paths = sorted(ROOT.glob("p*_strict*lagrangian*.py"))[:20]
    hits: list[dict[str, Any]] = []
    for path in sorted(set(candidate_paths + extra_paths)):
        rel = path.relative_to(REPO).as_posix()
        text = read_text(path)
        lowered = text.lower()
        count = sum(lowered.count(pattern.lower()) for pattern in GREP_PATTERNS)
        if count == 0:
            continue
        first_line = 0
        first_excerpt = ""
        for index, line in enumerate(text.splitlines(), start=1):
            if any(pattern.lower() in line.lower() for pattern in GREP_PATTERNS):
                first_line = index
                first_excerpt = line.strip()[:220]
                break
        hits.append({
            "path": rel,
            "pattern_hit_count": count,
            "first_hit_line": first_line,
            "first_hit_excerpt": first_excerpt,
        })
    hits.sort(key=lambda row: (-int(row["pattern_hit_count"]), row["path"]))
    return hits


def candidate_rows(artifacts: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    ranking = {
        "P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT": 100,
        "P1693_MULTISECTOR_SYMPY_BRIDGE": 90,
        "P1688_SECTOR_EOM_EXPORT": 85,
        "P1653_NONSKELETON_FULL_LAGRANGIAN": 80,
        "P1622_FULL_STRICT_LAGRANGIAN_EOM_EXPORT": 70,
        "P1866_SYMBOLIC_FULL_LAGRANGIAN_EXPORT": 68,
        "P1563_KERNEL_TO_EOM_EXPORT": 55,
        "P2315_SCHEMATIC_LAGRANGIAN_SPECTRUM": 35,
        "RELEASE7_SCHEMATIC_CANDIDATE_LAYER": 30,
    }
    for source_id, path in SOURCE_FILES.items():
        data = artifacts.get(source_id, {})
        if not path.exists():
            continue
        if path.suffix == ".md":
            status = "TEXT_SOURCE_PRESENT"
        else:
            status = data.get("status", "TEXT_SOURCE_PRESENT")
        candidate_form = first_present(
            data,
            [
                "full_lagrangian_density_nonskeleton",
                "full_lagrangian_density_anchor",
                "full_lagrangian_density",
                "L_total_context",
                "lagrangian_density",
                "full_lagrangian_non_skeleton",
            ],
            "schematic/text source",
        )
        if source_id == "P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT":
            candidate_form = data.get("eom_execution_results", {}).get("full_ltotal", "missing full_ltotal")
        if source_id == "P2315_SCHEMATIC_LAGRANGIAN_SPECTRUM":
            candidate_form = data.get("strict_schematic_lagrangian_eom_kernel_spectrum_probe", {}).get(
                "schematic_eom_export", {}
            )
        if source_id == "RELEASE7_SCHEMATIC_CANDIDATE_LAYER":
            candidate_form = "schematic nadsoliton core Lagrangian representative in Release 7 text"
        rows.append({
            "source_id": source_id,
            "path": path.relative_to(REPO).as_posix(),
            "sha256": sha256_file(path),
            "status": status,
            "strength_rank": ranking.get(source_id, 10),
            "candidate_form_excerpt": short_string(candidate_form),
        })
    rows.sort(key=lambda row: (-int(row["strength_rank"]), row["source_id"]))
    return rows


def termwise_audit_summary(p2086: dict[str, Any]) -> dict[str, Any]:
    results = p2086.get("eom_execution_results", {})
    terms = results.get("lagrangian_terms", {})
    residuals = results.get("symbolic_recomposition_residual", {})
    numeric = results.get("numeric_probe_residual", {})
    fields = sorted(results.get("eom_full", {}).keys())
    return {
        "source_id": "P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT",
        "loaded": p2086.get("packet_id") == "P2086",
        "result_kind": p2086.get("result_kind"),
        "term_count": len(terms),
        "field_count": len(fields),
        "fields": fields,
        "symbolic_zero_fields": sorted(k for k, v in residuals.items() if v in ("Integer(0)", "0", 0)),
        "numeric_zero_fields": sorted(k for k, v in numeric.items() if str(v) == "0"),
        "all_symbolic_residual_zero": bool(residuals) and all(v in ("Integer(0)", "0", 0) for v in residuals.values()),
        "all_numeric_residual_zero": bool(numeric) and all(str(v) == "0" for v in numeric.values()),
        "termwise_execution_status": results.get("termwise_execution_status"),
        "theorem_grade_tensor_status": "OPEN_C3_STILL_OPEN__NON_THEOREM_GRADE",
    }


def gap_summary(artifacts: dict[str, dict[str, Any]]) -> dict[str, Any]:
    gap_ids = [
        "P2030_TENSOR_PROJECTION_AUDIT",
        "P2031_B1_TENSOR_COMPONENT_SCAFFOLD",
        "P2032_METRIC_GAUGE_PROJECTION_RULE_AUDIT",
        "P2033_CURVED_B1_ANSATZ_NONAVAILABILITY",
    ]
    blockers = []
    for source_id in gap_ids:
        data = artifacts.get(source_id, {})
        blockers.append({
            "source_id": source_id,
            "status": data.get("status"),
            "result_kind": data.get("result_kind"),
            "gatekeeper_false_or_missing_ready_flags": sorted(
                key for key, value in data.get("gatekeeper_checks", {}).items()
                if (key.endswith("ready") or key.endswith("exported") or key.endswith("claimed")) and value is False
            ),
        })
    return {
        "tensor_gap_packet_count": len(blockers),
        "all_gap_packets_loaded": all(bool(artifacts.get(source_id, {}).get("status")) for source_id in gap_ids),
        "blockers": blockers,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {source_id: load_json(path) for source_id, path in SOURCE_FILES.items() if path.suffix == ".json"}
    artifacts["RELEASE7_SCHEMATIC_CANDIDATE_LAYER"] = {"status": "TEXT_SOURCE_PRESENT"}

    rows = candidate_rows(artifacts)
    search_hits = collect_repo_search_hits()
    p2086_summary = termwise_audit_summary(artifacts["P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT"])
    gaps = gap_summary(artifacts)

    p1653 = artifacts["P1653_NONSKELETON_FULL_LAGRANGIAN"]
    p1693 = artifacts["P1693_MULTISECTOR_SYMPY_BRIDGE"]
    p2314 = artifacts["P2314_FULL_EOM_LAGRANGIAN_INVENTORY"]
    p2315 = artifacts["P2315_SCHEMATIC_LAGRANGIAN_SPECTRUM"]

    strongest_form = {
        "selected_working_form_id": "P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT",
        "selected_structural_anchor_id": "P1653_NONSKELETON_FULL_LAGRANGIAN",
        "selected_covariant_sector_anchor_id": "P1693_MULTISECTOR_SYMPY_BRIDGE",
        "reason": "P2086 is the strongest current computational artifact because it executes a termwise Euler-Lagrange recomposition check for the reduced L_total; P1653/P1693 remain the stronger structural/covariant-sector anchors for the named full L_total decomposition.",
        "canonical_nonskeleton_decomposition": p1653.get("full_lagrangian_density_nonskeleton", {}).get(
            "L_total", "L_strict_scalar + L_SM_gauge + L_SM_fermions + L_SM_higgs + L_GR + L_mix"
        ),
        "covariant_sector_terms": p1693.get("full_lagrangian_density_anchor", {}),
        "reduced_computational_ltotal_srepr": artifacts["P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT"].get(
            "eom_execution_results", {}
        ).get("full_ltotal", ""),
        "not_selected_as_best_full_lagrangian": [
            {
                "source_id": "RELEASE7_SCHEMATIC_CANDIDATE_LAYER",
                "reason": "schematic candidate core layer only; P2315 derives only schematic Z12 EOM and spectrum",
            },
            {
                "source_id": "P1622_FULL_STRICT_LAGRANGIAN_EOM_EXPORT",
                "reason": "older compact full-density export superseded for coverage auditing by nonskeleton/sector/termwise artifacts",
            },
        ],
    }

    current_limitations = {
        "best_current_status": "BEST_WORKING_LTOTAL_IDENTIFIED__FULL_TASK3_THEOREM_STILL_OPEN",
        "full_eom_answer": "NO_FULL_TENSOR_RESOLVED_EOM_FOR_TASK3_YET",
        "full_lagrangian_answer": "YES_SECTOR_DECOMPOSED_WORKING_LTOTAL_EXISTS__NO_THEOREM_GRADE_GLOBAL_FULL_LAGRANGIAN_FOR_TASK3",
        "computational_coverage": p2086_summary,
        "tensor_and_metric_gaps": gaps,
        "p2314_status": p2314.get("status"),
        "p2315_status": p2315.get("status"),
        "g1_g3_update": "NONE__P2282_G1_G3_REMAIN_OPEN_UNCHANGED",
    }

    theorem_export = {
        "theorem_name": "P2316 current best strict Lagrangian/EOM coverage audit",
        "claim": "The current repository exports a best working strict L_total as a sector-decomposed/nonskeleton plus reduced termwise-computational object, but it does not export the full tensor-resolved theorem-grade EOM/Lagrangian needed to update Task-3 G1/G3.",
        "proof_witnesses": [
            "P1653 exports the nonskeleton L_total sector decomposition.",
            "P1693 exports covariant sector anchors and a reduced multisector Sympy EOM bridge.",
            "P2086 verifies termwise Euler-Lagrange recomposition residuals for the reduced L_total.",
            "P2030-P2033 keep tensor projection, component table, projection rule, and curved B1 ansatz unavailable.",
            "P2314/P2315 answer the full-EOM/Lagrangian and schematic-spectrum questions without claiming selector or G1/G3 closure.",
        ],
        "blocker": "C3/tensor-resolved metric-background theorem remains open on the EOM/Lagrangian track; selector and policy-margin sources are separate parallel blockers and are not prerequisites for continuing EOM export.",
        "strict_guardrails": {
            "no_legacy_kernel_role_transfer": True,
            "no_selector_premise_added": True,
            "no_qw2191_discharge_claimed": True,
            "no_toe_closure_claimed": True,
        },
    }

    probe = {
        "probe_id": "P2316_S1266_STRICT_CURRENT_BEST_LAGRANGIAN_EOM_COVERAGE_AUDIT",
        "repo_grep_audit": {
            "search_patterns": list(GREP_PATTERNS),
            "hit_count": len(search_hits),
            "top_hits": search_hits[:20],
        },
        "candidate_lagrangian_sources_ranked": rows,
        "strongest_current_lagrangian_form": strongest_form,
        "current_limitations_and_eom_coverage": current_limitations,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "repo_lagrangian_grep_hits_found": len(search_hits) >= 5,
        "candidate_sources_loaded": len(rows) >= 10,
        "p1653_nonskeleton_anchor_loaded": p1653.get("checkpoint") == "P1653_S603",
        "p1693_multisector_bridge_loaded": p1693.get("checkpoint") == "P1693_S643",
        "p2086_termwise_audit_loaded": p2086_summary["loaded"],
        "p2086_symbolic_recomposition_zero": p2086_summary["all_symbolic_residual_zero"],
        "p2086_numeric_probe_zero": p2086_summary["all_numeric_residual_zero"],
        "p2086_non_theorem_grade_preserved": p2086_summary["termwise_execution_status"] == "TERMWISE_EXECUTED_NON_THEOREM_GRADE",
        "tensor_gap_packets_loaded": gaps["all_gap_packets_loaded"],
        "p2314_inventory_loaded": p2314.get("packet_id") == "P2314",
        "p2315_schematic_spectrum_loaded": p2315.get("packet_id") == "P2315",
        "best_working_ltotal_identified": strongest_form["selected_working_form_id"] == "P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT",
        "full_task3_theorem_not_claimed": True,
        "strict_g1_g3_not_updated": True,
        "no_legacy_kernel_role_transfer": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2316_s1266_v1",
        "packet_id": "P2316",
        "stage_id": "S1266",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_BEST_WORKING_LTOTAL_IDENTIFIED_FULL_TASK3_EOM_THEOREM_STILL_MISSING",
        "result_kind": "STRICT_REPO_GREP_AND_COMPUTATIONAL_LAGRANGIAN_EOM_COVERAGE_AUDIT_NO_G1_G3_UPDATE",
        "strict_current_best_lagrangian_eom_coverage_audit_probe": probe,
        "recommended_next_honest_step": "Continue the EOM/Lagrangian track independently: export tensor-resolved metric/background equations for the selected L_total, or prove a nonexistence theorem for that tensor/background lift. Keep selector/QW-2191 closure as a separate parallel track.",
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")

    md_lines = [
        "# P2316 S1266: current best strict Lagrangian/EOM coverage audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Repository grep result",
        f"- Lagrangian/EOM search hits recorded: `{len(search_hits)}`.",
        "- Best current computational working form: `P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT`.",
        "- Best current structural sector anchor: `P1653_NONSKELETON_FULL_LAGRANGIAN` plus `P1693_MULTISECTOR_SYMPY_BRIDGE`.",
        "",
        "## Selected current L_total",
        f"`{strongest_form['canonical_nonskeleton_decomposition']}`",
        "",
        "## Computational check",
        f"- P2086 term count: `{p2086_summary['term_count']}`.",
        f"- P2086 fields varied: `{', '.join(p2086_summary['fields'])}`.",
        f"- Symbolic recomposition residual zero: `{p2086_summary['all_symbolic_residual_zero']}`.",
        f"- Numeric probe residual zero: `{p2086_summary['all_numeric_residual_zero']}`.",
        f"- Termwise status: `{p2086_summary['termwise_execution_status']}`.",
        "",
        "## Honest limitation",
        f"- Full tensor-resolved EOM answer: `{current_limitations['full_eom_answer']}`.",
        f"- Full Lagrangian answer: `{current_limitations['full_lagrangian_answer']}`.",
        f"- Tensor gap packet count: `{gaps['tensor_gap_packet_count']}`.",
        "- Task-3 G1/G3 update: `NONE`; the closure rows remain open.",
        "",
        "## Theorem fingerprint",
        f"`{probe['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
