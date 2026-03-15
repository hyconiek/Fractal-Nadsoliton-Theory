#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = GENERATED / "p444_current_strict_t168_vpsi_g4_g6_value_provider_repo_scan_probe.json"
OUT_SUMMARY = GENERATED / "p444_current_strict_t168_vpsi_g4_g6_value_provider_repo_scan_probe_summary.json"

# Avoid self-referential scan noise (always excluded, even for full scans).
ALWAYS_EXCLUDE_GLOBS = [
    "fundamental_action_reconstruction/P444_CURRENT_STRICT_T168_VPSI_G4_G6_VALUE_PROVIDER_REPO_SCAN_PROBE.md",
    "fundamental_action_reconstruction/p444_current_strict_t168_vpsi_g4_g6_value_provider_repo_scan_probe.py",
    "fundamental_action_reconstruction/generated/p444_current_strict_t168_vpsi_g4_g6_value_provider_repo_scan_probe.json",
    "fundamental_action_reconstruction/generated/p444_current_strict_t168_vpsi_g4_g6_value_provider_repo_scan_probe_summary.json",
]

# Default exclusions: huge vendor/cache dirs (repo-state hygiene scan, not dependency scan).
DEFAULT_EXCLUDE_GLOBS = [
    "**/.venv/**",
    "**/__pycache__/**",
    "**/.pytest_cache/**",
    "**/.mypy_cache/**",
    "**/.ruff_cache/**",
    "external_confirmatory_v2/**",
    "raw_strain_unfiltered/**",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "P444 repo-state scan: search the whole repo for any exported numeric value-provider instantiating "
            "T168 inputs vpsi[0..11], g4[0..11], g6[0..11] (length-12 arrays, vpsi nonzero)."
        )
    )
    p.add_argument(
        "--full-scan",
        action="store_true",
        help=(
            "Include vendor/cache directories (e.g. .venv, __pycache__). Still excludes P444 self artifacts to avoid scan noise."
        ),
    )
    p.add_argument(
        "--max-candidates",
        type=int,
        default=50,
        help="Max candidate JSON files to parse (after rg prefilter).",
    )
    return p.parse_args()


def run_rg_files(
    pattern: str,
    *,
    exclude_globs: list[str],
    multiline: bool = False,
    file_globs: list[str] | None = None,
) -> list[str]:
    globs = file_globs or ["*.json", "*.py", "*.md", "*.ipynb"]
    cmd = [
        "rg",
        "--hidden",
        *(["--multiline"] if multiline else []),
        *[arg for g in globs for arg in ("--glob", g)],
        *[arg for g in exclude_globs for arg in ("--glob", f"!{g}")],
        "-l",
        pattern,
        str(REPO),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    out = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    return out


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def load_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def classify_candidate(*, path: Path, obj: dict[str, Any]) -> dict[str, Any]:
    stage = obj.get("stage")
    status = obj.get("status")
    classification = obj.get("classification")

    s_stage = str(stage) if stage is not None else ""
    s_status = str(status) if status is not None else ""
    s_class = str(classification) if classification is not None else ""

    # Only treat as strict-derived if the artifact explicitly declares it.
    explicit_strict_derived = any(
        token in (s_status.upper(), s_class.lower())
        for token in ("STRICT_DERIVED", "strict_derived")
    )

    if "STRICT_EXTENSION_ONLY" in s_status.upper() or s_stage.upper().startswith("AX"):
        kind = "strict_extension_only_or_axiom_premise"
    elif "TOY" in s_status.upper() or obj.get("toy_instantiation") is True:
        kind = "toy_or_probe_instantiation"
    else:
        kind = "unclassified_numeric_instantiation"

    # Provide minimal context but avoid promotion.
    return {
        "path": str(path.relative_to(REPO)),
        "stage": stage,
        "status": status,
        "classification": classification,
        "kind": kind,
        "explicit_strict_derived_marker": bool(explicit_strict_derived),
    }


def extract_len12_numeric_list(obj: dict[str, Any], key: str) -> list[float] | None:
    arr = obj.get(key)
    if not (isinstance(arr, list) and len(arr) == 12):
        return None
    out: list[float] = []
    for x in arr:
        if not is_number(x):
            return None
        out.append(float(x))
    return out


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    args = parse_args()
    exclude_globs_used = list(ALWAYS_EXCLUDE_GLOBS) + ([] if args.full_scan else list(DEFAULT_EXCLUDE_GLOBS))

    number = r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"

    # Prefilter: find any JSON files with numeric array-valued vpsi/g4/g6 keys.
    files_vpsi = set(
        run_rg_files(rf"\"vpsi\"\s*:\s*\[\s*{number}", exclude_globs=exclude_globs_used, multiline=True, file_globs=["*.json"])
    )
    files_g4 = set(
        run_rg_files(rf"\"g4\"\s*:\s*\[\s*{number}", exclude_globs=exclude_globs_used, multiline=True, file_globs=["*.json"])
    )
    files_g6 = set(
        run_rg_files(rf"\"g6\"\s*:\s*\[\s*{number}", exclude_globs=exclude_globs_used, multiline=True, file_globs=["*.json"])
    )

    # Per-site key style (alternative instantiation forms).
    files_vpsi_i = set(
        run_rg_files(rf"\"vpsi\d+\"\s*:\s*{number}", exclude_globs=exclude_globs_used, file_globs=["*.json"])
    )
    files_g4_psi_i = set(
        run_rg_files(rf"\"g4_psi\d+\"\s*:\s*{number}", exclude_globs=exclude_globs_used, file_globs=["*.json"])
    )
    files_g6_psi_i = set(
        run_rg_files(rf"\"g6_psi\d+\"\s*:\s*{number}", exclude_globs=exclude_globs_used, file_globs=["*.json"])
    )

    candidate_json_files = sorted(files_vpsi & files_g4 & files_g6)
    alt_per_site_candidate_files = sorted(files_vpsi_i & files_g4_psi_i & files_g6_psi_i)

    parsed_candidates: list[dict[str, Any]] = []
    parse_errors: list[str] = []
    validated_t168_consumable: list[dict[str, Any]] = []

    for file_str in candidate_json_files[: args.max_candidates]:
        path = Path(file_str)
        obj = load_json(path)
        if obj is None or not isinstance(obj, dict):
            parse_errors.append(str(path.relative_to(REPO)))
            continue
        meta = classify_candidate(path=path, obj=obj)
        vpsi = extract_len12_numeric_list(obj, "vpsi")
        g4 = extract_len12_numeric_list(obj, "g4")
        g6 = extract_len12_numeric_list(obj, "g6")
        if vpsi is not None and g4 is not None and g6 is not None and all(x != 0.0 for x in vpsi):
            validated_t168_consumable.append(meta)
        parsed_candidates.append(meta)

    # Report alt per-site candidates (do not attempt to lift to arrays here; T168 expects arrays).
    alt_candidates: list[dict[str, Any]] = []
    for file_str in alt_per_site_candidate_files[: args.max_candidates]:
        path = Path(file_str)
        obj = load_json(path)
        if obj is None or not isinstance(obj, dict):
            continue
        alt_candidates.append(classify_candidate(path=path, obj=obj))

    any_strict_derived = any(c.get("explicit_strict_derived_marker") for c in validated_t168_consumable)
    decision_ready = bool(any_strict_derived)

    artifact = {
        "stage": "P444",
        "goal": "repo_scan_for_any_exported_numeric_value_provider_instantiating_T168_inputs_vpsi_g4_g6",
        "scan_scope": {
            "repo_root": str(REPO),
            "full_scan": bool(args.full_scan),
            "exclude_globs_used": exclude_globs_used,
        },
        "search_summary": {
            "json_files_with_numeric_array_vpsi": len(files_vpsi),
            "json_files_with_numeric_array_g4": len(files_g4),
            "json_files_with_numeric_array_g6": len(files_g6),
            "candidate_json_files_with_all_three_arrays": len(candidate_json_files),
            "candidate_json_files_with_all_three_arrays_sample": candidate_json_files[:10],
            "alt_per_site_candidate_json_files": len(alt_per_site_candidate_files),
            "alt_per_site_candidate_json_files_sample": alt_per_site_candidate_files[:10],
        },
        "candidates": {
            "parsed_candidate_json_files": parsed_candidates,
            "parse_errors": parse_errors,
            "validated_t168_consumable_candidates": validated_t168_consumable,
            "alt_per_site_style_candidates": alt_candidates,
        },
        "result": {
            "any_numeric_t168_consumable_candidate_found": bool(validated_t168_consumable),
            "any_explicit_strict_derived_t168_provider_found": bool(any_strict_derived),
            "decision_ready_from_repo_values": bool(decision_ready),
            "recommended_next_strict_target": "T168" if not decision_ready else None,
            "recommendation_reason": (
                "No explicitly strict-derived vpsi/g4/g6 value-provider object was found in a T168-consumable form; "
                "toy/probe/extension instantiations must not be promoted."
                if not decision_ready
                else "A strict-derived T168 provider appears to be present; manual provenance review still required."
            ),
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P444",
        "status": "PASS_SCAN_READY_NO_STRICT_DERIVED_T168_PROVIDER_FOUND" if not decision_ready else "PASS_FOUND_STRICT_DERIVED_T168_PROVIDER_CANDIDATE_REVIEW_REQUIRED",
        "decision_ready_from_repo_values": bool(decision_ready),
        "validated_t168_consumable_candidates_count": len(validated_t168_consumable),
        "any_explicit_strict_derived_t168_provider_found": bool(any_strict_derived),
        "recommended_next_strict_target": "T168" if not decision_ready else None,
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

