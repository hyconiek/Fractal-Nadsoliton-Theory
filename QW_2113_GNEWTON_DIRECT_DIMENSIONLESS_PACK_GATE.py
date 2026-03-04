#!/usr/bin/env python3
"""
QW-2113: strict direct-dimensionless G bridge candidate-pack gate.

Purpose:
- validate candidate direct external dimensionless bridge input,
- enforce QW-2108 strict contract (no backsolve from G_SI),
- emit strict-ready artifacts for downstream QW-2101/QW-2103/QW-2092.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent

QW2108_JSON = ROOT / "report_qw2108_gnewton_dimensionless_acquisition_spec.json"

CAND_JSON_DEFAULT = ROOT / "external_gnewton_bridge_qw2113_direct_dimensionless_candidate.json"
CAND_META_DEFAULT = ROOT / "external_gnewton_bridge_qw2113_direct_dimensionless_candidate.metadata.json"

TEMPLATE_JSON = ROOT / "external_gnewton_bridge_qw2113_direct_dimensionless_candidate.template.json"
TEMPLATE_META = ROOT / "external_gnewton_bridge_qw2113_direct_dimensionless_candidate.metadata.template.json"

READY_JSON = ROOT / "external_gnewton_bridge_qw2101.direct_dimensionless_ready.json"
READY_META = ROOT / "external_gnewton_bridge_qw2101.direct_dimensionless_ready.metadata.json"

OUT_JSON = ROOT / "report_qw2113_gnewton_direct_dimensionless_pack_gate.json"
OUT_MD = ROOT / "RAPORT_QW2113_GNEWTON_DIRECT_DIMENSIONLESS_PACK_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _is_https_url(url: str) -> bool:
    return str(url).strip().lower().startswith("https://")


def _parse_iso8601_utc(text: str) -> bool:
    s = str(text).strip()
    if not s:
        return False
    if s.endswith("Z"):
        s = s[:-1] + "+00:00"
    try:
        dt = datetime.fromisoformat(s)
    except ValueError:
        return False
    return dt.tzinfo is not None


def write_templates_if_missing(mu_ref: float, g_min: float, g_max: float) -> None:
    if not TEMPLATE_JSON.exists():
        payload = {
            "mu_ref_gev": mu_ref,
            "g_dimensionless_mu_ref": (g_min + g_max) / 2.0,
            "bridge_observable_origin": "external_dimensionless_observable",
            "notes": "Replace with directly measured external dimensionless observable value.",
        }
        TEMPLATE_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    if not TEMPLATE_META.exists():
        meta = {
            "source_dataset_name": "FILL_REAL_EXTERNAL_DATASET_NAME",
            "citation": "FILL_REAL_CITATION",
            "reference_url": "https://FILL_REAL_URL_OR_DOI",
            "source_version": "FILL_REAL_RELEASE_TAG",
            "provenance_anchor_free": True,
            "seeded_from_registry": False,
            "bridge_observable_origin": "external_dimensionless_observable",
            "acquired_utc": "2026-03-04T00:00:00+00:00",
            "artifact_sha256": "FILL_WITH_SHA256_OF_CANDIDATE_JSON",
            "acquisition_protocol": "external_dimensionless_observable_pipeline_v1",
            "acquisition_command": "FILL_WITH_EXACT_COMMAND_OR_PIPELINE_REF",
            "analyst": "Krzysztof Żuchowski",
            "analyst_email": "krzysztof.zuch@gmail.com",
            "analyst_website": "https://liderpasdom.pl",
            "json_keys": [
                "mu_ref_gev",
                "g_dimensionless_mu_ref",
                "bridge_observable_origin",
                "notes",
            ],
            "notes": "Fill all FILL_* fields before strict run.",
        }
        TEMPLATE_META.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2113 direct dimensionless G pack gate")
    p.add_argument("--candidate-json", default=str(CAND_JSON_DEFAULT), help="Candidate payload JSON.")
    p.add_argument("--candidate-meta", default=str(CAND_META_DEFAULT), help="Candidate metadata JSON.")
    args = p.parse_args()

    cand_json_path = Path(args.candidate_json).resolve()
    cand_meta_path = Path(args.candidate_meta).resolve()

    r2108 = load_json(QW2108_JSON)
    spec = r2108.get("bridge_spec", {})
    contract = r2108.get("strict_contract", {})
    mu_ref_target = float(spec.get("mu_ref_gev", 1.0))
    g_range = spec.get("g_dimensionless_acceptance_range", {})
    g_min = float(g_range.get("min"))
    g_max = float(g_range.get("max"))

    write_templates_if_missing(mu_ref_target, g_min, g_max)

    cand_json_exists = cand_json_path.exists()
    cand_meta_exists = cand_meta_path.exists()
    payload = {}
    meta = {}
    payload_error = ""
    meta_error = ""
    if cand_json_exists:
        try:
            payload = load_json(cand_json_path)
        except Exception as exc:  # noqa: BLE001
            payload_error = str(exc)
    if cand_meta_exists:
        try:
            meta = load_json(cand_meta_path)
        except Exception as exc:  # noqa: BLE001
            meta_error = str(exc)

    mu_ref = payload.get("mu_ref_gev")
    g_dim = payload.get("g_dimensionless_mu_ref")
    origin = str(payload.get("bridge_observable_origin", "")).strip()
    g_si = payload.get("g_si", None)

    flags = {
        "candidate_json_exists": cand_json_exists,
        "candidate_json_parse_ok": bool(cand_json_exists and not payload_error),
        "candidate_meta_exists": cand_meta_exists,
        "candidate_meta_parse_ok": bool(cand_meta_exists and not meta_error),
        "mu_ref_matches_target": bool(
            isinstance(mu_ref, (int, float))
            and math.isfinite(float(mu_ref))
            and abs(float(mu_ref) - mu_ref_target) <= 1e-12
        ),
        "g_dimensionless_positive": bool(
            isinstance(g_dim, (int, float)) and math.isfinite(float(g_dim)) and float(g_dim) > 0.0
        ),
        "g_dimensionless_in_qw2108_range": bool(
            isinstance(g_dim, (int, float)) and g_min <= float(g_dim) <= g_max
        ),
        "origin_external_dimensionless": bool(origin == "external_dimensionless_observable"),
        "g_si_optional_null": bool(g_si is None),
        "meta_source_fields_complete": bool(
            meta
            and all(
                bool(str(meta.get(k, "")).strip())
                for k in ["source_dataset_name", "citation", "reference_url", "source_version"]
            )
        ),
        "meta_reference_url_https": bool(meta and _is_https_url(meta.get("reference_url", ""))),
        "meta_anchor_free_true": bool(meta and meta.get("provenance_anchor_free") is True),
        "meta_seeded_from_registry_false": bool(meta and meta.get("seeded_from_registry") is False),
        "meta_origin_external_dimensionless": bool(
            meta and str(meta.get("bridge_observable_origin", "")).strip() == "external_dimensionless_observable"
        ),
        "meta_acquired_utc_parseable": bool(meta and _parse_iso8601_utc(meta.get("acquired_utc", ""))),
        "meta_json_keys_match": bool(
            meta
            and isinstance(meta.get("json_keys"), list)
            and set(meta.get("json_keys")) == set(payload.keys())
        ),
    }

    all_ready = all(bool(v) for v in flags.values())
    verdict = "GNEWTON_DIRECT_DIMENSIONLESS_PACK_READY" if all_ready else "GNEWTON_DIRECT_DIMENSIONLESS_PACK_PENDING"

    if all_ready:
        READY_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        READY_META.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "candidate_json": str(cand_json_path),
            "candidate_meta": str(cand_meta_path),
            "candidate_template_json": str(TEMPLATE_JSON),
            "candidate_template_meta": str(TEMPLATE_META),
        },
        "errors": {
            "payload_error": payload_error,
            "meta_error": meta_error,
        },
        "qw2108_contract": {
            "mu_ref_target_gev": mu_ref_target,
            "g_range": {"min": g_min, "max": g_max},
            "required_bridge_origin": contract.get("required_bridge_observable_origin"),
            "g_si_optional_must_be_null": contract.get("g_si_input_optional_must_be_null_for_strict"),
        },
        "candidate_summary": {
            "mu_ref_gev": mu_ref,
            "g_dimensionless_mu_ref": g_dim,
            "bridge_observable_origin": origin,
            "g_si": g_si,
        },
        "flags": flags,
        "pass_count": sum(1 for v in flags.values() if bool(v)),
        "total_flags": len(flags),
        "output_artifacts": {
            "ready_payload_json": str(READY_JSON) if all_ready else None,
            "ready_metadata_json": str(READY_META) if all_ready else None,
        },
        "verdict": verdict,
        "required_next_step": (
            "RUN_QW2101_STRICT_DIMENSIONLESS_ONLY_WITH_READY_FILE_THEN_QW2103_QW2092_QW2104_QW2105_QW2094"
            if all_ready
            else "PROVIDE_REAL_DIRECT_DIMENSIONLESS_BRIDGE_CANDIDATE_AND_RERUN_QW2113"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2113: GNEWTON DIRECT DIMENSIONLESS PACK GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{out['pass_count']}/{out['total_flags']}`",
        "",
        "## Candidate Summary",
        f"- mu_ref_gev: `{mu_ref}` (target `{mu_ref_target}`)",
        f"- g_dimensionless_mu_ref: `{g_dim}` (range `[{g_min}, {g_max}]`)",
        f"- bridge_observable_origin: `{origin}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2113] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2113] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2113] verdict={verdict} pass_count={out['pass_count']}/{out['total_flags']}")


if __name__ == "__main__":
    main()
