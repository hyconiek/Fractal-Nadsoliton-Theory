#!/usr/bin/env python3
"""
QW-2109: strict external evidence-manifest gate for T3/T4 inputs.

Purpose:
- enforce auditable evidence manifests for external artifacts,
- verify sidecar metadata integrity against on-disk payloads,
- reduce false-success risk from weak or unverifiable provenance.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple


ROOT = Path(__file__).resolve().parent

DEFAULT_HZ_CSV = ROOT / "external_hz_nodes_qw2099.csv"
DEFAULT_HZ_META = ROOT / "external_hz_nodes_qw2099.metadata.json"
DEFAULT_G_JSON = ROOT / "external_gnewton_bridge_qw2101.json"
DEFAULT_G_META = ROOT / "external_gnewton_bridge_qw2101.metadata.json"

OUT_JSON = ROOT / "report_qw2109_strict_external_evidence_manifest_gate.json"
OUT_MD = ROOT / "RAPORT_QW2109_STRICT_EXTERNAL_EVIDENCE_MANIFEST_GATE.md"


def _has_value(meta: Dict, key: str) -> bool:
    return bool(str(meta.get(key, "")).strip())


def _is_https_url(url: str) -> bool:
    u = str(url).strip().lower()
    return u.startswith("https://")


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


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _metadata_core_complete(meta: Dict) -> bool:
    required = [
        "source_dataset_name",
        "citation",
        "reference_url",
        "source_version",
    ]
    if not all(_has_value(meta, k) for k in required):
        return False
    low = " ".join(str(meta.get(k, "")).lower() for k in required)
    return "placeholder" not in low


def _metadata_evidence_complete(meta: Dict) -> bool:
    required = [
        "acquired_utc",
        "artifact_sha256",
        "acquisition_protocol",
        "acquisition_command",
        "analyst",
    ]
    return all(_has_value(meta, k) for k in required)


def _read_hz_csv(path: Path) -> Tuple[List[str], int]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        header = list(r.fieldnames or [])
        n = 0
        for _ in r:
            n += 1
    return header, n


def _g_payload_has_numeric_bridge(g_data: Dict) -> bool:
    mu = g_data.get("mu_ref_gev")
    g_dim = g_data.get("g_dimensionless_mu_ref")
    g_si = g_data.get("g_si")
    mu_ok = isinstance(mu, (int, float)) and math.isfinite(float(mu)) and float(mu) > 0.0
    dim_ok = isinstance(g_dim, (int, float)) and math.isfinite(float(g_dim)) and float(g_dim) > 0.0
    si_ok = isinstance(g_si, (int, float)) and math.isfinite(float(g_si)) and float(g_si) > 0.0
    return bool(mu_ok and (dim_ok or si_ok))


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2109 strict external evidence-manifest gate")
    p.add_argument("--hz-csv", default=str(DEFAULT_HZ_CSV), help="H(z) CSV artifact.")
    p.add_argument("--hz-meta", default=str(DEFAULT_HZ_META), help="H(z) metadata sidecar.")
    p.add_argument("--g-json", default=str(DEFAULT_G_JSON), help="G bridge JSON artifact.")
    p.add_argument("--g-meta", default=str(DEFAULT_G_META), help="G bridge metadata sidecar.")
    args = p.parse_args()

    hz_csv = Path(args.hz_csv).resolve()
    hz_meta_path = Path(args.hz_meta).resolve()
    g_json_path = Path(args.g_json).resolve()
    g_meta_path = Path(args.g_meta).resolve()

    hz_file_exists = hz_csv.exists()
    hz_meta_exists = hz_meta_path.exists()
    g_file_exists = g_json_path.exists()
    g_meta_exists = g_meta_path.exists()

    hz_meta: Dict = {}
    g_meta: Dict = {}
    g_data: Dict = {}
    hz_header: List[str] = []
    hz_rows = 0
    hz_error = ""
    g_error = ""
    hz_meta_error = ""
    g_meta_error = ""

    if hz_file_exists:
        try:
            hz_header, hz_rows = _read_hz_csv(hz_csv)
        except Exception as exc:  # noqa: BLE001
            hz_error = str(exc)

    if hz_meta_exists:
        try:
            hz_meta = _load_json(hz_meta_path)
        except Exception as exc:  # noqa: BLE001
            hz_meta_error = str(exc)

    if g_file_exists:
        try:
            g_data = _load_json(g_json_path)
        except Exception as exc:  # noqa: BLE001
            g_error = str(exc)

    if g_meta_exists:
        try:
            g_meta = _load_json(g_meta_path)
        except Exception as exc:  # noqa: BLE001
            g_meta_error = str(exc)

    hz_sha_actual = _sha256_file(hz_csv) if hz_file_exists and not hz_error else ""
    g_sha_actual = _sha256_file(g_json_path) if g_file_exists and not g_error else ""
    hz_sha_declared = str(hz_meta.get("artifact_sha256", "")).strip() if hz_meta else ""
    g_sha_declared = str(g_meta.get("artifact_sha256", "")).strip() if g_meta else ""

    expected_hz_cols = {"z", "h_km_s_mpc", "sigma_total"}
    hz_cols_meta = hz_meta.get("columns_schema", [])
    hz_rows_declared = hz_meta.get("record_count")

    g_keys_meta = g_meta.get("json_keys", [])
    g_origin_meta = str(g_meta.get("bridge_observable_origin", "")).strip()

    hz_flags = {
        "hz_file_exists": hz_file_exists,
        "hz_file_parse_ok": bool(hz_file_exists and not hz_error),
        "hz_meta_exists": hz_meta_exists,
        "hz_meta_parse_ok": bool(hz_meta_exists and not hz_meta_error),
        "hz_meta_core_complete": bool(hz_meta and _metadata_core_complete(hz_meta)),
        "hz_meta_evidence_complete": bool(hz_meta and _metadata_evidence_complete(hz_meta)),
        "hz_reference_url_https": bool(hz_meta and _is_https_url(hz_meta.get("reference_url", ""))),
        "hz_acquired_utc_parseable": bool(hz_meta and _parse_iso8601_utc(hz_meta.get("acquired_utc", ""))),
        "hz_sha_declared_present": bool(hz_sha_declared),
        "hz_sha_matches_payload": bool(hz_sha_declared and hz_sha_actual and hz_sha_declared == hz_sha_actual),
        "hz_seeded_from_registry_false": bool(hz_meta and (hz_meta.get("seeded_from_registry") is False)),
        "hz_schema_has_required_columns": bool(set(hz_header) >= expected_hz_cols),
        "hz_rows_positive": bool(hz_rows > 0),
        "hz_meta_record_count_matches": bool(
            isinstance(hz_rows_declared, int) and hz_rows_declared == hz_rows
        ),
        "hz_meta_columns_schema_matches": bool(
            isinstance(hz_cols_meta, list) and set(hz_cols_meta) == set(hz_header)
        ),
    }

    g_flags = {
        "g_file_exists": g_file_exists,
        "g_file_parse_ok": bool(g_file_exists and not g_error),
        "g_meta_exists": g_meta_exists,
        "g_meta_parse_ok": bool(g_meta_exists and not g_meta_error),
        "g_meta_core_complete": bool(g_meta and _metadata_core_complete(g_meta)),
        "g_meta_evidence_complete": bool(g_meta and _metadata_evidence_complete(g_meta)),
        "g_reference_url_https": bool(g_meta and _is_https_url(g_meta.get("reference_url", ""))),
        "g_acquired_utc_parseable": bool(g_meta and _parse_iso8601_utc(g_meta.get("acquired_utc", ""))),
        "g_sha_declared_present": bool(g_sha_declared),
        "g_sha_matches_payload": bool(g_sha_declared and g_sha_actual and g_sha_declared == g_sha_actual),
        "g_seeded_from_registry_false": bool(g_meta and (g_meta.get("seeded_from_registry") is False)),
        "g_payload_has_numeric_bridge": bool(g_data and _g_payload_has_numeric_bridge(g_data)),
        "g_meta_json_keys_match": bool(
            isinstance(g_keys_meta, list) and set(g_keys_meta) == set(g_data.keys())
        ),
        "g_origin_declared_nonempty": bool(g_origin_meta),
    }

    all_flags = {**hz_flags, **g_flags}
    pass_count = sum(1 for v in all_flags.values() if bool(v))
    total_flags = len(all_flags)
    all_pass = all(bool(v) for v in all_flags.values())

    verdict = (
        "STRICT_EXTERNAL_EVIDENCE_MANIFEST_GATE_PASS"
        if all_pass
        else "STRICT_EXTERNAL_EVIDENCE_MANIFEST_GATE_PENDING"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "hz_csv": str(hz_csv),
            "hz_meta": str(hz_meta_path),
            "g_json": str(g_json_path),
            "g_meta": str(g_meta_path),
        },
        "errors": {
            "hz_error": hz_error,
            "hz_meta_error": hz_meta_error,
            "g_error": g_error,
            "g_meta_error": g_meta_error,
        },
        "hz_manifest": {
            "sha256_declared": hz_sha_declared,
            "sha256_actual": hz_sha_actual,
            "header": hz_header,
            "rows": hz_rows,
            "metadata_record_count": hz_rows_declared,
            "metadata_columns_schema": hz_cols_meta,
        },
        "g_manifest": {
            "sha256_declared": g_sha_declared,
            "sha256_actual": g_sha_actual,
            "json_keys_actual": sorted(g_data.keys()) if isinstance(g_data, dict) else [],
            "metadata_json_keys": g_keys_meta,
            "metadata_bridge_observable_origin": g_origin_meta,
        },
        "hz_flags": hz_flags,
        "g_flags": g_flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "EVIDENCE_MANIFEST_READY_RERUN_QW2106_QW2099_QW2101_QW2102_QW2103_QW2104_QW2094"
            if all_pass
            else "FILL_STRICT_EVIDENCE_SIDECARS_SHA_ACQUIRED_UTC_PROTOCOL_AND_SCHEMA_THEN_RERUN_QW2109"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2109: STRICT EXTERNAL EVIDENCE MANIFEST GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## H(z) Manifest",
        f"- sha256 declared: `{hz_sha_declared}`",
        f"- sha256 actual:   `{hz_sha_actual}`",
        f"- rows: `{hz_rows}`",
        f"- header: `{hz_header}`",
        "",
        "## G-bridge Manifest",
        f"- sha256 declared: `{g_sha_declared}`",
        f"- sha256 actual:   `{g_sha_actual}`",
        f"- json_keys_actual: `{sorted(g_data.keys()) if isinstance(g_data, dict) else []}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2109] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2109] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2109] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
