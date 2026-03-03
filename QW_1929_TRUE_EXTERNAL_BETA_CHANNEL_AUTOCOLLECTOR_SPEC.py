#!/usr/bin/env python3
"""
QW-1929: True external beta-channel autocollector spec.

Creates an executable collection spec for a genuinely external beta-channel
package built directly from public raw archives (no internal proxy tables).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1929_true_external_beta_channel_autocollector_spec.json"
OUT_MD = ROOT / "RAPORT_QW1929_TRUE_EXTERNAL_BETA_CHANNEL_AUTOCOLLECTOR_SPEC.md"
SPEC_MD = ROOT / "external_confirmatory_v2" / "COLLECTION_SPEC_QW1929_TRUE_EXTERNAL_AUTOCOLLECTOR.md"


def main() -> None:
    d1925 = json.loads((ROOT / "report_qw1925_true_external_beta_channel_collection_spec.json").read_text(encoding="utf-8"))
    q1925 = d1925["spec"]

    spec = {
        "objective": (
            "Build beta-channel confirmatory package from independent external public datasets "
            "without using internal proxy/rebuild pair tables."
        ),
        "data_sources": {
            "nanograv_15yr_timing_archive": {
                "required_local_path": "NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz",
                "archive_contents_required": [
                    "NANOGrav15yr_PulsarTiming_v2.1.0/narrowband/tim/*.nb.tim",
                    "NANOGrav15yr_PulsarTiming_v2.1.0/narrowband/README.narrowband",
                ],
                "provenance_class": "independent_external_public_dataset",
            },
            "gwosc_event_catalog": {
                "api_url": "https://www.gw-openscience.org/eventapi/json/GWTC/",
                "cache_path": "external_confirmatory_v2/beta_channel_true_external/sources/gwosc_gwtc_eventapi.json",
                "required_fields": ["events", "GPS"],
                "provenance_class": "independent_external_public_dataset",
            },
        },
        "collector_rules": {
            "forbid_internal_pair_tables": True,
            "forbidden_inputs": [
                "external_confirmatory_v2/confirmatory_dataset_internal_proxy*/pta_v2_pairs.csv",
                "external_confirmatory_v2/confirmatory_dataset_external_source_rebuild*/pta_v2_pairs.csv",
            ],
            "tim_parser": {
                "format": "tempo2",
                "required_numeric_fields": ["freq_mhz", "mjd", "toa_err_us"],
                "optional_flag_fields": ["snr", "gof", "flux", "tobs", "chan", "subint", "pn", "pta"],
            },
            "quality_filters": {
                "pta_flag_equals_nanograv_if_present": True,
                "min_snr": 8.0,
                "max_toa_err_us": 12.0,
                "rows_target_after_subsample": 4000,
                "rows_min_hard_gate": 2 * int(q1925["power_targets"]["n_holdout_min_for_power_0p90"]),
            },
            "determinism": {
                "row_sampling": "min sha256(pair_id) across full parse",
                "split": "sha256(pair_id) parity",
                "intervention_map": "nearest GW event in MJD; regime=pre/post by signed delta",
            },
        },
        "feature_construction": {
            "theta_deg": "linear map of freq_mhz to [0,180]",
            "hxy": "robustly normalized weighted score from snr, toa_err, gof, flux, tobs",
            "group_scope": "obs_key=pulsar|pn|round(mjd,3)",
            "f_std": "std(hxy) within obs_key sequence",
            "f_autoc1": "lag-1 autocorrelation within obs_key sequence",
            "f_switch": "sign-switch rate of first difference within obs_key sequence",
            "f_slope": "linear slope dhxy/dtheta within obs_key sequence",
        },
        "output_package": {
            "base_dir": "external_confirmatory_v2/beta_channel_true_external",
            "required_files": list(q1925["hard_requirements"]["required_files"]),
            "required_roles_in_manifest": list(q1925["hard_requirements"]["required_roles_in_manifest"]),
            "required_pair_columns": list(q1925["dataset_schema"]["beta_channel_pairs_required_columns"]),
            "required_event_columns": list(q1925["dataset_schema"]["intervention_events_required_columns"]),
            "regime_allowed_values": list(q1925["dataset_schema"]["regime_allowed_values"]),
        },
        "externality_statement_requirements": {
            "must_include_tokens": list(q1925["hard_requirements"]["required_externality_statement_tokens"]),
            "must_not_include_phrases": ["internal proxy", "not independent"],
            "provider_forbidden_tokens": list(q1925["hard_requirements"]["forbidden_provider_tokens"]),
            "minimum_truth_claim": (
                "package assembled from independent external raw datasets and event catalogs, "
                "without internal proxy/rebuild pair tables"
            ),
        },
        "acceptance_gate": {
            "script": "QW_1927_TRUE_EXTERNAL_BETA_CHANNEL_READINESS_GATE.py",
            "required_flags_true": [
                "externality_ok",
                "schema_pairs_ok",
                "schema_events_ok",
                "regime_values_ok",
                "has_pre_and_post",
                "n_pairs_ge_power_target",
            ],
            "required_verdict": "TRUE_EXTERNAL_BETA_CHANNEL_READY",
        },
    }

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "based_on": [
            "report_qw1925_true_external_beta_channel_collection_spec.json",
            "QW_1927_TRUE_EXTERNAL_BETA_CHANNEL_READINESS_GATE.py",
        ],
        "spec": spec,
        "verdict": "TRUE_EXTERNAL_AUTOCOLLECTOR_SPEC_READY",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1929: TRUE EXTERNAL BETA-CHANNEL AUTOCOLLECTOR SPEC",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Core Requirement",
        "- Build package directly from public external raw data (NANOGrav archive + GWOSC event API),",
        "- without using internal proxy/rebuild pair tables.",
        "",
        "## Hard Acceptance",
        "- Run `QW_1927_TRUE_EXTERNAL_BETA_CHANNEL_READINESS_GATE.py`.",
        "- Required verdict: `TRUE_EXTERNAL_BETA_CHANNEL_READY`.",
        "- All hard flags must be `True` (including `externality_ok`).",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
        f"- Spec MD: `{SPEC_MD.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    SPEC_MD.parent.mkdir(parents=True, exist_ok=True)
    spec_lines = [
        "# COLLECTION SPEC QW-1929: TRUE EXTERNAL AUTOCOLLECTOR",
        "",
        "## Scope",
        "Build `external_confirmatory_v2/beta_channel_true_external` from:",
        "1. `NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz` raw timing archive,",
        "2. GWOSC GWTC event catalog (`https://www.gw-openscience.org/eventapi/json/GWTC/`).",
        "",
        "## Non-Negotiable Rule",
        "Do NOT read any internal proxy/rebuild pair tables.",
        "",
        "## Required Outputs",
        "- `manifest_beta_channel.json`",
        "- `beta_channel_pairs.csv`",
        "- `intervention_events.csv`",
        "- `protocol_freeze.json`",
        "",
        "## Pair Schema",
        "- pair_id, theta_deg, hxy, f_std, f_autoc1, f_switch, f_slope, intervention_id, regime",
        "",
        "## Externality Guardrail",
        "- provider must not contain INTERNAL / INTERNAL_PROXY,",
        "- statement must include `independent` and `external`,",
        "- statement must not include `not independent` or `internal proxy`.",
        "",
        "## Gate",
        "- run `python3 QW_1927_TRUE_EXTERNAL_BETA_CHANNEL_READINESS_GATE.py`",
        "- expected verdict: `TRUE_EXTERNAL_BETA_CHANNEL_READY`",
        "",
    ]
    SPEC_MD.write_text("\n".join(spec_lines), encoding="utf-8")

    print(f"[QW-1929] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1929] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1929] Saved spec: {SPEC_MD.name}")


if __name__ == "__main__":
    main()

