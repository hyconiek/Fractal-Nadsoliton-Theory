#!/usr/bin/env python3
"""
QW-2053: Independent multiteam confirmatory protocol lock.

Purpose:
- formalize a single immutable confirmatory contract for external teams,
- lock kernel + pass/fail criteria + artifact hashes,
- verify that prerequisite internal gates are already green.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent

LOCK_DIR = ROOT / "external_confirmatory_v2" / "independent_multiteam_lock_qw2053"
LOCK_JSON = LOCK_DIR / "protocol_lock_qw2053.json"
RUNBOOK_MD = LOCK_DIR / "RUNBOOK_QW2053.md"

OUT_JSON = ROOT / "report_qw2053_independent_multiteam_protocol_lock.json"
OUT_MD = ROOT / "RAPORT_QW2053_INDEPENDENT_MULTITEAM_PROTOCOL_LOCK.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def canonical_sha256(data: Dict) -> str:
    b = json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(b).hexdigest()


def main() -> None:
    LOCK_DIR.mkdir(parents=True, exist_ok=True)

    p2049 = ROOT / "report_qw2049_spectral_micro_stagec_intersection_gate.json"
    p2051 = ROOT / "report_qw2051_independent_rehearsal_gate.json"
    p2052 = ROOT / "report_qw2052_external_source_only_governance_gate.json"
    p2050_manifest = ROOT / "external_confirmatory_v2" / "independent_bundle_qw2050_spectral_micro_bridge" / "manifest_qw2050.json"
    p2050_runbook = ROOT / "external_confirmatory_v2" / "independent_bundle_qw2050_spectral_micro_bridge" / "RUNBOOK_QW2050.md"
    p_sources = ROOT / "DATA_SOURCES_EXTERNAL_DOWNLOADS.md"

    required_paths: List[Path] = [
        p2049,
        p2051,
        p2052,
        p2050_manifest,
        p2050_runbook,
        p_sources,
        ROOT / "QW_2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.py",
        ROOT / "QW_2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.py",
        ROOT / "QW_2051_INDEPENDENT_REHEARSAL_GATE.py",
        ROOT / "QW_2052_EXTERNAL_SOURCE_ONLY_GOVERNANCE_GATE.py",
    ]

    missing = [str(p.relative_to(ROOT)) for p in required_paths if not p.exists()]

    d2049 = load_json(p2049) if p2049.exists() else {}
    d2051 = load_json(p2051) if p2051.exists() else {}
    d2052 = load_json(p2052) if p2052.exists() else {}

    kernel = d2049.get("stagec_pool", {}).get("selected_kernel", {})
    ext_primary = d2049.get("external", {}).get("primary", {})
    ext_stress = d2049.get("external", {}).get("stress", {})
    gate_flags = d2049.get("flags", {})

    contract = {
        "lock_id": "QW2053_INDEPENDENT_MULTITEAM_PROTOCOL_LOCK",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "TRULY_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE",
        "fixed_kernel": {
            "omega": float(kernel.get("omega", 0.0)),
            "phi": float(kernel.get("phi", 0.0)),
            "beta": float(kernel.get("beta", 0.0)),
            "eta": float(kernel.get("eta", 0.0)),
        },
        "hard_no_retune_rules": {
            "forbid_sector_retune": True,
            "forbid_kernel_change": True,
            "forbid_threshold_change": True,
            "forbid_posthoc_selection": True,
            "forbid_internal_proxy_substitution": True,
        },
        "hard_acceptance_criteria": {
            "require_qw2049_verdict": "SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS",
            "require_qw2051_verdict": "INDEPENDENT_REHEARSAL_GATE_PASS",
            "require_qw2052_verdict": "EXTERNAL_SOURCE_ONLY_GOVERNANCE_PASS",
            "require_all_qw2049_flags_true": True,
            "external_primary": {
                "pearson_corr_min": float(ext_primary.get("permutation", {}).get("q95_corr", 0.0)),
                "rmse_gain_ratio_min": 0.0,
                "p_corr_max": 0.01,
                "p_rmse_gain_max": 0.01,
            },
            "external_stress": {
                "pearson_corr_min": float(ext_stress.get("permutation", {}).get("q95_corr", 0.0)),
                "rmse_gain_ratio_min": 0.0,
                "p_corr_max": 0.01,
                "p_rmse_gain_max": 0.01,
            },
        },
        "required_execution_order": [
            "python3 QW_2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.py",
            "python3 QW_2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.py",
            "python3 QW_2051_INDEPENDENT_REHEARSAL_GATE.py",
            "python3 QW_2052_EXTERNAL_SOURCE_ONLY_GOVERNANCE_GATE.py",
        ],
        "required_artifacts": [
            "external_confirmatory_v2/independent_bundle_qw2050_spectral_micro_bridge/manifest_qw2050.json",
            "external_confirmatory_v2/independent_bundle_qw2050_spectral_micro_bridge/RUNBOOK_QW2050.md",
            "DATA_SOURCES_EXTERNAL_DOWNLOADS.md",
            "report_qw2049_spectral_micro_stagec_intersection_gate.json",
            "report_qw2051_independent_rehearsal_gate.json",
            "report_qw2052_external_source_only_governance_gate.json",
        ],
    }

    artifact_hashes = []
    for p in required_paths:
        if p.exists():
            artifact_hashes.append(
                {
                    "path": str(p.relative_to(ROOT)),
                    "bytes": int(p.stat().st_size),
                    "sha256": sha256_file(p),
                }
            )

    lock_payload = {
        "contract": contract,
        "artifact_hashes": artifact_hashes,
        "source_reports": {
            "qw2049": "report_qw2049_spectral_micro_stagec_intersection_gate.json",
            "qw2051": "report_qw2051_independent_rehearsal_gate.json",
            "qw2052": "report_qw2052_external_source_only_governance_gate.json",
        },
    }
    lock_payload["lock_sha256"] = canonical_sha256(lock_payload)

    LOCK_JSON.write_text(json.dumps(lock_payload, ensure_ascii=False, indent=2), encoding="utf-8")

    runbook_lines = [
        "# RUNBOOK QW-2053: Independent Multiteam Protocol Lock",
        "",
        f"- Lock file: `{LOCK_JSON.relative_to(ROOT)}`",
        f"- lock_sha256: `{lock_payload['lock_sha256']}`",
        "",
        "## Non-Negotiable Rules",
        "- No kernel change and no sector retune.",
        "- No threshold edits and no post-hoc selection.",
        "- No internal proxy substitution for external confirmatory inputs.",
        "",
        "## Required Execution Order",
    ]
    for step in contract["required_execution_order"]:
        runbook_lines.append(f"- `{step}`")
    runbook_lines.extend(
        [
            "",
            "## Required Source/Bundle Material",
            "- `external_confirmatory_v2/independent_bundle_qw2050_spectral_micro_bridge/manifest_qw2050.json`",
            "- `external_confirmatory_v2/independent_bundle_qw2050_spectral_micro_bridge/RUNBOOK_QW2050.md`",
            "- `DATA_SOURCES_EXTERNAL_DOWNLOADS.md`",
            "",
            "## Pass Condition",
            "- Final pass only if QW-2049, QW-2051 and QW-2052 verdicts match locked contract and all hard flags remain true.",
        ]
    )
    RUNBOOK_MD.write_text("\n".join(runbook_lines) + "\n", encoding="utf-8")

    flags = {
        "required_artifacts_present": len(missing) == 0,
        "spectral_gate_pass": d2049.get("verdict") == "SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS",
        "rehearsal_pass": d2051.get("verdict") == "INDEPENDENT_REHEARSAL_GATE_PASS",
        "source_only_governance_pass": d2052.get("verdict") == "EXTERNAL_SOURCE_ONLY_GOVERNANCE_PASS",
        "qw2049_all_hard_flags_true": bool(gate_flags) and all(bool(v) for v in gate_flags.values()),
        "lock_file_written": LOCK_JSON.exists(),
        "runbook_written": RUNBOOK_MD.exists(),
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    if pass_count == total_flags:
        verdict = "INDEPENDENT_MULTITEAM_PROTOCOL_LOCK_READY"
        readiness = "READY_FOR_TRUE_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE"
    elif pass_count >= total_flags - 2:
        verdict = "INDEPENDENT_MULTITEAM_PROTOCOL_LOCK_PARTIAL"
        readiness = "PARTIAL"
    else:
        verdict = "INDEPENDENT_MULTITEAM_PROTOCOL_LOCK_FAIL"
        readiness = "NOT_READY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "lock_file": str(LOCK_JSON.relative_to(ROOT)),
        "lock_sha256": lock_payload["lock_sha256"],
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "missing_artifacts": missing,
        "verdict": verdict,
        "readiness": readiness,
        "required_next_step": "HAND_OFF_QW2050_QW2053_LOCKED_PACKAGE_TO_TRULY_INDEPENDENT_MULTITEAM",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2053: INDEPENDENT MULTITEAM PROTOCOL LOCK",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- lock_sha256: `{out['lock_sha256']}`",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    if missing:
        lines.append("")
        lines.append("## Missing Artifacts")
        for m in missing:
            lines.append(f"- {m}")
    lines.extend(
        [
            "",
            "## Artifacts",
            f"- lock: `{LOCK_JSON.name}`",
            f"- runbook: `{RUNBOOK_MD.name}`",
            f"- json: `{OUT_JSON.name}`",
            "",
            "## Required Next Step",
            f"- {out['required_next_step']}",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2053] Saved lock: {LOCK_JSON.relative_to(ROOT)}")
    print(f"[QW-2053] Saved runbook: {RUNBOOK_MD.relative_to(ROOT)}")
    print(f"[QW-2053] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2053] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2053] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
