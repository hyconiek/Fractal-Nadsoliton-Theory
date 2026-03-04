#!/usr/bin/env python3
"""
QW-2050: Freeze bundle for spectral micro bridge external confirmatory run.

Creates a deterministic handoff package with:
- selected kernel and gate results,
- required scripts/reports,
- SHA256 manifest,
- runbook for independent multiteam execution.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
BUNDLE_DIR = ROOT / "external_confirmatory_v2" / "independent_bundle_qw2050_spectral_micro_bridge"
MANIFEST = BUNDLE_DIR / "manifest_qw2050.json"
RUNBOOK = BUNDLE_DIR / "RUNBOOK_QW2050.md"
OUT_JSON = ROOT / "report_qw2050_spectral_micro_bridge_freeze_bundle.json"
OUT_MD = ROOT / "RAPORT_QW2050_SPECTRAL_MICRO_BRIDGE_FREEZE_BUNDLE.md"
DATA_SOURCES_DOC = "DATA_SOURCES_EXTERNAL_DOWNLOADS.md"


REQUIRED_FILES = [
    "QW_1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.py",
    "QW_2038_DERIVATION_COMPATIBLE_KERNEL_REFREEZE_SCAN.py",
    "QW_2039_DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE.py",
    "QW_2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.py",
    "QW_2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.py",
    "report_qw2038_derivation_compatible_kernel_refreeze_scan.json",
    "report_qw2039_derivation_compatible_refrozen_kernel_gate.json",
    "report_qw2048_spectral_phase_locked_pointwise_derivation.json",
    "report_qw2049_spectral_micro_stagec_intersection_gate.json",
    "RAPORT_QW2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.md",
    "RAPORT_QW2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.md",
    "external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv",
    "external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg/pta_v2_pairs.csv",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    d2049 = json.loads((ROOT / "report_qw2049_spectral_micro_stagec_intersection_gate.json").read_text(encoding="utf-8"))

    BUNDLE_DIR.mkdir(parents=True, exist_ok=True)

    missing = []
    manifest_rows: List[Dict[str, object]] = []
    for rel in REQUIRED_FILES:
        p = ROOT / rel
        if not p.exists():
            missing.append(rel)
            continue
        manifest_rows.append(
            {
                "path": rel,
                "bytes": int(p.stat().st_size),
                "sha256": sha256_file(p),
            }
        )

    selected_kernel = d2049["stagec_pool"]["selected_kernel"]

    bundle_meta = {
        "bundle_id": "QW2050_SPECTRAL_MICRO_BRIDGE_FREEZE",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_gate": "report_qw2049_spectral_micro_stagec_intersection_gate.json",
        "source_verdict": d2049.get("verdict"),
        "source_readiness": d2049.get("readiness"),
        "external_sources_document": DATA_SOURCES_DOC,
        "selected_kernel": selected_kernel,
        "required_files_count": len(REQUIRED_FILES),
        "present_files_count": len(manifest_rows),
        "missing_files": missing,
        "manifest": manifest_rows,
    }
    MANIFEST.write_text(json.dumps(bundle_meta, ensure_ascii=False, indent=2), encoding="utf-8")

    runbook_lines = [
        "# RUNBOOK QW-2050: Spectral Micro Bridge Independent Confirmatory",
        "",
        f"- Bundle generated UTC: {bundle_meta['generated_utc']}",
        f"- Source verdict: {bundle_meta['source_verdict']}",
        f"- Source readiness: {bundle_meta['source_readiness']}",
        "",
        "## Fixed Kernel",
        (
            "- omega/phi/beta/eta: "
            f"{selected_kernel['omega']:.6f} / {selected_kernel['phi']:.6f} / "
            f"{selected_kernel['beta']:.6f} / {selected_kernel['eta']:.6f}"
        ),
        "",
        "## Integrity Step",
        "1. Verify all SHA256 entries from `manifest_qw2050.json`.",
        "2. Reject bundle if any hash mismatch is found.",
        "",
        "## External Data Sources (Not Frozen In Git)",
        "- Large public raw archives are intentionally not included in this bundle.",
        f"- Source list and download commands: `{DATA_SOURCES_DOC}`",
        "- Reproducibility is based on fixed scripts/manifests plus external source provenance.",
        "",
        "## Independent Execution",
        "1. Run `python3 QW_2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.py`.",
        "2. Run `python3 QW_2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.py`.",
        "3. Confirm `SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS`.",
        "",
        "## Decision Rule",
        "- PASS only if all QW-2049 flags are True and external primary/stress criteria pass.",
        "- Any change in kernel vector invalidates this bundle and requires new freeze.",
        "",
        "## Notes",
        "- No sector retune is allowed.",
        "- Use only files listed in manifest.",
    ]
    RUNBOOK.write_text("\n".join(runbook_lines) + "\n", encoding="utf-8")

    flags = {
        "source_gate_pass": bool(d2049.get("verdict") == "SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS"),
        "all_required_files_present": bool(len(missing) == 0),
        "data_sources_doc_present": bool((ROOT / DATA_SOURCES_DOC).exists()),
        "manifest_written": bool(MANIFEST.exists()),
        "runbook_written": bool(RUNBOOK.exists()),
    }
    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    if pass_count == total_flags:
        verdict = "SPECTRAL_MICRO_BRIDGE_FREEZE_BUNDLE_READY"
        readiness = "READY_FOR_TRUE_INDEPENDENT_MULTITEAM_RUN"
    elif pass_count >= 3:
        verdict = "SPECTRAL_MICRO_BRIDGE_FREEZE_BUNDLE_PARTIAL"
        readiness = "PARTIAL"
    else:
        verdict = "SPECTRAL_MICRO_BRIDGE_FREEZE_BUNDLE_FAIL"
        readiness = "NOT_READY"

    out = {
        "generated_utc": bundle_meta["generated_utc"],
        "bundle_dir": str(BUNDLE_DIR.relative_to(ROOT)),
        "manifest": str(MANIFEST.relative_to(ROOT)),
        "runbook": str(RUNBOOK.relative_to(ROOT)),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "missing_files": missing,
        "selected_kernel": selected_kernel,
        "verdict": verdict,
        "readiness": readiness,
        "required_next_step": "RUN_TRULY_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2050: SPECTRAL MICRO BRIDGE FREEZE BUNDLE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Bundle",
        f"- dir: `{out['bundle_dir']}`",
        f"- manifest: `{out['manifest']}`",
        f"- runbook: `{out['runbook']}`",
        "",
        "## Fixed Kernel",
        (
            f"- omega/phi/beta/eta: {selected_kernel['omega']:.6f} / {selected_kernel['phi']:.6f} / "
            f"{selected_kernel['beta']:.6f} / {selected_kernel['eta']:.6f}"
        ),
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    if missing:
        lines.append("")
        lines.append("## Missing Files")
        for m in missing:
            lines.append(f"- {m}")

    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {out['required_next_step']}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2050] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2050] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2050] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
