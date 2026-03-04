#!/usr/bin/env python3
"""
QW-2051: Independent rehearsal gate for QW-2050 freeze bundle.

Performs a local dry-run of what an external team would do, but in an isolated
/tmp workspace:
1) verify manifest hashes,
2) rerun QW-2048 and QW-2049 in isolation,
3) compare core outputs against baseline reports in repository.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
BUNDLE_DIR = ROOT / "external_confirmatory_v2" / "independent_bundle_qw2050_spectral_micro_bridge"
MANIFEST_PATH = BUNDLE_DIR / "manifest_qw2050.json"
OUT_JSON = ROOT / "report_qw2051_independent_rehearsal_gate.json"
OUT_MD = ROOT / "RAPORT_QW2051_INDEPENDENT_REHEARSAL_GATE.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_tree_files(manifest: Dict, target_root: Path) -> None:
    for row in manifest["manifest"]:
        rel = Path(row["path"])
        src = ROOT / rel
        dst = target_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def run_py(script_rel: str, cwd: Path) -> Dict[str, object]:
    p = subprocess.run(
        ["python3", script_rel],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "returncode": int(p.returncode),
        "stdout": p.stdout[-6000:],
        "stderr": p.stderr[-6000:],
    }


def kernel_equal(k1: Dict[str, float], k2: Dict[str, float], tol: float = 1e-12) -> bool:
    keys = ["omega", "phi", "beta", "eta"]
    for k in keys:
        if abs(float(k1[k]) - float(k2[k])) > tol:
            return False
    return True


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    # 1) Hash verification against current repo files.
    hash_mismatches = []
    for row in manifest["manifest"]:
        p = ROOT / row["path"]
        got = sha256_file(p)
        if got != row["sha256"]:
            hash_mismatches.append(
                {
                    "path": row["path"],
                    "expected": row["sha256"],
                    "got": got,
                }
            )

    # 2) Isolated rerun.
    with tempfile.TemporaryDirectory(prefix="qw2051_rehearsal_") as td:
        tmp = Path(td)
        copy_tree_files(manifest, tmp)

        r2048 = run_py("QW_2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.py", cwd=tmp)
        r2049 = run_py("QW_2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.py", cwd=tmp)

        rerun2048 = None
        rerun2049 = None
        if (tmp / "report_qw2048_spectral_phase_locked_pointwise_derivation.json").exists():
            rerun2048 = json.loads(
                (tmp / "report_qw2048_spectral_phase_locked_pointwise_derivation.json").read_text(encoding="utf-8")
            )
        if (tmp / "report_qw2049_spectral_micro_stagec_intersection_gate.json").exists():
            rerun2049 = json.loads(
                (tmp / "report_qw2049_spectral_micro_stagec_intersection_gate.json").read_text(encoding="utf-8")
            )

        # 3) Compare with baseline in repo.
        base2048 = json.loads((ROOT / "report_qw2048_spectral_phase_locked_pointwise_derivation.json").read_text(encoding="utf-8"))
        base2049 = json.loads((ROOT / "report_qw2049_spectral_micro_stagec_intersection_gate.json").read_text(encoding="utf-8"))

        metrics_stable = False
        kernel_stable = False
        if rerun2049 is not None:
            rb = rerun2049["external"]
            bb = base2049["external"]
            keys = [
                ("primary", "holdout_metrics", "pearson_corr"),
                ("primary", "holdout_metrics", "rmse_gain_ratio"),
                ("stress", "holdout_metrics", "pearson_corr"),
                ("stress", "holdout_metrics", "rmse_gain_ratio"),
            ]
            deltas = []
            for a, b, c in keys:
                deltas.append(abs(float(rb[a][b][c]) - float(bb[a][b][c])))
            metrics_stable = bool(max(deltas) <= 1e-12)

            kernel_stable = kernel_equal(
                rerun2049["stagec_pool"]["selected_kernel"],
                base2049["stagec_pool"]["selected_kernel"],
                tol=1e-12,
            )

        flags = {
            "manifest_hashes_match": bool(len(hash_mismatches) == 0),
            "rerun_qw2048_exit0": bool(r2048["returncode"] == 0),
            "rerun_qw2048_pass": bool(rerun2048 is not None and rerun2048.get("verdict") == "SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION_PASS"),
            "rerun_qw2049_exit0": bool(r2049["returncode"] == 0),
            "rerun_qw2049_pass": bool(rerun2049 is not None and rerun2049.get("verdict") == "SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS"),
            "selected_kernel_stable": bool(kernel_stable),
            "key_metrics_stable": bool(metrics_stable),
        }

        pass_count = int(sum(1 for v in flags.values() if v))
        total_flags = int(len(flags))

        if pass_count == total_flags:
            verdict = "INDEPENDENT_REHEARSAL_GATE_PASS"
            readiness = "REHEARSAL_CONFIRMS_BUNDLE_REPRODUCIBILITY"
        elif pass_count >= 5:
            verdict = "INDEPENDENT_REHEARSAL_GATE_PARTIAL"
            readiness = "PARTIAL"
        else:
            verdict = "INDEPENDENT_REHEARSAL_GATE_FAIL"
            readiness = "NOT_READY"

        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "bundle_manifest": str(MANIFEST_PATH.relative_to(ROOT)),
            "tmp_rehearsal_dir": str(tmp),
            "hash_mismatches": hash_mismatches,
            "runs": {
                "qw2048": r2048,
                "qw2049": r2049,
            },
            "flags": flags,
            "pass_count": pass_count,
            "total_flags": total_flags,
            "verdict": verdict,
            "readiness": readiness,
            "required_next_step": "RUN_TRULY_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE",
        }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2051: INDEPENDENT REHEARSAL GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Readiness: **{out['readiness']}**",
        f"- pass_count: {out['pass_count']}/{out['total_flags']}",
        "",
        "## Manifest",
        f"- {out['bundle_manifest']}",
        "",
        "## Flags",
    ]
    for k, v in out["flags"].items():
        lines.append(f"- {k}: {v}")

    if out["hash_mismatches"]:
        lines.extend(["", "## Hash Mismatches"])
        for mm in out["hash_mismatches"]:
            lines.append(f"- {mm['path']}")

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

    print(f"[QW-2051] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2051] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2051] verdict={out['verdict']} pass_count={out['pass_count']}/{out['total_flags']}")


if __name__ == "__main__":
    main()
