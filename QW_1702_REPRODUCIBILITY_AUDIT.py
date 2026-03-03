#!/usr/bin/env python3
"""
QW-1702: Reproducibility audit (date-aware execution log).

Runs selected cornerstone scripts twice and checks:
1) return codes
2) runtime stability
3) byte-level reproducibility of output reports
4) selected metric extraction from produced JSON files
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1702_reproducibility_audit.json"
OUT_MD = ROOT / "RAPORT_QW1702_REPRODUCIBILITY_AUDIT.md"


@dataclass
class ScriptTask:
    script: str
    outputs: List[str]


TASKS = [
    ScriptTask(
        script="116_ALGEBRAIC_STRUCTURE_VERIFICATION.py",
        outputs=["report_116_algebraic_structure_verification.json"],
    ),
    ScriptTask(
        script="117_TOPOLOGICAL_CHARGES_AND_FAMILIES.py",
        outputs=["report_117_topological_charges_and_families.json"],
    ),
    ScriptTask(
        script="118_COMPOSITE_HIGGS_AND_EMERGENT_MASSES.py",
        outputs=["report_118_composite_higgs_and_emergent_masses.json"],
    ),
    ScriptTask(
        script="verify_values.py",
        outputs=[],
    ),
]


def sha256_file(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def normalize_json(obj):
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            kl = str(k).lower()
            if any(tok in kl for tok in ("time", "date", "timestamp", "generated", "created")):
                continue
            out[k] = normalize_json(v)
        return out
    if isinstance(obj, list):
        return [normalize_json(x) for x in obj]
    return obj


def semantic_json_hash(path: Path) -> Optional[str]:
    if not path.exists() or path.suffix.lower() != ".json":
        return None
    try:
        obj = json.loads(path.read_text(encoding="utf-8", errors="ignore"))
        norm = normalize_json(obj)
        raw = json.dumps(norm, sort_keys=True, ensure_ascii=False).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()
    except Exception:
        return None


def flatten_json(obj, prefix="") -> Dict[str, float]:
    out: Dict[str, float] = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            out.update(flatten_json(v, key))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            key = f"{prefix}[{i}]"
            out.update(flatten_json(v, key))
    elif isinstance(obj, (int, float)):
        out[prefix] = float(obj)
    return out


def run_once(script: Path) -> Dict[str, object]:
    t0 = time.perf_counter()
    proc = subprocess.run(
        ["python3", script.name],
        cwd=str(ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=1800,
    )
    dt = time.perf_counter() - t0
    return {
        "returncode": proc.returncode,
        "runtime_sec": dt,
        "stdout_tail": proc.stdout[-4000:],
        "stderr_tail": proc.stderr[-4000:],
    }


def extract_metric_sample(json_path: Path) -> Dict[str, float]:
    if not json_path.exists():
        return {}
    try:
        obj = json.loads(json_path.read_text(encoding="utf-8", errors="ignore"))
        flat = flatten_json(obj)
        selected = {
            k: v
            for k, v in flat.items()
            if any(tok in k.lower() for tok in ("error", "closure", "unit", "jacobi", "mass", "vev", "ckm", "winding"))
        }
        # Keep output compact and deterministic
        return dict(sorted(selected.items())[:40])
    except Exception:
        return {}


def main() -> None:
    start_utc = datetime.now(timezone.utc).isoformat()
    results = []

    for task in TASKS:
        script_path = ROOT / task.script
        if not script_path.exists():
            results.append(
                {
                    "script": task.script,
                    "exists": False,
                    "error": "script not found",
                }
            )
            continue

        # Run #1
        run1 = run_once(script_path)
        out_hashes_run1 = {}
        out_semantic_hashes_run1 = {}
        metrics_run1 = {}
        for out in task.outputs:
            p = ROOT / out
            out_hashes_run1[out] = sha256_file(p)
            out_semantic_hashes_run1[out] = semantic_json_hash(p)
            metrics_run1[out] = extract_metric_sample(p)

        # Run #2
        run2 = run_once(script_path)
        out_hashes_run2 = {}
        out_semantic_hashes_run2 = {}
        metrics_run2 = {}
        for out in task.outputs:
            p = ROOT / out
            out_hashes_run2[out] = sha256_file(p)
            out_semantic_hashes_run2[out] = semantic_json_hash(p)
            metrics_run2[out] = extract_metric_sample(p)

        deterministic_raw = (
            run1["returncode"] == 0
            and run2["returncode"] == 0
            and all(out_hashes_run1.get(k) == out_hashes_run2.get(k) for k in out_hashes_run1)
        )
        deterministic_semantic = (
            run1["returncode"] == 0
            and run2["returncode"] == 0
            and all(
                (
                    out_semantic_hashes_run1.get(k) == out_semantic_hashes_run2.get(k)
                    if out_semantic_hashes_run1.get(k) is not None and out_semantic_hashes_run2.get(k) is not None
                    else out_hashes_run1.get(k) == out_hashes_run2.get(k)
                )
                for k in out_hashes_run1
            )
        )

        results.append(
            {
                "script": task.script,
                "exists": True,
                "run1": run1,
                "run2": run2,
                "outputs": task.outputs,
                "output_hashes_run1": out_hashes_run1,
                "output_hashes_run2": out_hashes_run2,
                "output_semantic_hashes_run1": out_semantic_hashes_run1,
                "output_semantic_hashes_run2": out_semantic_hashes_run2,
                "metrics_run1": metrics_run1,
                "metrics_run2": metrics_run2,
                "deterministic_outputs_raw": deterministic_raw,
                "deterministic_outputs_semantic": deterministic_semantic,
            }
        )

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "started_utc": start_utc,
        "tasks": results,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1702: REPRODUCIBILITY AUDIT",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Liczba zadań: {len(results)}",
        "",
        "## Wyniki",
    ]

    for r in results:
        if not r.get("exists"):
            lines.append(f"- `{r['script']}`: ❌ brak pliku")
            continue
        rc1 = r["run1"]["returncode"]
        rc2 = r["run2"]["returncode"]
        dt1 = r["run1"]["runtime_sec"]
        dt2 = r["run2"]["runtime_sec"]
        det_raw = r["deterministic_outputs_raw"]
        det_sem = r["deterministic_outputs_semantic"]
        lines.append(
            f"- `{r['script']}`: rc1={rc1}, rc2={rc2}, t1={dt1:.3f}s, t2={dt2:.3f}s, deterministic_raw={det_raw}, deterministic_semantic={det_sem}"
        )

    lines.extend(
        [
            "",
            "## Artefakty",
            f"- JSON szczegółowy: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1702] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1702] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
