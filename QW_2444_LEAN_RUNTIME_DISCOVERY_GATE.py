#!/usr/bin/env python3
"""QW-2444: Lean runtime discovery gate.

Strict environment diagnostics for Lean machine-check runtime:
- detects candidate executables,
- verifies executability and version probe,
- emits explicit blocker state if runtime unavailable.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
LOCAL_ELAN_HOME = ROOT / ".elan"
LOCAL_HOME = ROOT / ".home_lean"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def which_all(cmd: str) -> list[Path]:
    out: list[Path] = []
    seen: set[str] = set()
    for p in os.environ.get("PATH", "").split(os.pathsep):
        if not p:
            continue
        c = Path(p) / cmd
        k = str(c.resolve()) if c.exists() else str(c)
        if c.exists() and os.access(c, os.X_OK) and k not in seen:
            seen.add(k)
            out.append(c)
    return out


def probe(path: Path) -> dict[str, Any]:
    env = os.environ.copy()
    if str(path).startswith(str(LOCAL_ELAN_HOME)):
        env["ELAN_HOME"] = str(LOCAL_ELAN_HOME)
        env["HOME"] = str(LOCAL_HOME)
    try:
        proc = subprocess.run([str(path), "--version"], capture_output=True, text=True, check=False, env=env)
        return {
            "path": str(path),
            "exists": path.exists(),
            "executable": os.access(path, os.X_OK),
            "exit_code": proc.returncode,
            "stdout": proc.stdout.strip(),
            "stderr": proc.stderr.strip(),
            "ok": proc.returncode == 0,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "path": str(path),
            "exists": path.exists(),
            "executable": os.access(path, os.X_OK) if path.exists() else False,
            "exit_code": 127,
            "stdout": "",
            "stderr": str(exc),
            "ok": False,
        }


def main() -> None:
    candidates: list[Path] = []
    candidates.extend(which_all("lean"))

    for fixed in [
        LOCAL_ELAN_HOME / "bin/lean",
        Path.home() / ".elan/bin/lean",
        Path("/usr/bin/lean"),
        Path("/usr/local/bin/lean"),
        Path("/opt/lean/bin/lean"),
    ]:
        if fixed not in candidates:
            candidates.append(fixed)

    probes = [probe(p) for p in candidates]
    valid = [p for p in probes if p.get("ok")]
    selected = valid[0]["path"] if valid else None

    flags = {
        "candidate_scan_completed": True,
        "at_least_one_candidate_path_checked": len(candidates) > 0,
        "lean_runtime_available": selected is not None,
        "selected_runtime_probe_ok": selected is not None,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if selected is not None:
        verdict = "LEAN_RUNTIME_DISCOVERY_GATE_PASS_RUNTIME_AVAILABLE"
        required_next_step = "RERUN_DUAL_SINGLE_FOUNDATION_PROVIDER_EXECUTION_WITH_SELECTED_RUNTIME"
    else:
        verdict = "LEAN_RUNTIME_DISCOVERY_GATE_PASS_PARTIAL_BLOCKED_BY_NO_RUNTIME"
        required_next_step = "ATTACH_OR_INSTALL_LEAN_RUNTIME_AND_RERUN"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "candidates": [str(c) for c in candidates],
        "probes": probes,
        "selected_runtime": selected,
        "scope_boundary": {
            "environment_diagnostics_only": True,
            "theorem_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2444_lean_runtime_discovery.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "selected_runtime": selected,
        "n_candidates": len(candidates),
        "n_valid": len(valid),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2444_lean_runtime_discovery_gate.json"
    out_md = ROOT / "RAPORT_QW2444_LEAN_RUNTIME_DISCOVERY_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2444: LEAN RUNTIME DISCOVERY GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_candidates: `{len(candidates)}`",
                f"- n_valid: `{len(valid)}`",
                f"- selected_runtime: `{selected}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "selected_runtime": selected, "n_candidates": len(candidates)}))


if __name__ == "__main__":
    main()
