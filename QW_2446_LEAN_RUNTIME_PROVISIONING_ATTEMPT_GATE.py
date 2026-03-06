#!/usr/bin/env python3
"""QW-2446: Lean runtime provisioning attempt gate.

Attempts local (workspace-scoped) Lean runtime provisioning via elan.
Records network/runtime blockers explicitly.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
ELAN_DIR = ROOT / ".elan"
SCRIPT_PATH = Path("/tmp/elan-init.sh")
ELAN_SCRIPT_URL = "https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def run(cmd: list[str]) -> dict[str, Any]:
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=False)
    return {
        "cmd": cmd,
        "exit_code": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
    }


def main() -> None:
    q2444 = json.loads((ROOT / "report_qw2444_lean_runtime_discovery_gate.json").read_text(encoding="utf-8"))

    runtime_already_available = q2444.get("verdict") == "LEAN_RUNTIME_DISCOVERY_GATE_PASS_RUNTIME_AVAILABLE"
    q2444_runtime_unavailable_confirmed = (
        q2444.get("verdict") == "LEAN_RUNTIME_DISCOVERY_GATE_PASS_PARTIAL_BLOCKED_BY_NO_RUNTIME"
    )

    curl_run: dict[str, Any] | None = None
    install_run: dict[str, Any] | None = None
    elan_ver_run: dict[str, Any] | None = None
    lean_ver_run: dict[str, Any] | None = None

    if not runtime_already_available:
        curl_run = run(["curl", "-L", ELAN_SCRIPT_URL, "-o", str(SCRIPT_PATH), "--max-time", "30"])
        if curl_run["exit_code"] == 0:
            install_run = run(
                [
                    "bash",
                    str(SCRIPT_PATH),
                    "-y",
                    "--default-toolchain",
                    "stable",
                    "--no-modify-path",
                    "--elan-dir",
                    str(ELAN_DIR),
                ]
            )
            if install_run["exit_code"] == 0:
                elan_ver_run = run([str(ELAN_DIR / "bin/elan"), "--version"])
                lean_ver_run = run([str(ELAN_DIR / "bin/lean"), "--version"])

    curl_err = ""
    if curl_run is not None:
        curl_err = (curl_run.get("stdout", "") + "\n" + curl_run.get("stderr", "")).lower()
    dns_failure = "could not resolve host" in curl_err

    runtime_installed_locally = bool(
        install_run is not None
        and install_run.get("exit_code") == 0
        and elan_ver_run is not None
        and elan_ver_run.get("exit_code") == 0
        and lean_ver_run is not None
        and lean_ver_run.get("exit_code") == 0
    )

    flags = {
        "runtime_available_before_attempt": runtime_already_available,
        "q2444_runtime_available": runtime_already_available,
        "q2444_runtime_unavailable_confirmed": q2444_runtime_unavailable_confirmed,
        "provisioning_attempt_executed": not runtime_already_available,
        "provisioning_skipped_due_to_existing_runtime": runtime_already_available,
        "curl_download_succeeded": bool(curl_run is not None and curl_run["exit_code"] == 0),
        "dns_resolution_failed": dns_failure,
        "runtime_installed_locally": runtime_installed_locally,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if runtime_already_available:
        verdict = "LEAN_RUNTIME_PROVISIONING_ATTEMPT_GATE_PASS_SKIPPED_RUNTIME_ALREADY_AVAILABLE"
        required_next_step = "USE_AVAILABLE_RUNTIME_FOR_PROVIDER_EXECUTION_GATES"
    elif runtime_installed_locally:
        verdict = "LEAN_RUNTIME_PROVISIONING_ATTEMPT_GATE_PASS_RUNTIME_INSTALLED_LOCALLY"
        required_next_step = "RERUN_QW2444_AND_QW2445_WITH_LOCAL_RUNTIME"
    elif dns_failure:
        verdict = "LEAN_RUNTIME_PROVISIONING_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_NETWORK_DNS"
        required_next_step = "PROVIDE_OFFLINE_LEAN_RUNTIME_OR_ENABLE_DNS_FOR_RAW_GITHUB"
    else:
        verdict = "LEAN_RUNTIME_PROVISIONING_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_ENVIRONMENT"
        required_next_step = "COLLECT_PROVISIONING_ERROR_AND_ATTACH_RUNTIME_MANUALLY"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2444_lean_runtime_discovery_gate.json",
        "attempt": {
            "elan_script_url": ELAN_SCRIPT_URL,
            "elan_dir": str(ELAN_DIR),
            "runtime_already_available_before_attempt": runtime_already_available,
            "curl": curl_run,
            "install": install_run,
            "elan_version": elan_ver_run,
            "lean_version": lean_ver_run,
        },
        "scope_boundary": {
            "environment_provisioning_only": True,
            "theorem_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2446_lean_runtime_provisioning_attempt.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": proof_obj["source"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "runtime_available_before_attempt": runtime_already_available,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2446_lean_runtime_provisioning_attempt_gate.json"
    out_md = ROOT / "RAPORT_QW2446_LEAN_RUNTIME_PROVISIONING_ATTEMPT_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    md_lines = [
        "# RAPORT QW-2446: LEAN RUNTIME PROVISIONING ATTEMPT GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        f"- runtime_available_before_attempt: `{flags['runtime_available_before_attempt']}`",
        f"- provisioning_attempt_executed: `{flags['provisioning_attempt_executed']}`",
        f"- provisioning_skipped_due_to_existing_runtime: `{flags['provisioning_skipped_due_to_existing_runtime']}`",
        f"- curl_download_succeeded: `{flags['curl_download_succeeded']}`",
        f"- dns_resolution_failed: `{flags['dns_resolution_failed']}`",
        f"- runtime_installed_locally: `{flags['runtime_installed_locally']}`",
    ]
    if runtime_already_available:
        md_lines.extend(
            [
                "",
                "## Interpretacja",
                "- Runtime Lean był już dostępny przed uruchomieniem tej bramki.",
                "- `runtime_installed_locally=false` oznacza brak nowej instalacji w tej bramce, nie brak runtime.",
            ]
        )
    out_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "dns_resolution_failed": flags["dns_resolution_failed"]}))


if __name__ == "__main__":
    main()
