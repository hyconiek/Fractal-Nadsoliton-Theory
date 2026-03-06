#!/usr/bin/env python3
"""QW-2449: non-axiomatic dual canonical-export provider derivation attempt gate.

Runs strict lexical/definition scan for dual canonical-export provider symbols:
- RG_CanonicalAction_to_WellPosedness_EXPORT
- QFT_CanonicalAction_to_Positivity_EXPORT

No theorem-level/full-closure claim is allowed in this step.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent

RG_SYMBOL = "RG_CanonicalAction_to_WellPosedness_EXPORT"
QFT_SYMBOL = "QFT_CanonicalAction_to_Positivity_EXPORT"


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def run_rg(pattern: str) -> list[str]:
    proc = subprocess.run(
        ["rg", "-n", "--hidden", "--glob", "*.lean", pattern, str(ROOT)],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode not in (0, 1):
        raise RuntimeError(proc.stderr.strip() or "rg failed")
    return [ln for ln in proc.stdout.splitlines() if ln.strip()]


def parse_rg_hit(hit: str) -> tuple[Path, str]:
    # format: /abs/path/file.lean:line:content
    path_s, _, content = hit.partition(":")
    _, _, content = content.partition(":")
    return Path(path_s), content


def file_is_strict_non_axiomatic(path: Path) -> bool:
    if not path.exists():
        return False
    text = path.read_text(encoding="utf-8")
    has_axiom_decl = bool(re.search(r"^\s*axiom\s+", text, flags=re.MULTILINE))
    hits_derived_pending = "_DerivedOrPending" in text
    return (not has_axiom_decl) and (not hits_derived_pending)


def is_non_axiomatic_definition_hit(hit: str, symbol: str) -> bool:
    # strict gate: direct non-axiom declaration + source file must be axiom-token free
    path, line = parse_rg_hit(hit)
    def_pat = re.compile(rf"^\s*(theorem|lemma|def|abbrev|constant)\s+{re.escape(symbol)}\b")
    ax_pat = re.compile(rf"^\s*axiom\s+{re.escape(symbol)}\b")
    if ax_pat.search(line):
        return False
    if not def_pat.search(line):
        return False
    return file_is_strict_non_axiomatic(path)


def main() -> None:
    q2448 = load("report_qw2448_dual_single_foundation_v2_minimal_blocker_cut_gate.json")
    q2445 = load("report_qw2445_dual_nadsoliton_single_foundation_provider_execution_status_v2_gate.json")

    rg_refs = run_rg(rf"{RG_SYMBOL}|{QFT_SYMBOL}")
    rg_def_lines = run_rg(rf"^\s*(axiom|theorem|lemma|def|abbrev|constant)\s+{RG_SYMBOL}\b")
    qft_def_lines = run_rg(rf"^\s*(axiom|theorem|lemma|def|abbrev|constant)\s+{QFT_SYMBOL}\b")

    rg_non_axiomatic_defs = [ln for ln in rg_def_lines if is_non_axiomatic_definition_hit(ln, RG_SYMBOL)]
    qft_non_axiomatic_defs = [ln for ln in qft_def_lines if is_non_axiomatic_definition_hit(ln, QFT_SYMBOL)]

    flags = {
        "q2448_minimal_cut_ready": q2448.get("verdict")
        == "DUAL_SINGLE_FOUNDATION_V2_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED",
        "q2445_runtime_backed_execution_present": q2445.get("flags", {}).get("execution_attempt_performed") is True,
        "lexical_references_present": len(rg_refs) > 0,
        "rg_non_axiomatic_definition_exists": len(rg_non_axiomatic_defs) > 0,
        "qft_non_axiomatic_definition_exists": len(qft_non_axiomatic_defs) > 0,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if flags["q2448_minimal_cut_ready"] and flags["q2445_runtime_backed_execution_present"]:
        if flags["rg_non_axiomatic_definition_exists"] and flags["qft_non_axiomatic_definition_exists"]:
            verdict = "NON_AXIOMATIC_DUAL_CANONICAL_EXPORT_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PACKET_READY_FOR_MACHINE_DISCHARGE"
            required_next_step = "RUN_MACHINE_DISCHARGE_FOR_DUAL_NON_AXIOMATIC_EXPORT_PROVIDERS"
        else:
            verdict = "NON_AXIOMATIC_DUAL_CANONICAL_EXPORT_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_NO_NON_AXIOMATIC_PROVIDER_DEFINITION"
            required_next_step = "INTRODUCE_STRICT_NON_AXIOMATIC_PROVIDER_DEFINITIONS_FOR_BOTH_EXPORT_SYMBOLS"
    else:
        verdict = "NON_AXIOMATIC_DUAL_CANONICAL_EXPORT_PROVIDER_DERIVATION_ATTEMPT_GATE_FAIL"
        required_next_step = "REPAIR_QW2445_QW2448_PRECONDITIONS"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2448": "report_qw2448_dual_single_foundation_v2_minimal_blocker_cut_gate.json",
            "q2445": "report_qw2445_dual_nadsoliton_single_foundation_provider_execution_status_v2_gate.json",
        },
        "scan": {
            "reference_hits": rg_refs,
            "rg_definition_lines": rg_def_lines,
            "qft_definition_lines": qft_def_lines,
            "rg_non_axiomatic_definitions": rg_non_axiomatic_defs,
            "qft_non_axiomatic_definitions": qft_non_axiomatic_defs,
        },
        "scope_boundary": {
            "lexical_definition_scan_only": True,
            "provider_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    proof_path = ROOT / "proof_object_qw2449_non_axiomatic_dual_canonical_export_provider_derivation_attempt.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "counts": {
            "n_reference_hits": len(rg_refs),
            "n_rg_definition_lines": len(rg_def_lines),
            "n_qft_definition_lines": len(qft_def_lines),
            "n_rg_non_axiomatic_definitions": len(rg_non_axiomatic_defs),
            "n_qft_non_axiomatic_definitions": len(qft_non_axiomatic_defs),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2449_non_axiomatic_dual_canonical_export_provider_derivation_attempt_gate.json"
    out_md = ROOT / "RAPORT_QW2449_NON_AXIOMATIC_DUAL_CANONICAL_EXPORT_PROVIDER_DERIVATION_ATTEMPT_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2449: NON AXIOMATIC DUAL CANONICAL EXPORT PROVIDER DERIVATION ATTEMPT GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_reference_hits: `{out['counts']['n_reference_hits']}`",
                f"- n_rg_non_axiomatic_definitions: `{out['counts']['n_rg_non_axiomatic_definitions']}`",
                f"- n_qft_non_axiomatic_definitions: `{out['counts']['n_qft_non_axiomatic_definitions']}`",
                "",
                "## Wniosek rygorystyczny",
                "- Brak full-closure claim; granica theorem-level pozostaje jawna.",
                "- Warstwa provider wymaga jawnych, nieaksjomatycznych definicji eksportu (RG/QFT).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "counts": out["counts"]}))


if __name__ == "__main__":
    main()
