#!/usr/bin/env python3
"""QW-2265: RG canonical export bridge availability gate.

Strictly checks whether unresolved RG export-symbol refs detected by QW-2261
exist in explicit axiomatic bridge layer, without claiming non-axiomatic closure.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
BRIDGE_FILE = "FIN_L12_CANONICAL_EXPORT_AXIOMATIC_BRIDGE.lean"


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def lean_symbols(text: str) -> set[str]:
    names: set[str] = set()
    for pat in [r"^axiom\s+([A-Za-z0-9_']+)\s*:", r"^theorem\s+([A-Za-z0-9_']+)\s*:"]:
        names.update(re.findall(pat, text, flags=re.M))
    return names


def main() -> None:
    p2261 = load("proof_object_qw2261_rg_active_reference_locality_integrity.json")

    unresolved: set[str] = set()
    for row in p2261.get("active_instances_integrity", []):
        for r in row.get("unresolved_refs", []):
            unresolved.add(r)

    bridge_path = ROOT / BRIDGE_FILE
    bridge_exists = bridge_path.exists()
    bridge_text = bridge_path.read_text(encoding="utf-8") if bridge_exists else ""
    bridge_symbols = lean_symbols(bridge_text)

    bridged_unresolved = sorted([r for r in unresolved if r in bridge_symbols])
    unresolved_not_bridged = sorted([r for r in unresolved if r not in bridge_symbols])

    flags = {
        "q2261_unresolved_refs_present": len(unresolved) > 0,
        "bridge_file_present": bridge_exists,
        "canonical_export_symbol_present_in_bridge": "RG_CanonicalAction_to_WellPosedness_EXPORT" in bridge_symbols,
        "all_unresolved_refs_have_bridge_symbol": len(unresolved_not_bridged) == 0,
        "bridge_layer_is_axiomatic_scoped": "axiomatic layer" in bridge_text.lower() or "axiom" in bridge_text,
        "non_axiomatic_closure_not_claimed": True,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = flags["bridge_file_present"] and flags["canonical_export_symbol_present_in_bridge"]

    verdict = (
        "RG_CANONICAL_EXPORT_BRIDGE_AVAILABILITY_GATE_PASS_PARTIAL_AXIOMATIC_BRIDGE_AVAILABLE"
        if core_ok
        else "RG_CANONICAL_EXPORT_BRIDGE_AVAILABILITY_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "p2261": "proof_object_qw2261_rg_active_reference_locality_integrity.json",
            "bridge_file": BRIDGE_FILE,
        },
        "unresolved_refs_from_q2261": sorted(unresolved),
        "bridge_symbols": sorted(bridge_symbols),
        "bridged_unresolved_refs": bridged_unresolved,
        "unresolved_refs_not_bridged": unresolved_not_bridged,
        "scope_boundary": {
            "bridge_is_axiomatic_only": True,
            "non_axiomatic_closure_claimed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2265_rg_canonical_export_bridge_availability.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "n_unresolved_refs": len(unresolved),
        "n_bridged_unresolved_refs": len(bridged_unresolved),
        "n_unresolved_refs_not_bridged": len(unresolved_not_bridged),
        "bridged_unresolved_refs": bridged_unresolved,
        "unresolved_refs_not_bridged": unresolved_not_bridged,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "EITHER_FORMALIZE_LOCAL_IMPORT_CHAIN_OR_REPLACE_AXIOMATIC_BRIDGE_BY_NON_AXIOMATIC_DERIVATION",
    }

    out_json = ROOT / "report_qw2265_rg_canonical_export_bridge_availability_gate.json"
    out_md = ROOT / "RAPORT_QW2265_RG_CANONICAL_EXPORT_BRIDGE_AVAILABILITY_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2265: RG CANONICAL EXPORT BRIDGE AVAILABILITY GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- unresolved_refs_from_q2261: `{sorted(unresolved)}`",
                f"- bridged_unresolved_refs: `{bridged_unresolved}`",
                f"- unresolved_refs_not_bridged: `{unresolved_not_bridged}`",
                "- Scope: bridge is axiomatic-only, non-axiomatic closure still open.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "n_unresolved_refs": len(unresolved),
                "n_unresolved_refs_not_bridged": len(unresolved_not_bridged),
            }
        )
    )


if __name__ == "__main__":
    main()
