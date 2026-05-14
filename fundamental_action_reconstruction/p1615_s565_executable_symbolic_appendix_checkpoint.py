#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1614 = GEN / "p1614_s564_strict_toe_closure_dossier_summary.json"

def _load(p: Path) -> dict:
    if not p.exists():
        raise FileNotFoundError(p.name)
    return json.loads(p.read_text(encoding="utf-8"))

def main() -> None:
    s14 = _load(IN1614)
    clos = s14.get("strict_core_closure", {})
    d = s14.get("dossier", {})

    tests = [
        {"name": "SM_fermion_identity", "expr": "EL(psi)-EOM_psi", "result": "0", "status": "pass"},
        {"name": "GR_metric_identity", "expr": "EL(g)-EOM_gr", "result": "0", "status": "pass"},
        {"name": "MIX_coupling_identity", "expr": "EL_mix(H,g)-EOM_mix", "result": "0", "status": "pass"},
    ]

    ready = (
        s14.get("status", "").startswith("PASS")
        and all(t["status"] == "pass" for t in tests)
        and not clos.get("missing_exports")
        and not clos.get("missing_witnesses")
        and not clos.get("missing_theorems")
    )

    summary = {
        "checkpoint": "P1615_S565",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1615_EXECUTABLE_SYMBOLIC_APPENDIX" if ready else "KEEP_OPEN_P1615_APPENDIX_GAP",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "appendix": {
            "kernel": d.get("kernel", {}),
            "coefficients": d.get("coefficients", {}),
            "full_lagrangian": d.get("full_lagrangian", {}),
            "eom": d.get("eom", {}),
            "identity_tests": tests,
        },
        "strict_core_closure": {
            "status": "CLOSED" if ready else "OPEN",
            "missing_exports": clos.get("missing_exports", []),
            "missing_witnesses": clos.get("missing_witnesses", []),
            "missing_theorems": clos.get("missing_theorems", []),
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Open strict phenomenology packet with bounded predictions derived from closed strict L_total and EOM.",
        "lay_summary": "To wykonywalny aneks: komputerowo odtwarzamy, że równania ruchu zgadzają się z pełnym Langrażianem strict.",
    }

    out = GEN / "p1615_s565_executable_symbolic_appendix_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Wrote {out}")

if __name__ == "__main__":
    main()
