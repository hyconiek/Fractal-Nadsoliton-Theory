#!/usr/bin/env python3
"""P1837 S787 strict full-Lagrangian explicit anchor registry checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1747 = load("p1747_s697_strict_full_lagrangian_non_skeleton_and_bidirectional_witness_bundle_checkpoint.json")
    p1836 = load("p1836_s786_strict_full_lagrangian_non_skeleton_manifest_checkpoint.json")

    anchor = p1747.get("full_lagrangian_density_explicit_anchor", {})
    manifest = p1836.get("full_lagrangian_non_skeleton_manifest", {}).get("L_total", {})

    out = {
        "packet_id": "P1837",
        "stage_id": "S787",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "explicit_anchor_registry": {
            "reduced_non_skeleton_anchor_present": bool(anchor),
            "reduced_anchor_payload": anchor,
            "manifest_sector_count": len(manifest),
            "manifest_sector_keys": list(manifest.keys()),
        },
        "promotion_policy": "reduced_anchor_is_reference_only_not_full_nonproxy_sector_export",
        "technical_progress": "Existing explicit reduced non-skeleton Lagrangian anchor is registered and linked to the current full-sector manifest.",
        "proven": "Repo contains a concrete reduced explicit Lagrangian anchor, but it is insufficient for full nonproxy sector closure.",
        "open": "Sector-by-sector nonproxy term exports and EOM witnesses are still required for strict full-chain closure.",
        "false_pass_risk": "Confusing reduced explicit anchor with full nonproxy L_total export would be a category error.",
        "next_honest_step": "Use reduced anchor as consistency reference while delivering full sector exports through P1834 gates.",
        "lay_explanation": "Mamy już konkretny przykład pełniejszego wzoru w wersji zredukowanej, ale to jeszcze nie jest komplet całej teorii.",
    }

    path = GEN / "p1837_s787_strict_full_lagrangian_explicit_anchor_registry_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
