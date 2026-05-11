#!/usr/bin/env python3
"""P1230: check minimal section completeness and symmetry-state discipline for W1."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"
REQUIRED_SECTIONS = ["S1_statement", "S2_assumptions", "S3_scope_boundary"]


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=Path, default=GEN / "p1229_w1_minimal_proof_packet_index.json")
    parser.add_argument("--p1228", type=Path, default=GEN / "p1228_w1_real_candidate_evidence_note_summary.json")
    parser.add_argument("--p1223-input", type=Path, default=GEN / "p1223_w1_input_payload.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1230_w1_section_completeness_and_symmetry_summary.json")
    args = parser.parse_args()

    idx = _load(args.index)
    p1228 = _load(args.p1228)
    p1223_input = _load(args.p1223_input)

    ref = p1228.get("candidate_ref", "")
    entry = idx.get("entries", {}).get(ref, {}) if isinstance(idx.get("entries", {}), dict) else {}
    sections = entry.get("sections", []) if isinstance(entry.get("sections", []), list) else []
    missing_sections = [s for s in REQUIRED_SECTIONS if s not in sections]

    available = p1223_input.get("available_inputs", {})
    symmetry_breaking_axiom_declared = bool(available.get("explicit_symmetry_breaking_axiom_declared", False))
    scoped_non_strict_tag = bool(available.get("scoped_non_strict_tag_if_axiom_used", False))

    if symmetry_breaking_axiom_declared and scoped_non_strict_tag:
        symmetry_state = "AXIOM_DECLARED_NON_STRICT"
    else:
        symmetry_state = "UNBROKEN_IN_STRICT_CORE"

    out = {
        "packet": "P1230",
        "as_of": "2026-05-11",
        "candidate_ref": ref,
        "required_sections": REQUIRED_SECTIONS,
        "present_sections": sections,
        "missing_sections": missing_sections,
        "section_completeness_ok": len(missing_sections) == 0,
        "symmetry_state": symmetry_state,
        "symmetry_breaking_axiom_declared": symmetry_breaking_axiom_declared,
        "scoped_non_strict_tag_if_axiom_used": scoped_non_strict_tag,
        "kernel_track_state": "strict_gate_operational_track_active; no legacy-role transfer claim",
        "nadsoliton_ontology_state": "nadsoliton_is_primordial_information_no_substrate_below",
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Completeness + symmetry discipline checkpoint; no strict-core closure claim.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1230] completeness_ok={out['section_completeness_ok']} symmetry_state={symmetry_state} wrote {args.out}")


if __name__ == "__main__":
    main()
