"""P3151/S2101: representation-content source-selection audit for R_SM^ax.

P3150 derived the hypercharge ray only after the one-family SM-like field
pattern had already been installed.  This packet attacks the next missing
object: can the representation *content* itself be selected by current strict
objects?  We construct a finite shape scanner over small SU(3)xSU(2) slots and
show that the local algebraic/Yukawa/anomaly filters admit multiple shapes, so
R_SM^ax is not uniquely source-selected without importing SM family data.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3151_s2101_rsm_representation_content_source_selection_audit.json"
MD = GEN / "p3151_s2101_rsm_representation_content_source_selection_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3150": GEN / "p3150_s2100_rsm_hypercharge_ray_source_selection_audit.json",
    "P3149": GEN / "p3149_s2099_brs_ltotal_interface_invariance_audit.json",
    "P3148": GEN / "p3148_s2098_sm_representation_registry_completion_audit.json",
    "P1962": GEN / "p1962_s912_strict_matter_higgs_brst_extension_registry_audit.json",
}

FIELDS = ["q", "u", "d", "l", "e", "h"]
# Small strict-local receiver alphabet: singlet/fundamental/antifundamental for SU3
# and singlet/doublet for SU2.  t3 is triality mod 3; dim3 is complex dimension.
COLOR = {
    "1": {"dim3": 1, "t3": 0},
    "3": {"dim3": 3, "t3": 1},
    "bar3": {"dim3": 3, "t3": -1},
}
WEAK = {"1": {"dim2": 1}, "2": {"dim2": 2}}
SM_SHAPE = {
    "q": ("3", "2"),
    "u": ("bar3", "1"),
    "d": ("bar3", "1"),
    "l": ("1", "2"),
    "e": ("1", "1"),
    "h": ("1", "2"),
}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def contains_su2_singlet(reps: tuple[str, ...]) -> bool:
    # Only the 1/2 spin alphabet is scanned.  A tensor of SU2 doublets contains
    # a singlet iff the number of doublet factors is even.
    return sum(1 for r in reps if r == "2") % 2 == 0


def yukawa_shape_ok(shape: dict[str, tuple[str, str]]) -> bool:
    triples = [("q", "h", "u"), ("q", "h", "d"), ("l", "h", "e")]
    for a, b, c in triples:
        t_sum = sum(COLOR[shape[x][0]]["t3"] for x in (a, b, c)) % 3
        if t_sum != 0:
            return False
        if not contains_su2_singlet(tuple(shape[x][1] for x in (a, b, c))):
            return False
    return True


def coarse_anomaly_shape_ok(shape: dict[str, tuple[str, str]]) -> bool:
    # Shape-only constraints that do not import hypercharge values: no net SU3
    # chirality among non-Higgs Weyl slots and even number of SU2 doublets after
    # color multiplicity (Witten parity).  These are necessary consistency gates,
    # not a full physical source theorem.
    su3_chirality = sum(COLOR[shape[f][0]]["t3"] * WEAK[shape[f][1]]["dim2"] for f in FIELDS if f != "h")
    witten_doublets = sum(COLOR[shape[f][0]]["dim3"] for f in FIELDS if f != "h" and shape[f][1] == "2")
    return su3_chirality % 3 == 0 and witten_doublets % 2 == 0


def shape_to_record(shape: dict[str, tuple[str, str]]) -> dict[str, str]:
    return {f: f"({shape[f][0]},{shape[f][1]})" for f in FIELDS}


def scan_shapes() -> dict[str, Any]:
    options = list(itertools.product(COLOR, WEAK))
    total = 0
    yukawa_pass: list[dict[str, tuple[str, str]]] = []
    full_pass: list[dict[str, tuple[str, str]]] = []
    for combo in itertools.product(options, repeat=len(FIELDS)):
        total += 1
        shape = dict(zip(FIELDS, combo, strict=True))
        if yukawa_shape_ok(shape):
            yukawa_pass.append(shape)
            if coarse_anomaly_shape_ok(shape):
                full_pass.append(shape)
    sm_matches = [s for s in full_pass if s == SM_SHAPE]
    distinct_dimension_patterns = sorted({tuple((COLOR[s[f][0]]["dim3"], WEAK[s[f][1]]["dim2"]) for f in FIELDS) for s in full_pass})
    return {
        "alphabet_slots_per_field": len(options),
        "total_shapes_scanned": total,
        "yukawa_shape_pass_count": len(yukawa_pass),
        "coarse_anomaly_pass_count": len(full_pass),
        "sm_shape_present": len(sm_matches) == 1,
        "unique_shape_selected": len(full_pass) == 1,
        "distinct_dimension_pattern_count": len(distinct_dimension_patterns),
        "sample_passing_shapes": [shape_to_record(s) for s in full_pass[:12]],
        "sm_shape": shape_to_record(SM_SHAPE),
    }


def build_payload() -> dict[str, Any]:
    scan = scan_shapes()
    candidates = [
        {"candidate": "finite small-representation Yukawa/anomaly shape scanner", "selects_sm_shape": scan["unique_shape_selected"] and scan["sm_shape_present"], "strict_source_law": False, "noncircular": True, "verdict": "constructive obstruction: SM shape is allowed but not unique"},
        {"candidate": "P3150 hypercharge ray", "selects_sm_shape": False, "strict_source_law": False, "noncircular": True, "verdict": "depends on representation slots; cannot source the slots it assumes"},
        {"candidate": "P3148 R_SM^ax registry", "selects_sm_shape": True, "strict_source_law": False, "noncircular": False, "verdict": "installed target registry, circular as a source"},
        {"candidate": "P1960/P1961 local Lie/BRST interface", "selects_sm_shape": False, "strict_source_law": False, "noncircular": True, "verdict": "tests compatibility after fields are supplied, not field content selection"},
    ]
    accepted = [r for r in candidates if r["selects_sm_shape"] and r["strict_source_law"] and r["noncircular"]]
    counts = {**{k: scan[k] for k in ["total_shapes_scanned", "yukawa_shape_pass_count", "coarse_anomaly_pass_count", "distinct_dimension_pattern_count"]}, "candidate_source_rows": len(candidates), "strict_accepted_source_rows": len(accepted)}
    return {
        "status": "P3151_RSM_REPRESENTATION_CONTENT_SOURCE_SELECTION_MULTI_WITNESS_NO_STRICT_SOURCE",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "constructed_object": {"name": "R_shape^scan finite representation-content source-selection scanner", "classification": "bounded_multi_witness_obstruction", "scope": "small SU3 singlet/fundamental/antifundamental and SU2 singlet/doublet receiver alphabet for six P3148 slots"},
        "finite_scan": scan,
        "candidate_source_rows": candidates,
        "finite_theorem": {"name": "P3151_T1_representation_shape_nonselection_obstruction", "statement": "On the scanned small receiver alphabet, the local Yukawa singlet conditions and shape-only anomaly/parity gates do not uniquely select the P3148 one-family representation content.  The SM shape is present, but the pass set contains multiple shapes and multiple dimension patterns.  Therefore current algebraic compatibility checks are receivers, not a strict nadsoliton source law for the five chiral multiplets plus Higgs.", "finite_counts": counts},
        "decision": {"bounded_result": "P3151 constructs the requested representation-content object and proves a finite non-uniqueness obstruction: compatibility is not selection.", "why_not_strict": "The only row that selects the exact SM shape is the installed R_SM^ax registry itself, which is circular.  The noncircular local scanner admits many alternatives and no current strict object supplies an external minimization, chirality, generation, or source law that picks the SM row.", "next_honest_step": "Pivot to one missing physical bridge not already closed by P3150/P3151: audit charge-unit normalization for the hypercharge ray, or audit GR/EH nonproxy coupling for the axiom branch.  The most proof-grade next step is P3152: construct a finite charge-unit normalization obstruction/witness for Y_SM^ray, testing whether any current strict object fixes the scale Y(H)=1/2 without importing the SM electric-charge convention.", "negative_export_flags": {"strict_representation_content_source_exported": False, "strict_SM_generation_exported": False, "unit_bearing_L_total_exported": False, "global_BV_BRST_exported": False, "strict_GR_generation_exported": False, "selector_closure_exported": False, "bridge_or_role_transfer_exported": False, "ToE_closure_exported": False}},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    scan = payload["finite_scan"]
    lines = ["# P3151/S2101 R_SM representation-content source-selection audit", "", f"Status: `{payload['status']}`", "", "## Constructed object", f"- `{payload['constructed_object']['name']}`", f"- Classification: `{payload['constructed_object']['classification']}`", f"- Scope: `{payload['constructed_object']['scope']}`", "", "## Finite theorem", f"`{th['name']}`: {th['statement']}", "", "## Finite counts"]
    for k, v in th["finite_counts"].items():
        lines.append(f"- `{k}`: `{v}`")
    lines.extend(["", "## SM row", json.dumps(scan["sm_shape"], sort_keys=True), "", "## Sample passing alternatives"])
    for row in scan["sample_passing_shapes"]:
        lines.append(f"- `{json.dumps(row, sort_keys=True)}`")
    lines.extend(["", "## Candidate source rows"])
    for row in payload["candidate_source_rows"]:
        lines.append(f"- `{row['candidate']}`: selects SM shape `{row['selects_sm_shape']}`, strict source `{row['strict_source_law']}`, noncircular `{row['noncircular']}`; {row['verdict']}")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Why this is not strict", payload["decision"]["why_not_strict"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3151/S2101 R_SM representation-content source-selection audit", "## P3151/S2101 R_SM representation-content source-selection audit\n\n`P3151/S2101` constructs `R_shape^scan`, a finite representation-content source-selection scanner for the P3148/P3150 SM axiom branch.  Over a small `SU(3)` singlet/fundamental/antifundamental and `SU(2)` singlet/doublet alphabet, the SM representation shape is present, but local Yukawa singlet gates plus shape-only anomaly/parity gates leave multiple passing shapes and dimension patterns.  Thus current compatibility checks do not source-select the five chiral multiplets plus Higgs, and no strict SM generation, unit-bearing `L_total`, GR/EH closure, bridge/role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3151/S2101 representation content remains installed", "## P3151/S2101 representation content remains installed\n\n`P3151/S2101` shows that `R_SM^ax` remains an installed representation registry: the finite local algebraic receiver tests admit the SM shape but do not uniquely select it.  A physical Lagrangian/EOM package would still need a noncircular source of field content, charge normalization, unit measure, and global BV/BRST plus GR/EH coupling.\n")
    append_once(AGENTS, "Current R_SM representation-content source-selection guardrail (P3151/S2101, 2026-07-13)", "## Current R_SM representation-content source-selection guardrail (P3151/S2101, 2026-07-13)\n\n- P3151 constructs `R_shape^scan`, a finite source-selection scanner for the P3148/P3150 representation content.\n- The scanned small `SU(3)xSU(2)` receiver alphabet contains the SM one-family shape, but Yukawa singlet gates plus shape-only anomaly/parity gates admit multiple shapes and multiple dimension patterns.\n- This proves compatibility is not source-selection: current strict artifacts do not noncircularly select the five chiral multiplets plus Higgs.\n- Do not promote `R_shape^scan`, `R_SM^ax`, or `Y_SM^ray` to strict SM generation, unit-bearing `L_total`, global BV/BRST, GR/EH closure, selector closure, bridge/role transfer, or ToE.\n- Next honest move: pivot to charge-unit normalization for `Y_SM^ray` or GR/EH nonproxy coupling; the most finite proof-grade option is a P3152 charge-unit normalization obstruction/witness for fixing `Y(H)=1/2` without importing the SM electric-charge convention.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
