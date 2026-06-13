#!/usr/bin/env python3
from __future__ import annotations

import hashlib, itertools, json, subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2670_s1620_higher_order_boolean_cross_invariant_lift_audit.json"
MD = GEN / "p2670_s1620_higher_order_boolean_cross_invariant_lift_audit.md"
P2669 = GEN / "p2669_s1619_boundary_cycle_cross_invariant_boolean_ansatz_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
NEGATIVE_EXPORT_FLAGS = [
    "higher_order_boolean_cross_invariant_source_exported", "auxiliary_lift_origin_exported",
    "pair_variable_physical_origin_exported", "boundary_sector_variable_physical_origin_exported",
    "boundary_phase_bit_target_exported_unconditionally", "intrinsic_entropy_level_exported",
    "target_independent_beta_source_exported", "q_w_2191_discharged", "role_bearing_ltotal_reenabled", "toe_closure_claimed",
]

def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None

def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()

def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}

def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(["rg", "-n", pattern, ".", "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json", "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:120]}

def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "higher_order_boolean_content": "higher-order invariant|higher order invariant|Boolean lift|auxiliary lift|ANF|GF\\(2\\)",
        "physical_origin_content": "physical origin|coding convention|not from labels|pair variable|boundary-sector variable|source theorem",
        "boundary_cycle_sector_content": "boundary-cycle|boundary cycle|sector swap|square holonomy|non-exact sector|boundary sector",
        "selector_content": "selector|pair2_positive|pair1/pair2|w_break|QW-2191",
        "nonclosure_guard_content": "L_total|ToE closure|beta source|nonexport|no false pass",
    }
    return {"tool": "rg", "mode": "content-first higher-order Boolean cross-invariant physical-origin antiduplication audit", "patterns": {k: rg_count(v) for k, v in patterns.items()}}

def upstream_consistency() -> dict[str, Any]:
    p2669 = load_json(P2669)
    d = p2669.get("closure_decision", {})
    return {
        "p2669_mathematical_candidates_exist": d.get("mathematical_candidates_exist") is True,
        "p2669_no_mixed_candidate": d.get("mixed_candidates_exist") is False,
        "p2669_no_passing_source": d.get("passing_cross_invariant_count") == 0,
        "p2669_no_convention_free_source": d.get("convention_free_source_exported_for_any_candidate") is False,
    }

def truth_from_mask(mask: int, p: int, s: int, u: int) -> int:
    return (mask >> ((p << 2) | (s << 1) | u)) & 1

def mobius_anf_terms(mask: int) -> list[str]:
    vals = [truth_from_mask(mask, (i>>2)&1, (i>>1)&1, i&1) for i in range(8)]
    coeff = vals[:]
    for bit in range(3):
        for i in range(8):
            if i & (1 << bit): coeff[i] ^= coeff[i ^ (1 << bit)]
    names = ["1", "u", "sector1", "sector1*u", "pair2", "pair2*u", "pair2*sector1", "pair2*sector1*u"]
    return [names[i] for i, c in enumerate(coeff) if c]

def higher_order_witness() -> dict[str, Any]:
    rows=[]
    for mask in range(256):
        odd = all(truth_from_mask(mask,p,s,u) ^ truth_from_mask(mask,p,1-s,u) == 1 for p,s,u in itertools.product([0,1], repeat=3))
        ties = all(truth_from_mask(mask,1,1,u)==1 and truth_from_mask(mask,1,0,u)==0 for u in [0,1])
        branch_sensitive = any(truth_from_mask(mask,0,s,u) != truth_from_mask(mask,1,s,u) for s,u in itertools.product([0,1], repeat=2))
        depends_u = any(truth_from_mask(mask,p,s,0) != truth_from_mask(mask,p,s,1) for p,s in itertools.product([0,1], repeat=2))
        terms = mobius_anf_terms(mask)
        contains_higher = any("u" in t or "pair2*sector1" in t for t in terms)
        source = False
        passes = odd and ties and branch_sensitive and contains_higher and source
        rows.append({"mask": mask, "anf_terms": terms, "sector_swap_odd_for_all_auxiliary_values": odd, "ties_pair2_positive_to_sector1_for_all_auxiliary_values": ties, "branch_sensitive": branch_sensitive, "depends_on_auxiliary_lift": depends_u, "contains_higher_order_or_auxiliary_term": contains_higher, "convention_free_physical_origin_exported": source, "passes_now": passes})
    candidates=[r for r in rows if r["sector_swap_odd_for_all_auxiliary_values"] and r["ties_pair2_positive_to_sector1_for_all_auxiliary_values"] and r["branch_sensitive"]]
    higher=[r for r in candidates if r["contains_higher_order_or_auxiliary_term"]]
    return {"statement": "P2670 exhaustively enumerates all 256 Boolean functions f(pair2, sector1, auxiliary_lift) to test whether adding one higher-order/auxiliary bit can improve on P2669. Higher-order mathematical forms exist, but every passing mathematical form still requires a convention-free physical origin for pair2, sector1, and any auxiliary lift; none is exported by current bridge-completed nadsoliton dynamics.", "total_functions": len(rows), "candidate_count": len(candidates), "higher_order_candidate_count": len(higher), "passing_source_count": sum(r["passes_now"] for r in rows), "auxiliary_dependent_candidate_count": sum(1 for r in candidates if r["depends_on_auxiliary_lift"]), "sample_candidates": candidates[:12], "convention_free_physical_origin_exported_for_any_candidate": any(r["convention_free_physical_origin_exported"] for r in rows)}

def closure_decision(w: dict[str, Any]) -> dict[str, Any]:
    return {"decision": "P2670_HIGHER_ORDER_BOOLEAN_CROSS_INVARIANT_LIFT_AUDIT__NO_PHYSICAL_ORIGIN_SOURCE", "professorial_verdict": "P2670 is a sharper finite no-false-pass extension of P2669. A one-bit auxiliary/higher-order Boolean lift produces more mathematical candidate functions, including auxiliary-dependent candidates, but it does not derive the coding of pair2, boundary sector 1, or the auxiliary lift from bridge-completed nadsoliton boundary dynamics. Thus higher-order Boolean freedom does not source the entropy-bit anchor.", "next_honest_step": "Stop adding formal Boolean lifts unless a concrete bridge-completed boundary-dynamics observable is supplied. The next honest step is a physical-origin audit of candidate observables that could define pair2, sector1, and any auxiliary bit without imported labels; if no observable passes, promote P2669/P2670 to a no-go for Boolean cross-invariant entropy-bit sourcing.", "mathematical_higher_order_candidates_exist": w["higher_order_candidate_count"] > 0, "passing_source_count": w["passing_source_count"], "boundary_phase_bit_target_exported_now": False, "unconditional_uv_unit_selected_now": False, "beta_source_exported_now": False, "qw2191_discharged_now": False, "role_bearing_ltotal_now": False, "toe_closure_now": False}

def write_markdown(payload: dict[str, Any]) -> None:
    w=payload["higher_order_witness"]; d=payload["closure_decision"]
    lines=["# P2670/S1620 higher-order Boolean cross-invariant lift audit", "", f"Status: `{payload['status']}`", "", "## Content-first repo grep"]
    for k,v in payload["semantic_rg_antiduplication_audit"]["patterns"].items(): lines.append(f"- `{k}`: {v['count']} hits")
    lines += ["", "## Higher-order witness", w["statement"], f"Total functions: `{w['total_functions']}`.", f"Candidate count: `{w['candidate_count']}`.", f"Higher-order candidate count: `{w['higher_order_candidate_count']}`.", f"Auxiliary-dependent candidate count: `{w['auxiliary_dependent_candidate_count']}`.", f"Passing source count: `{w['passing_source_count']}`.", "", "## Verdict", d["professorial_verdict"], "", "## Next honest step", d["next_honest_step"], "", "## Negative exports"]
    for k,v in payload["negative_export_flags"].items(): lines.append(f"- `{k}`: `{v}`")
    MD.write_text("\n".join(lines)+"\n", encoding="utf-8")

def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    witness=higher_order_witness(); decision=closure_decision(witness)
    payload={"status":"P2670_HIGHER_ORDER_BOOLEAN_CROSS_INVARIANT_LIFT_AUDIT_NO_FALSE_PASS", "semantic_rg_antiduplication_audit": semantic_rg_audit(), "source_hashes":{"P2669": sha256_file(P2669), "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET), "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT)}, "upstream_consistency": upstream_consistency(), "higher_order_witness": witness, "closure_decision": decision, "negative_export_flags": {k: False for k in NEGATIVE_EXPORT_FLAGS}}
    payload["payload_sha256"] = sha256_json({k:v for k,v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False)+"\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2670/S1620 higher-order Boolean lift guard", "## P2670/S1620 higher-order Boolean lift guard\n\n`P2670/S1620` exhaustively audits all `256` Boolean functions `f(pair2, sector1, auxiliary_lift)`.  Higher-order/auxiliary mathematical candidates exist, but none supplies a convention-free physical origin for the pair variable, boundary-sector variable, or auxiliary lift from bridge-completed nadsoliton boundary dynamics.  Therefore no boundary-phase bit target, intrinsic UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2670/S1620 higher-order Boolean lift Ltotal guard", "## P2670/S1620 higher-order Boolean lift Ltotal guard\n\n`P2670/S1620` keeps `L_total` closed to higher-order Boolean cross-invariant entropy terms.  Adding an auxiliary Boolean lift is not a variational source unless the auxiliary bit and the pair/sector coding are derived from strict nadsoliton boundary dynamics rather than imported labels.\n")
    return payload

if __name__ == "__main__":
    main()
