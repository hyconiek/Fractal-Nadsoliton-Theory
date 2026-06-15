#!/usr/bin/env python3
"""P2770/S1720: kernel characteristic expressivity audit.

The user asked whether the current kernel formulas fully express the intended
FIN ontology: the nadsoliton as the single primordial foundation, a solitonic
self-learning neural-network-like system, geometrically self-coupled to itself.
This bounded audit treats that request as an expressivity/provenance question,
not as permission to promote closure.  It compares the explicit formula tokens
in K_legacy_ont and K_strict_gate against the required typed characteristics and
then scans current artifacts for theorem-level exports that could bridge the
missing semantics.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2770_s1720_kernel_characteristic_expressivity_audit.json"
MD = GEN / "p2770_s1720_kernel_characteristic_expressivity_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
K1 = ROOT / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md"
K2 = ROOT / "K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md"
F2 = ROOT / "F2_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET.md"
F3 = ROOT / "F3_CURRENT_FAR_FRONTIER_KERNEL_ARTIFACT_SENSITIVITY_CLASSIFICATION_PACKET.md"
S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"
P2769 = GEN / "p2769_s1719_combined_normalization_orbit_transitivity_no_go.json"

KERNELS = {
    "K_legacy_ont": {
        "formula": "alpha_geo*cos(omega*d+phi)/(1+beta_tors*d)",
        "tokens": {"amplitude", "phase", "oscillation", "linear_damping", "distance_dependence"},
        "parameters": ["alpha_geo", "omega", "phi", "beta_tors", "d"],
    },
    "K_strict_gate": {
        "formula": "cos(omega*d+phi)/(1+beta*d^eta)",
        "tokens": {"phase", "oscillation", "nonlinear_damping", "compression", "distance_dependence"},
        "parameters": ["omega", "phi", "beta", "eta", "d"],
    },
}

CHARACTERISTICS = [
    {
        "id": "single_primordial_foundation",
        "description": "nadsoliton is the sole primordial informational foundation, with no lower information layer",
        "required_formula_tokens": {"ontology_source"},
        "required_theorem_patterns": [r"nadsoliton itself is the primordial information|no separate informational layer|single primordial foundation"],
    },
    {
        "id": "solitonic_state",
        "description": "stable solitonic state or self-maintaining localized carrier",
        "required_formula_tokens": {"oscillation", "damping"},
        "required_theorem_patterns": [r"solitonic state|nadsoliton|source-topology kernel"],
    },
    {
        "id": "self_learning_neural_network_character",
        "description": "self-learning neural-network-like update law, adaptive weights, or loss/feedback dynamics",
        "required_formula_tokens": {"adaptive_update", "feedback_weights", "learning_rule"},
        "required_theorem_patterns": [r"self-learning|neural|Hebbian|adaptive boundary reweight|feedback"],
    },
    {
        "id": "geometric_self_coupling",
        "description": "kernel couples the nadsoliton to itself by an exported geometric self-coupling law",
        "required_formula_tokens": {"self_coupling", "geometric_coupling"},
        "required_theorem_patterns": [r"geometric self-coupling|samosprz|self-coupl|self coupling|coupled to itself"],
    },
    {
        "id": "complete_kernel_expression",
        "description": "one kernel formula fully expresses all above ontology rather than acting as projection/control",
        "required_formula_tokens": {"ontology_source", "adaptive_update", "self_coupling", "geometric_coupling", "completion_bridge"},
        "required_theorem_patterns": [r"completion bridge|bridge-completion|K_legacy_ont -> K_strict_gate|fully express"],
    },
]

NEGATIVE_EXPORT_FLAGS = [
    "kernel_fully_expresses_nadsoliton_characteristics",
    "self_learning_kernel_theorem_exported",
    "geometric_self_coupling_kernel_theorem_exported",
    "legacy_strict_bridge_closed",
    "physical_role_transfer_started",
    "role_bearing_ltotal_promoted",
    "selector_closure_exported",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def run_rg(pattern: str) -> list[str]:
    cmd = ["rg", "-n", "--glob", "!generated/*.json", pattern, "AGENTS.md", "fundamental_action_reconstruction"]
    proc = subprocess.run(cmd, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"rg failed for {pattern!r}: {proc.stderr}")
    return [line for line in proc.stdout.splitlines() if line.strip()]


def normalized_tokens(kernel: dict[str, Any]) -> set[str]:
    tokens = set(kernel["tokens"])
    if "linear_damping" in tokens or "nonlinear_damping" in tokens:
        tokens.add("damping")
    return tokens


def formula_expressivity_matrix() -> dict[str, Any]:
    rows = []
    for characteristic in CHARACTERISTICS:
        row = {
            "characteristic_id": characteristic["id"],
            "description": characteristic["description"],
            "required_formula_tokens": sorted(characteristic["required_formula_tokens"]),
            "kernels": {},
        }
        for name, kernel in KERNELS.items():
            tokens = normalized_tokens(kernel)
            missing = sorted(characteristic["required_formula_tokens"] - tokens)
            present = sorted(characteristic["required_formula_tokens"] & tokens)
            row["kernels"][name] = {
                "formula": kernel["formula"],
                "present_required_tokens": present,
                "missing_required_tokens": missing,
                "formula_level_pass": not missing,
            }
        row["any_kernel_formula_level_pass"] = any(k["formula_level_pass"] for k in row["kernels"].values())
        rows.append(row)
    return {
        "row_count": len(rows),
        "kernel_count": len(KERNELS),
        "rows": rows,
        "all_characteristics_formula_level_pass": all(r["any_kernel_formula_level_pass"] for r in rows),
        "failed_characteristics": [r["characteristic_id"] for r in rows if not r["any_kernel_formula_level_pass"]],
    }


def content_theorem_scan() -> dict[str, Any]:
    rows = []
    for characteristic in CHARACTERISTICS:
        hits_by_pattern = {}
        total = 0
        samples = []
        for pattern in characteristic["required_theorem_patterns"]:
            hits = run_rg(pattern)
            hits_by_pattern[pattern] = len(hits)
            total += len(hits)
            samples.extend(hits[:4])
        rows.append({
            "characteristic_id": characteristic["id"],
            "pattern_hit_count": total,
            "hits_by_pattern": hits_by_pattern,
            "sample_hits": samples[:8],
        })
    closure_hits = run_rg(r"accepted_as_.*(self_learning|geometric_self_coupling|kernel_fully_expresses|bridge_closed).*True|kernel_fully_expresses_nadsoliton_characteristics.*True")
    return {
        "row_count": len(rows),
        "rows": rows,
        "all_characteristics_have_context_hits": all(r["pattern_hit_count"] > 0 for r in rows),
        "theorem_level_closure_hit_count": len(closure_hits),
        "theorem_level_closure_sample_hits": closure_hits[:12],
    }


def finite_feature_rank_witness(matrix: dict[str, Any]) -> dict[str, Any]:
    # Binary coverage over the five characteristic rows.  This finite witness is
    # deliberately simple: it shows formula syntax covers only the solitonic row
    # (oscillation+damping) and cannot encode update/self-coupling/completion data.
    coverage_vectors = {}
    for name in KERNELS:
        vector = []
        for row in matrix["rows"]:
            vector.append(1 if row["kernels"][name]["formula_level_pass"] else 0)
        coverage_vectors[name] = vector
    combined = [max(coverage_vectors[name][i] for name in coverage_vectors) for i in range(len(CHARACTERISTICS))]
    return {
        "characteristic_order": [c["id"] for c in CHARACTERISTICS],
        "coverage_vectors_by_kernel": coverage_vectors,
        "combined_formula_coverage_vector": combined,
        "combined_formula_coverage_count": sum(combined),
        "required_coverage_count": len(CHARACTERISTICS),
        "coverage_defect_count": len(CHARACTERISTICS) - sum(combined),
        "coverage_defect_ids": [c["id"] for c, v in zip(CHARACTERISTICS, combined) if v == 0],
        "finite_conclusion": "The explicit kernel formulas cover oscillation/damping solitonic shape only; they do not encode ontology-source status, learning/update dynamics, geometric self-coupling, or a legacy-to-strict completion theorem.",
    }


def acceptance_matrix(formula_matrix: dict[str, Any], scan: dict[str, Any], witness: dict[str, Any], p2769: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2769_no_closure_boundary_present": p2769.get("status") == "P2769_COMBINED_NORMALIZATION_ORBIT_TRANSITIVITY_NO_GO_NO_CLOSURE",
        "formula_covers_all_required_characteristics": formula_matrix["all_characteristics_formula_level_pass"],
        "formula_coverage_defect_zero": witness["coverage_defect_count"] == 0,
        "theorem_level_closure_export_found": scan["theorem_level_closure_hit_count"] > 0,
        "self_learning_update_token_present_in_kernel": any("adaptive_update" in normalized_tokens(k) for k in KERNELS.values()),
        "geometric_self_coupling_token_present_in_kernel": any("geometric_coupling" in normalized_tokens(k) for k in KERNELS.values()),
        "completion_bridge_token_present_in_kernel": any("completion_bridge" in normalized_tokens(k) for k in KERNELS.values()),
    }
    return {
        "facts": facts,
        "accepted_as_full_kernel_characterization": False,
        "accepted_as_self_learning_kernel_theorem": False,
        "accepted_as_geometric_self_coupling_kernel_theorem": False,
        "accepted_as_legacy_strict_bridge_closure": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "Current kernel formulas are scalar distance kernels with oscillation and damping/compression terms.  They are compatible with a solitonic shape reading, but by formula syntax and current theorem exports they do not fully express self-learning neural update dynamics, geometric self-coupling, sole-foundation ontology, or the legacy-to-strict completion bridge.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["finite_feature_rank_witness"]
    lines = [
        "# P2770/S1720 kernel characteristic expressivity audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite formula coverage witness",
        f"- characteristic_order={witness['characteristic_order']}",
        f"- coverage_vectors_by_kernel={witness['coverage_vectors_by_kernel']}",
        f"- combined_formula_coverage_count={witness['combined_formula_coverage_count']}/{witness['required_coverage_count']}",
        f"- coverage_defect_ids={witness['coverage_defect_ids']}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2769 = read_json(P2769)
    formula_matrix = formula_expressivity_matrix()
    scan = content_theorem_scan()
    witness = finite_feature_rank_witness(formula_matrix)
    acceptance = acceptance_matrix(formula_matrix, scan, witness, p2769)
    payload = {
        "status": "P2770_KERNEL_CHARACTERISTIC_EXPRESSIVITY_AUDIT_NO_FULL_EXPRESSION_NO_CLOSURE",
        "input_hashes": {"K1": sha(K1), "K2": sha(K2), "F2": sha(F2), "F3": sha(F3), "S2": sha(S2), "P2769": sha(P2769)},
        "input_statuses": {"P2769": p2769.get("status")},
        "audited_question": "Do the current kernel formulas fully express the nadsoliton as single primordial information, self-learning neural-network-like system, and geometrically self-coupled sole foundation?",
        "kernels": {name: {**data, "tokens": sorted(data["tokens"])} for name, data in KERNELS.items()},
        "formula_expressivity_matrix": formula_matrix,
        "content_theorem_scan": scan,
        "finite_feature_rank_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not claim that K_legacy_ont or K_strict_gate fully express the self-learning, geometrically self-coupled nadsoliton ontology.  The next proof-grade move should introduce exactly one explicit typed term/theorem: either a geometric self-coupling operator C_geo[K] with a finite invariance/closure test, or a self-learning update law for kernel parameters with source/provenance and a bounded convergence/consistency witness.  After that, rerun this P2770 expressivity matrix as the acceptance test; otherwise preserve the P2697-P2770 no-full-expression/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2770/S1720 kernel characteristic expressivity audit", "## P2770/S1720 kernel characteristic expressivity audit\n\n`P2770/S1720` audits whether the explicit kernel formulas fully express the intended FIN characteristics: single primordial informational nadsoliton, solitonic state, self-learning neural-network-like dynamics, and geometric self-coupling.  The finite formula coverage matrix finds that `K_legacy_ont` and `K_strict_gate` encode oscillation plus damping/compression shape data, so they support only a bounded solitonic-shape reading at formula level.  They do not themselves encode an ontology-source theorem, a self-learning update law, a geometric self-coupling operator, or the completed `K_legacy_ont -> K_strict_gate` bridge.  No physical-role transfer, role-bearing `L_total`, selector closure, bridge closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2770/S1720 kernel-expressivity Ltotal guard", "## P2770/S1720 kernel-expressivity Ltotal guard\n\n`P2770/S1720` adds no variational source term.  It is an expressivity/no-closure audit: the current scalar distance kernels are compatible with oscillatory damped solitonic shape, but do not supply self-learning neural dynamics, geometric self-coupling, or full ontology-to-Lagrangian closure.\n")
    append_once(AGENTS, "Current kernel characteristic expressivity audit guardrail (P2770/S1720, 2026-06-15)", "## Current kernel characteristic expressivity audit guardrail (P2770/S1720, 2026-06-15)\n\n- P2770 audits the user's core foundation question: whether `K_legacy_ont` or `K_strict_gate` fully express the nadsoliton as the single primordial informational, self-learning neural-network-like, geometrically self-coupled foundation.\n- The finite formula-coverage matrix finds only bounded solitonic-shape coverage from oscillation plus damping/compression; the formulas do not themselves encode an ontology-source theorem, self-learning update law, geometric self-coupling operator, or completed `K_legacy_ont -> K_strict_gate` bridge.\n- Do not promote kernel expressivity to physical-role transfer, role-bearing `L_total`, selector closure, bridge closure, or ToE closure.  A next admissible move must introduce exactly one explicit typed term/theorem, preferably a geometric self-coupling operator `C_geo[K]` with finite invariance/closure test or a self-learning kernel-parameter update law with bounded convergence/provenance witness, and then rerun the P2770 acceptance matrix.\n")
    return payload


if __name__ == "__main__":
    main()
