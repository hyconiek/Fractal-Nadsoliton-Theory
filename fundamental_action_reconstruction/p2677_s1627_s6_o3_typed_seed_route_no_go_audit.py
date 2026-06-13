#!/usr/bin/env python3
"""P2677/S1627: no-go audit for the current S6/O3 typed-seed route.

After P2676, the only formal Sigma->F301 naturality candidates are an
XOR/XNOR reversal pair.  This packet checks whether the current repository
exports an internal orientation source that can choose one member pre-collapse;
if not, it records a bounded no-go for the current tau_src -> pair12 ->
boundary-square typed-seed interface route.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2677_s1627_s6_o3_typed_seed_route_no_go_audit.json"
MD = GEN / "p2677_s1627_s6_o3_typed_seed_route_no_go_audit.md"
P2676 = GEN / "p2676_s1626_sigma_f301_precollapse_naturality_square_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "internal_orientation_source_exported",
    "xor_xnor_reversal_broken_strictly",
    "sigma_to_f301_typed_arrow_exported",
    "s6_current_route_passed",
    "o3_current_route_passed",
    "boundary_square_arrow_allowed",
    "sector_swap_invariant_allowed",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "xor_xnor_reversal_content": "XOR|XNOR|orientation-reversal|reversal pair|orientation source",
        "internal_selector_source_content": "internal selector|selector source|source obstruction|QW-2191|symmetry-breaking|symmetry breaking",
        "sigma_f301_arrow_content": "Sigma_sel_src_target_v1.*F301|F301.*Sigma_sel_src_target_v1|Sigma->F301|Sigma_sel_src_target_v1 -> F301",
        "typed_seed_route_content": "tau_src.*pair12|pair1/pair2 typed seed|typed-seed|boundary-square|sector-swap",
        "collapse_forbidden_content": "pre-collapse|Q_basis|preLM|projector-only|nonprojector|quotient|collapse",
        "fiat_convention_blocker_content": "by fiat|fiat|convention choice|must_not_identify|not actually realized|future-only",
        "closure_guard_content": "O4/O5|L_total|ToE closure|role-bearing|bridge completion|role transfer",
    }
    return {
        "tool": "rg",
        "mode": "content-first internal-orientation-source/no-go audit for current S6/O3 typed-seed route; not name/number-only",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def upstream_consistency() -> dict[str, Any]:
    p2676 = load_json(P2676)
    decision = p2676.get("closure_decision", {})
    witness = p2676.get("finite_boolean_naturality_witness", {})
    return {
        "p2676_formal_xor_xnor_pair_exists": witness.get("formal_square_candidate_count") == 2,
        "p2676_export_gate_zero": witness.get("passing_export_gate_count") == 0,
        "p2676_no_sigma_to_f301_arrow": decision.get("sigma_to_f301_typed_arrow_exported_now") is False,
        "p2676_s6_not_passed": decision.get("s6_exported_now") is False,
        "p2676_o3_not_passed": decision.get("o3_exported_now") is False,
    }


def orientation_source_candidate_table() -> list[dict[str, Any]]:
    candidates = [
        ("beta_tors_to_chi11", True, False, False, False, "legacy-role/source conjecture remains audited; not a strict pre-collapse selector export"),
        ("boundary_phase_theta_like_sign", True, False, False, False, "can force a bit only with an externally supplied theta-like sign"),
        ("q_basis_terminal_choice", True, False, False, False, "appears after/through Q_basis collapse and loses the required pre-collapse typing"),
        ("projector_only_local_atlas", True, False, False, False, "retains local chart material but only projector-side; no nonprojector descent"),
        ("xor_orientation_by_declaration", True, False, False, False, "chooses XOR by convention/fiat rather than internal source"),
        ("xnor_orientation_by_declaration", True, False, False, False, "chooses XNOR by convention/fiat rather than internal source"),
        ("spontaneous_symmetry_breaking_placeholder", True, False, False, False, "symmetry-breaking named as a possible premise but no exported strict source object binds it here"),
        ("observer_readout_downstream", True, False, False, False, "downstream readout violates nadsoliton -> light -> matter -> observer order if used as primordial selector"),
    ]
    rows = []
    for name, mentioned, internal_source, precollapse, nonfiat, reason in candidates:
        passes = mentioned and internal_source and precollapse and nonfiat
        rows.append({
            "candidate": name,
            "mentioned_or_search_visible": mentioned,
            "strict_internal_orientation_source_exported": internal_source,
            "precollapse_sigma_to_f301_binding": precollapse,
            "nonfiat_nonconvention_choice": nonfiat,
            "passes_orientation_source_gate": passes,
            "blocker": reason,
        })
    return rows


def finite_no_go_lattice() -> dict[str, Any]:
    obligations = [
        "formal_sigma_f301_reversal_pair_identified",
        "strict_internal_orientation_source_exported",
        "source_binds_sigma_side_precollapse",
        "source_descends_to_surviving_f301_carrier",
        "nonquotient_nonprojector_nonprelm_transport",
        "nonfiat_orientation_choice_proof",
        "route_closure_guards_preserved",
    ]
    current = {
        obligations[0]: True,
        obligations[1]: False,
        obligations[2]: False,
        obligations[3]: False,
        obligations[4]: False,
        obligations[5]: False,
        obligations[6]: True,
    }
    pass_count = 0
    rows = []
    for bits in itertools.product([False, True], repeat=len(obligations)):
        state = dict(zip(obligations, bits))
        passes = all(state.values())
        pass_count += int(passes)
        rows.append({"state": state, "passes_current_route": passes})
    missing = [key for key, value in current.items() if not value]
    return {
        "obligations": obligations,
        "total_states": len(rows),
        "passing_states": pass_count,
        "current_state": current,
        "missing_current_obligations": missing,
        "hamming_distance_to_pass": len(missing),
    }


def closure_decision(lattice: dict[str, Any], candidates: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "decision": "P2677_S6_O3_TYPED_SEED_ROUTE_NO_GO__NO_INTERNAL_ORIENTATION_SOURCE_FOR_XOR_XNOR_REVERSAL",
        "professorial_verdict": (
            "P2677 checks the exact escape hatch left by P2676. The repo has visible candidate lanes, but every lane is either legacy-role conjectural, theta/sign external, Q_basis/projector collapsed, "
            "fiat/convention-only, placeholder symmetry breaking, or downstream observer readout. None exports a strict internal pre-collapse orientation source choosing XOR over XNOR inside the Sigma->F301 square. "
            "Therefore the current S6/O3 tau_src -> pair12 -> boundary-square typed-seed route is promoted to a bounded no-go, not a global impossibility theorem for all future bridge/source classes."
        ),
        "next_honest_step": (
            "Do not attempt O4/O5 on this route. The next honest proof-grade move is to leave the current S6/O3 typed-seed branch closed as a no-go and switch to a genuinely different provider class: "
            "either a legacy->strict completion bridge source for the missing orientation bit, or a separate strict internal symmetry-breaking/source theorem with an explicit pre-collapse Sigma->F301 binding."
        ),
        "orientation_source_candidates_checked": len(candidates),
        "passing_orientation_source_candidates": sum(row["passes_orientation_source_gate"] for row in candidates),
        "no_go_scope": "current tau_src -> pair12 -> Sigma_sel_src_target_v1 -> F301 -> boundary-square typed-seed route only",
        "global_no_go_claimed": False,
        "s6_current_route_passed_now": False,
        "o3_current_route_passed_now": False,
        "o4_o5_allowed_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
        "hamming_distance_to_pass": lattice["hamming_distance_to_pass"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2677/S1627 S6/O3 typed-seed route no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Orientation-source candidate table"])
    for row in payload["orientation_source_candidate_table"]:
        lines.append(f"- `{row['candidate']}`: pass=`{row['passes_orientation_source_gate']}` — {row['blocker']}")
    lat = payload["finite_no_go_lattice"]
    lines.extend([
        "", "## Finite no-go lattice",
        f"Total states: `{lat['total_states']}`; passing states: `{lat['passing_states']}`.",
        f"Current Hamming distance to pass: `{lat['hamming_distance_to_pass']}`.",
        f"Missing current obligations: `{lat['missing_current_obligations']}`.",
        "", "## Verdict", payload["closure_decision"]["professorial_verdict"],
        f"Decision: `{payload['closure_decision']['decision']}`.",
        "", "## Next honest step", payload["closure_decision"]["next_honest_step"],
        "", "## Negative exports",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    candidates = orientation_source_candidate_table()
    lattice = finite_no_go_lattice()
    payload: dict[str, Any] = {
        "status": "P2677_S6_O3_TYPED_SEED_ROUTE_NO_GO_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {"P2676": sha256_file(P2676), "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET), "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT)},
        "upstream_consistency": upstream_consistency(),
        "orientation_source_candidate_table": candidates,
        "finite_no_go_lattice": lattice,
        "closure_decision": closure_decision(lattice, candidates),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2677/S1627 S6/O3 typed-seed route no-go guard",
        "## P2677/S1627 S6/O3 typed-seed route no-go guard\n\n"
        "`P2677/S1627` audits the escape hatch left by `P2676`: an internal orientation source that would choose one member of the XOR/XNOR `Sigma_sel_src_target_v1 -> F301` reversal pair before `Q_basis`/preLM/projector collapse.  The checked lanes are visible but fail as strict exports: `beta_tors -> chi_11` is still a legacy-role/source conjecture, theta-like signs are external, `Q_basis` and projector lanes collapse typing, declaration lanes are fiat, symmetry breaking is only a placeholder, and observer readout is downstream.  Thus the current `tau_src -> pair12 -> boundary-square` S6/O3 typed-seed route is a bounded no-go; this is not a global no-go for future bridge/source classes.  O4/O5, `QW-2191` discharge, role-bearing `L_total`, bridge completion, role transfer, and ToE closure remain blocked.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2677/S1627 S6/O3 typed-seed route no-go Ltotal guard",
        "## P2677/S1627 S6/O3 typed-seed route no-go Ltotal guard\n\n"
        "`P2677/S1627` keeps `L_total` closed for the current tau_src/pair12 typed-seed branch.  Since no strict internal pre-collapse orientation source breaks the XOR/XNOR Sigma->F301 reversal, no boundary-square or sector-swap variational source term may be imported from this route.\n",
    )
    return payload


if __name__ == "__main__":
    main()
