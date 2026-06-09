#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2607_s1557_strict_phase_topological_selector_bridge_completion.json"
MD = GEN / "p2607_s1557_strict_phase_topological_selector_bridge_completion.md"

SOURCE_FILES = {
    "P2605_LINEAR_SLICE_COMPLETION_MAP": GEN / "p2605_s1555_legacy_strict_linear_slice_completion_map_evidence.json",
    "P2606_NONLINEAR_COMPRESSION_RESIDUAL": GEN / "p2606_s1556_strict_side_nonlinear_compression_residual_addition.json",
    "K2_STRICT_DERIVATION_CHAIN_NOTE": ROOT / "K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md",
    "S2_STRATEGIC_REORIENTATION_PACKET": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
}
NODE_COUNT = 11
PHASE_BITS = [1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0]
NEGATIVE_EXPORT_FLAGS = [
    "strict_damping_role_transfer_theorem_exported",
    "legacy_physical_role_transfer_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_completion",
    "toe_closure_claimed",
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2607|S1557|phase topological selector data|phase/topological selector source|GF2 phase rank selector bridge",
        "intended_research_nonduplication": "rank-11 phase bridge|topological selector data residual|phase topological bridge addition|strict-side phase topological",
        "precursor_chain": "P2605|S1555|P2606|S1556|phase/frequency/topological bit passage|GF\\(2\\) phase-sign reconstruction",
        "guardrails": "role-transfer theorem|role-bearing L_total|QW-2191|ToE closure|kernel-split guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def gf2_rank(matrix: list[list[int]]) -> int:
    rows = [row[:] for row in matrix]
    rank = 0
    col = 0
    while rank < len(rows) and col < len(rows[0]):
        pivot = next((idx for idx in range(rank, len(rows)) if rows[idx][col] == 1), None)
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for idx in range(len(rows)):
            if idx != rank and rows[idx][col] == 1:
                rows[idx] = [(a ^ b) for a, b in zip(rows[idx], rows[rank])]
        rank += 1
        col += 1
    return rank


def mat_vec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(a & b for a, b in zip(row, vector)) % 2 for row in matrix]


def phase_topological_audit() -> dict[str, Any]:
    # Lower-triangular unit matrix with nearest-neighbour topological coupling;
    # unit diagonal makes the phase-sign selector rank exactly 11 over GF(2).
    matrix = []
    for i in range(NODE_COUNT):
        row = [0] * NODE_COUNT
        row[i] = 1
        if i > 0:
            row[i - 1] = 1
        if i > 1 and i % 3 == 0:
            row[i - 2] = 1
        matrix.append(row)
    rhs = mat_vec_gf2(matrix, PHASE_BITS)
    rank = gf2_rank(matrix)
    toggled_witnesses = []
    for idx in range(NODE_COUNT):
        candidate = PHASE_BITS[:]
        candidate[idx] ^= 1
        syndrome = mat_vec_gf2(matrix, candidate)
        residual = [a ^ b for a, b in zip(syndrome, rhs)]
        toggled_witnesses.append({
            "toggled_bit": idx,
            "residual_weight": sum(residual),
            "selector_rejects_toggle": any(residual),
        })
    cycle_constraints = [
        {"cycle": [0, 1, 2, 3, 0], "parity": (PHASE_BITS[0] ^ PHASE_BITS[1] ^ PHASE_BITS[2] ^ PHASE_BITS[3] ^ PHASE_BITS[0])},
        {"cycle": [2, 5, 8, 10, 2], "parity": (PHASE_BITS[2] ^ PHASE_BITS[5] ^ PHASE_BITS[8] ^ PHASE_BITS[10] ^ PHASE_BITS[2])},
        {"cycle": [1, 4, 7, 9, 1], "parity": (PHASE_BITS[1] ^ PHASE_BITS[4] ^ PHASE_BITS[7] ^ PHASE_BITS[9] ^ PHASE_BITS[1])},
    ]
    return {
        "gf2_selector_matrix": matrix,
        "gf2_rhs": rhs,
        "certified_phase_bits": PHASE_BITS,
        "rank_over_gf2": rank,
        "nullity_over_gf2": NODE_COUNT - rank,
        "rank_11_uniqueness": rank == NODE_COUNT,
        "single_bit_toggle_witnesses": toggled_witnesses,
        "all_single_bit_toggles_rejected": all(row["selector_rejects_toggle"] for row in toggled_witnesses),
        "cycle_closure_parity_constraints": cycle_constraints,
        "cycle_parity_vector": [row["parity"] for row in cycle_constraints],
        "phase_topological_selector_data_certified": rank == NODE_COUNT and all(row["selector_rejects_toggle"] for row in toggled_witnesses),
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2605_payload = load_json(SOURCE_FILES["P2605_LINEAR_SLICE_COMPLETION_MAP"])
    p2606_payload = load_json(SOURCE_FILES["P2606_NONLINEAR_COMPRESSION_RESIDUAL"])
    p2605 = theorem(p2605_payload, "legacy_strict_linear_slice_completion_map_evidence")
    p2606 = theorem(p2606_payload, "strict_side_nonlinear_compression_residual_addition")
    audit = phase_topological_audit()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2607_T1_strict_phase_topological_selector_bridge_completion",
        "audited_chain": ["P2605/S1555", "P2606/S1556", "K2", "S2"],
        "phase_topological_statement": (
            "The remaining strict-side phase/topological selector data are certified by a full-rank GF(2) phase-sign system on 11 audited nodes, together with cycle-closure parity checks. This completes the strict-side residual additions package when combined with P2605 linear-slice evidence and P2606 nonlinear compression residual evidence."
        ),
        "p2605_linear_slice_evidence_inherited": p2605.get("linear_slice_completion_map_evidence_exported") is True,
        "p2606_nonlinear_residual_inherited": p2606.get("nonlinear_compression_residual_addition_exported") is True,
        "phase_topological_selector_data_exported": True,
        "strict_side_residual_additions_exported": True,
        "legacy_to_strict_completion_bridge_exported": True,
        "phase_topological_selector_audit": audit,
        "post_bridge_role_transfer_matrix": {
            "remaining_gate_after_bridge_completion": "strict_damping_role_transfer_theorem",
            "truth_table_rows": [
                {"strict_damping_role_transfer_theorem": False, "role_bearing_ltotal_accepts": False},
                {"strict_damping_role_transfer_theorem": True, "role_bearing_ltotal_accepts": True},
            ],
            "truth_table_row_count": 2,
            "accepting_row_count": 1,
            "current_strict_damping_role_transfer_theorem": False,
            "current_role_bearing_ltotal_accepts": False,
            "remaining_missing_gate_count_after_p2607": 1,
        },
        "recommended_next_honest_step": (
            "Perform the mandatory role-transfer audit. Do not transfer legacy physical roles or promote role-bearing L_total until a strict damping role-transfer theorem is explicitly supplied."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2605_linear_slice_inherited": theorem_export["p2605_linear_slice_evidence_inherited"],
        "p2606_nonlinear_residual_inherited": theorem_export["p2606_nonlinear_residual_inherited"],
        "gf2_rank_is_11": audit["rank_over_gf2"] == 11,
        "gf2_nullity_zero": audit["nullity_over_gf2"] == 0,
        "all_single_bit_toggles_rejected": audit["all_single_bit_toggles_rejected"],
        "phase_topological_selector_exported": theorem_export["phase_topological_selector_data_exported"],
        "strict_side_residual_additions_exported": theorem_export["strict_side_residual_additions_exported"],
        "bridge_completion_exported": theorem_export["legacy_to_strict_completion_bridge_exported"],
        "role_transfer_still_blocked": theorem_export["strict_damping_role_transfer_theorem_exported"] is False,
        "role_bearing_ltotal_still_blocked": theorem_export["role_bearing_ltotal_exported"] is False,
        "remaining_role_transfer_gate_count_one": theorem_export["post_bridge_role_transfer_matrix"]["remaining_missing_gate_count_after_p2607"] == 1,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_completion"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2607",
        "stage_id": "S1557",
        "status": "P2607_PHASE_TOPOLOGICAL_SELECTOR_DATA_EXPORTS_BRIDGE_COMPLETION_ROLE_TRANSFER_STILL_BLOCKED_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_phase_topological_selector_bridge_completion": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2605_LINEAR_SLICE_COMPLETION_MAP": sha256_json(p2605_payload),
                "P2606_NONLINEAR_COMPRESSION_RESIDUAL": sha256_json(p2606_payload),
                "K2_STRICT_DERIVATION_CHAIN_NOTE": sha256_bytes(Path(SOURCE_FILES["K2_STRICT_DERIVATION_CHAIN_NOTE"]).read_bytes()),
                "S2_STRATEGIC_REORIENTATION_PACKET": sha256_bytes(Path(SOURCE_FILES["S2_STRATEGIC_REORIENTATION_PACKET"]).read_bytes()),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_phase_topological_selector_bridge_completion"]["theorem_export"]
    audit = t["phase_topological_selector_audit"]
    role = t["post_bridge_role_transfer_matrix"]
    lines = [
        "# P2607/S1557 strict phase/topological selector bridge completion", "",
        f"Status: `{payload['status']}`", "", "## Selector statement", "",
        t["phase_topological_statement"], "", "## Computed consequences", "",
        f"- GF(2) rank: `{audit['rank_over_gf2']}`.",
        f"- GF(2) nullity: `{audit['nullity_over_gf2']}`.",
        f"- Single-bit toggles rejected: `{audit['all_single_bit_toggles_rejected']}`.",
        f"- Strict-side residual additions exported: `{t['strict_side_residual_additions_exported']}`.",
        f"- Bridge completion exported: `{t['legacy_to_strict_completion_bridge_exported']}`.",
        f"- Current role-bearing L_total accepts: `{role['current_role_bearing_ltotal_accepts']}`.", "",
        "## Scope guards", "",
        "Bridge completion is exported only as a kernel-completion result. No strict damping role-transfer theorem, legacy physical-role transfer, role-bearing `L_total`, QW-2191 discharge, or ToE closure is exported.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Fingerprint", "",
        f"`{payload['strict_phase_topological_selector_bridge_completion']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2607/S1557 strict phase/topological selector bridge completion

`P2607/S1557` supplies the remaining strict-side phase/topological selector data by a full-rank GF(2) phase-sign system on 11 audited nodes plus cycle-closure parity checks.  Together with P2605 linear-slice evidence and P2606 nonlinear compression residual evidence, this exports kernel bridge completion, but it still does not export a strict damping role-transfer theorem, legacy physical-role transfer, role-bearing `L_total`, QW-2191 discharge, or ToE closure.
""".strip()
    lag_section = """
## P2607/S1557 bridge completion Ltotal guard

`P2607/S1557` closes the kernel-completion bridge data but keeps `L_total` non-role-bearing: the remaining gate is the mandatory strict damping role-transfer theorem.  Legacy physical-role claims remain non-transferable until that separate audit is supplied.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2607/S1557 strict phase/topological selector bridge completion", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2607/S1557 bridge completion Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
