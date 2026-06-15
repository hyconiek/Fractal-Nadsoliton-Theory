#!/usr/bin/env python3
"""P2760/S1710: foundation -> kernel/coupling/Lagrangian gap matrix.

This audit answers the professorial question: where do the theory foundations
(nadsoliton ontology, legacy/strict kernel split, selector guardrail) fail to be
exported as mathematically licensed formulas in the coupling kernels and
Lagrangian/EOM chain?  It combines content-first provenance checks with a small
finite numerical kernel comparison and a formula-obligation matrix.  It does
not claim bridge closure, role transfer, selector closure, L_total closure, or
ToE closure.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
P1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"
P1866 = GEN / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json"
P2759 = GEN / "p2759_s1709_post_p2758_no_new_live_frontier_reconciliation.json"
OUT = GEN / "p2760_s1710_foundation_kernel_lagrangian_gap_matrix.json"
MD = GEN / "p2760_s1710_foundation_kernel_lagrangian_gap_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "kernel_split_foundation": r"K_legacy_ont|K_strict_gate|legacy kernel|strict gate kernel|kernel split",
    "ontology_foundation": r"nadsoliton itself|primordial information|no separate informational layer|nadsoliton -> light -> matter",
    "bridge_obligations": r"completion bridge|amplitude/normalization|phase/frequency|damping/compression|role-transfer theorem",
    "selector_boundary": r"QW-2191|selector closure|orientation/polarity|P2721|lambda/P2721",
    "lagrangian_chain": r"L_total|Lagrangian|Euler-Lagrange|effective couplings|kernel_strict -> coefficients",
    "no_toe_boundary": r"ToE closure|toe_closed|role-bearing L_total|no-new-live-frontier",
}

NEGATIVE_EXPORT_FLAGS = [
    "legacy_strict_bridge_closed",
    "legacy_role_transfer_started",
    "selector_source_exported",
    "qw2191_discharged",
    "ltotal_promoted",
    "toe_closure_exported",
    "foundation_to_lagrangian_gap_closed",
]

LEGACY = {
    "alpha_geo": 4.0 * math.log(2.0),
    "omega": math.pi / 4.0,
    "phi": math.pi / 6.0,
    "beta_tors": 0.01,
}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}


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


def evidence_scan() -> dict[str, Any]:
    rows = []
    for name, pattern in CONTENT_PATTERNS.items():
        hits = run_rg(pattern)
        rows.append({"lane": name, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return {"row_count": len(rows), "rows": rows, "hit_counts": {r["lane"]: r["hit_count"] for r in rows}, "all_patterns_have_hits": all(r["hit_count"] > 0 for r in rows)}


def k_legacy(d: float) -> float:
    return LEGACY["alpha_geo"] * math.cos(LEGACY["omega"] * d + LEGACY["phi"]) / (1.0 + LEGACY["beta_tors"] * d)


def k_strict(d: float) -> float:
    return math.cos(STRICT["omega"] * d + STRICT["phi"]) / (1.0 + STRICT["beta"] * (d**STRICT["eta"]))


def finite_kernel_comparison() -> dict[str, Any]:
    rows = []
    max_abs_delta = 0.0
    max_abs_delta_d = None
    sign_mismatch_count = 0
    for d in range(0, 13):
        legacy = k_legacy(float(d))
        strict = k_strict(float(d)) if d != 0 else math.cos(STRICT["phi"])
        delta = legacy - strict
        max_abs_delta = max(max_abs_delta, abs(delta))
        if abs(delta) == max_abs_delta:
            max_abs_delta_d = d
        sign_mismatch = (legacy > 0) != (strict > 0) if legacy != 0 and strict != 0 else legacy != strict
        sign_mismatch_count += int(sign_mismatch)
        rows.append({"d": d, "K_legacy_ont": legacy, "K_strict_gate": strict, "delta": delta, "sign_mismatch": sign_mismatch})
    return {
        "legacy_formula": "alpha_geo*cos((pi/4)*d+pi/6)/(1+0.01*d)",
        "strict_formula": "cos(0.18575*d+0.16250)/(1+1.0*d^1.8)",
        "sample_domain": list(range(13)),
        "rows": rows,
        "max_abs_delta": max_abs_delta,
        "max_abs_delta_at_d": max_abs_delta_d,
        "sign_mismatch_count": sign_mismatch_count,
        "amplitude_ratio_at_d0_legacy_over_strict": rows[0]["K_legacy_ont"] / rows[0]["K_strict_gate"],
        "same_formula_on_sample_domain": all(abs(r["delta"]) < 1e-12 for r in rows),
    }


def formula_gap_matrix(p1562: dict[str, Any], p1563: dict[str, Any], p1866: dict[str, Any]) -> dict[str, Any]:
    p1562_claims_closure = bool(p1562.get("qw2191_closed")) or bool(p1562.get("toe_closed"))
    p1563_blocks_closure = (p1563.get("qw2191_closed") is False) and (p1563.get("toe_closed") is False)
    p1866_open = p1866.get("status") == "OPEN_OBSTRUCTION_WITH_TRACE"
    rows = [
        {
            "gap_id": "G1_ontology_to_kernel_measure",
            "foundation_obligation": "derive the shell/distance measure d and kernel carrier from nadsoliton-as-primordial-information, not as a fitted later gate coordinate",
            "current_formula_layer": "K_legacy_ont(d) and K_strict_gate(d) both use d, but K2/F2 classify the strict tuple as later operational gate data",
            "blocking_evidence": "No bridge theorem identifies the strict gate coordinate with the ontological source-layer distance.",
            "closed": False,
        },
        {
            "gap_id": "G2_amplitude_normalization_passage",
            "foundation_obligation": "map alpha_geo=4 ln 2 and any amplitude absorption into strict normalization without role transfer",
            "current_formula_layer": "K_legacy_ont carries alpha_geo multiplicatively; K_strict_gate has no alpha_geo factor and P2754-P2758 do not export selector authority from entropy",
            "blocking_evidence": "No role-safe amplitude absorption/source theorem is exported.",
            "closed": False,
        },
        {
            "gap_id": "G3_phase_frequency_topological_source",
            "foundation_obligation": "source omega, phi, and the certified phase/topological data from a non-premise selector/orientation theorem",
            "current_formula_layer": "legacy uses pi/4, pi/6; strict uses 0.18575, 0.16250 from gate/refreeze chain",
            "blocking_evidence": "QW-2191/P2721/P2759 still block non-premise selector or polarity closure.",
            "closed": False,
        },
        {
            "gap_id": "G4_damping_compression_bridge",
            "foundation_obligation": "derive nonlinear strict beta*d^eta compression from legacy beta_tors*d and specify residual strict-side additions",
            "current_formula_layer": "legacy denominator is 1+beta_tors*d; strict denominator is 1+beta*d^eta",
            "blocking_evidence": "No target-independent positive beta/Z_beta source theorem or beta_tors->strict compression map is exported.",
            "closed": False,
        },
        {
            "gap_id": "G5_kernel_moments_to_physical_couplings",
            "foundation_obligation": "prove that moment-derived lambda_sm_eff/kappa_gr_eff/epsilon_mix_eff are physical couplings with units, sign conventions, and variational provenance",
            "current_formula_layer": "P1562 maps strict moments to effective couplings; P1563/P1866 encode Lagrangian/EOM candidates",
            "blocking_evidence": "P1563 lists missing selector/full proof exports; P1866 remains OPEN_OBSTRUCTION_WITH_TRACE.",
            "closed": False,
        },
        {
            "gap_id": "G6_lagrangian_reverse_closure",
            "foundation_obligation": "show EOM/residual tables reverse-constrain the admissible kernel tuple and close full nonproxy variational equations",
            "current_formula_layer": "P1866 gives symbolic L_total decomposition and proxy EOM exports",
            "blocking_evidence": "P2685-P2687/P1866 require tensor residual/source classes; L_total is not role-bearing closure.",
            "closed": False,
        },
        {
            "gap_id": "G7_closure_flag_consistency",
            "foundation_obligation": "prevent stale closure flags from overriding later no-go guardrails",
            "current_formula_layer": "P1562 has qw2191_closed/toe_closed true while P1563 and current guardrails keep both false/open",
            "blocking_evidence": "Treat P1562 closure booleans as stale/overruled by P1563/P1866/P2759 unless revalidated by a new theorem.",
            "closed": False,
            "stale_flag_detected": p1562_claims_closure and p1563_blocks_closure and p1866_open,
        },
    ]
    return {
        "rows": rows,
        "gap_count": len(rows),
        "closed_gap_count": sum(1 for row in rows if row["closed"]),
        "open_gap_count": sum(1 for row in rows if not row["closed"]),
        "stale_closure_flag_detected": any(row.get("stale_flag_detected", False) for row in rows),
    }


def acceptance_matrix(scan: dict[str, Any], kernel_cmp: dict[str, Any], gaps: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "foundation_and_formula_evidence_present": scan["all_patterns_have_hits"],
        "legacy_and_strict_kernels_differ_on_finite_sample": not kernel_cmp["same_formula_on_sample_domain"],
        "gap_matrix_has_open_gaps": gaps["open_gap_count"] > 0,
        "stale_closure_flag_detected_and_quarantined": gaps["stale_closure_flag_detected"],
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_or_p2721_source_exported": False,
        "ltotal_reverse_closure_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_foundation_to_lagrangian_closure": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The finite audit finds real formula content, but the bridge from ontology to strict kernel, amplitude/phase/damping passage, selector source, physical coupling provenance, and Lagrangian reverse closure remain open.  Stale P1562 closure flags are inconsistent with later P1563/P1866/P2759 boundaries and must not promote closure.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cmp = payload["finite_kernel_comparison"]
    gaps = payload["formula_gap_matrix"]
    lines = [
        "# P2760/S1710 foundation-kernel-Lagrangian gap matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite kernel comparison",
        f"- legacy_formula={cmp['legacy_formula']}",
        f"- strict_formula={cmp['strict_formula']}",
        f"- sample_domain={cmp['sample_domain']}",
        f"- same_formula_on_sample_domain={cmp['same_formula_on_sample_domain']}",
        f"- max_abs_delta={cmp['max_abs_delta']}",
        f"- max_abs_delta_at_d={cmp['max_abs_delta_at_d']}",
        f"- sign_mismatch_count={cmp['sign_mismatch_count']}",
        f"- amplitude_ratio_at_d0_legacy_over_strict={cmp['amplitude_ratio_at_d0_legacy_over_strict']}",
        "",
        "## Open gap matrix",
        f"- gap_count={gaps['gap_count']}",
        f"- open_gap_count={gaps['open_gap_count']}",
        f"- stale_closure_flag_detected={gaps['stale_closure_flag_detected']}",
        "",
    ]
    for row in gaps["rows"]:
        lines.append(f"- {row['gap_id']}: closed={row['closed']} — {row['blocking_evidence']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p1562 = read_json(P1562)
    p1563 = read_json(P1563)
    p1866 = read_json(P1866)
    p2759 = read_json(P2759)
    scan = evidence_scan()
    kernel_cmp = finite_kernel_comparison()
    gaps = formula_gap_matrix(p1562, p1563, p1866)
    acceptance = acceptance_matrix(scan, kernel_cmp, gaps)
    payload = {
        "status": "P2760_FOUNDATION_KERNEL_LAGRANGIAN_GAP_MATRIX_NO_CLOSURE",
        "input_hashes": {
            "P1562": sha(P1562),
            "P1563": sha(P1563),
            "P1866": sha(P1866),
            "P2759": sha(P2759),
        },
        "input_statuses": {
            "P1562": p1562.get("status"),
            "P1563": p1563.get("status"),
            "P1866": p1866.get("status"),
            "P2759": p2759.get("status"),
        },
        "content_evidence_scan": scan,
        "finite_kernel_comparison": kernel_cmp,
        "formula_gap_matrix": gaps,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote current kernel-derived coupling/Lagrangian formulas to foundation-level closure.  The next proof-grade move should target exactly one gap with a machine-checkable theorem: first choice is G5, the kernel-moments-to-physical-couplings provenance theorem with units/sign/variational normalization and explicit quarantine of stale P1562 closure flags; absent that theorem, preserve the P2697-P2760 no-new-live-frontier/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2760/S1710 foundation-kernel-Lagrangian gap matrix", "## P2760/S1710 foundation-kernel-Lagrangian gap matrix\n\n`P2760/S1710` audits the theoretical foundations against their formula-level expression in kernels, effective couplings, and the Lagrangian/EOM chain.  On the finite sample `d=0..12`, `K_legacy_ont(d)=alpha_geo*cos((pi/4)d+pi/6)/(1+0.01d)` and `K_strict_gate(d)=cos(0.18575d+0.16250)/(1+d^1.8)` are not the same formula; the maximum sampled absolute delta is nonzero and the amplitude/phase/damping structures differ.  The gap matrix leaves seven open obligations: ontology-to-kernel measure, amplitude normalization passage, phase/frequency/topological source, damping/compression bridge, kernel-moments-to-physical-couplings provenance, Lagrangian reverse closure, and stale closure-flag consistency.  In particular, P1562 stale `qw2191_closed/toe_closed` booleans are quarantined by P1563/P1866/current guardrails.  No bridge closure, role transfer, selector closure, role-bearing `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2760/S1710 foundation-kernel-Lagrangian gap Ltotal guard", "## P2760/S1710 foundation-kernel-Lagrangian gap Ltotal guard\n\n`P2760/S1710` is a gap matrix, not a new variational source term.  It confirms that kernel-derived effective coupling and Lagrangian formulas remain blocked by missing ontology-to-kernel, amplitude/phase/damping bridge, physical-coupling provenance, selector/P2721, and reverse-closure theorems.  It also quarantines stale P1562 closure flags against later P1563/P1866/current guardrails.  Therefore it cannot promote role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.\n")
    append_once(AGENTS, "Current foundation-kernel-Lagrangian gap matrix guardrail (P2760/S1710, 2026-06-15)", "## Current foundation-kernel-Lagrangian gap matrix guardrail (P2760/S1710, 2026-06-15)\n\n- P2760 audits the foundations against the explicit kernel/coupling/Lagrangian formulas and finds seven open gaps: ontology-to-kernel measure, amplitude normalization, phase/frequency/topological source, damping/compression bridge, kernel-moments-to-physical-couplings provenance, Lagrangian reverse closure, and stale closure-flag consistency.\n- The finite comparison confirms that `K_legacy_ont` and `K_strict_gate` remain formula-distinct on `d=0..12`; P1562 stale `qw2191_closed/toe_closed` flags are quarantined by P1563/P1866/current guardrails.\n- Do not promote kernel-derived effective couplings or `L_total` to foundation-level closure, role transfer, selector closure, bridge closure, or ToE closure until one named gap is closed by a machine-checkable theorem.  The next preferred target is G5: kernel-moments-to-physical-couplings provenance with units/sign/variational normalization.\n")
    return payload


if __name__ == "__main__":
    main()
