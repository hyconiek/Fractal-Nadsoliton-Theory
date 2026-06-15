#!/usr/bin/env python3
"""P2775/S1725: full-Laplacian-spectrum pair discriminator.

P2774 left a concrete degeneracy: H=ln(16)=4 ln 2 plus equal 4-regular
Laplacian trace cannot distinguish torus_4x4 from circulant_pm1_pm2.  This
bounded follow-up tests exactly one stronger metric principle recommended there:
the full graph-Laplacian spectrum.  The result is a finite positive discriminator
for the P2774 pair, but not a theorem-level canonical nadsoliton geometry source:
it is only a pair-local spectral witness and has no exported ontological/source
law selecting why this spectrum, graph class, or target geometry is mandatory.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2774_s1724_entropy_laplacian_trace_geometry_degeneracy import graph_summary, graph_edges, N, ALPHA_GEO

GEN = ROOT / "generated"
P2774 = GEN / "p2774_s1724_entropy_laplacian_trace_geometry_degeneracy.json"
OUT = GEN / "p2775_s1725_full_laplacian_spectrum_pair_discriminator.json"
MD = GEN / "p2775_s1725_full_laplacian_spectrum_pair_discriminator.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NEGATIVE_EXPORT_FLAGS = [
    "canonical_geometry_source_exported",
    "full_spectrum_nadsoliton_geometry_theorem_exported",
    "kernel_geometry_closure_exported",
    "kernel_fully_expresses_nadsoliton_characteristics",
    "role_bearing_ltotal_promoted",
    "bridge_closure_exported",
    "selector_closure_exported",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def laplacian_matrix(kind: str) -> np.ndarray:
    mat = np.zeros((N, N), dtype=float)
    for a, b in graph_edges(kind):
        mat[a, b] = mat[b, a] = -1.0
        mat[a, a] += 1.0
        mat[b, b] += 1.0
    return mat


def rounded_spectrum(kind: str, places: int = 10) -> list[float]:
    eigs = np.linalg.eigvalsh(laplacian_matrix(kind))
    return [round(float(x), places) for x in eigs]


def multiplicities(values: list[float]) -> dict[str, int]:
    out: dict[str, int] = {}
    for value in values:
        key = f"{value:.10f}"
        out[key] = out.get(key, 0) + 1
    return out


def spectral_witness() -> dict[str, Any]:
    kinds = ["torus_4x4", "circulant_pm1_pm2"]
    rows = []
    for kind in kinds:
        summary = graph_summary(kind)
        spectrum = rounded_spectrum(kind)
        rows.append({
            "geometry": kind,
            "node_count": N,
            "edge_count": summary["edge_count"],
            "laplacian_trace": summary["laplacian_trace"],
            "distance_histogram": summary["distance_histogram"],
            "laplacian_spectrum_rounded": spectrum,
            "laplacian_spectrum_multiplicities": multiplicities(spectrum),
            "spectral_sum": round(sum(spectrum), 10),
            "spectral_zero_multiplicity": sum(1 for value in spectrum if abs(value) < 1e-10),
        })
    signatures = {json.dumps(row["laplacian_spectrum_rounded"]) for row in rows}
    return {
        "entropy_object": "uniform distribution on 16 support points",
        "shannon_entropy_nats": math.log(N),
        "alpha_geo_4_ln_2": ALPHA_GEO,
        "entropy_matches_alpha_geo": abs(math.log(N) - ALPHA_GEO) < 1e-12,
        "candidate_extra_principle": "full graph-Laplacian spectrum on the P2774 pair",
        "geometry_rows": rows,
        "geometry_row_count": len(rows),
        "same_laplacian_trace": len({row["laplacian_trace"] for row in rows}) == 1,
        "distinct_full_laplacian_spectrum_count": len(signatures),
        "full_spectrum_breaks_p2774_pair_degeneracy": len(signatures) == len(rows),
        "finite_positive_statement": "The full Laplacian spectrum distinguishes the P2774 torus_4x4 and circulant_pm1_pm2 pair even though entropy and Laplacian trace do not.",
    }


def acceptance_matrix(witness: dict[str, Any], p2774: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2774_entropy_laplacian_trace_degeneracy_present": p2774.get("status") == "P2774_ENTROPY_LAPLACIAN_TRACE_GEOMETRY_DEGENERACY_NO_CLOSURE",
        "full_spectrum_principle_supplied": True,
        "full_spectrum_breaks_p2774_pair_degeneracy": witness["full_spectrum_breaks_p2774_pair_degeneracy"],
        "sourced_nadsoliton_spectral_law_exported": False,
        "global_graph_class_uniqueness_audited": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_pair_local_degeneracy_breaking_witness": facts["full_spectrum_breaks_p2774_pair_degeneracy"],
        "accepted_as_full_spectrum_nadsoliton_geometry_theorem": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_kernel_full_expression_theorem": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The full Laplacian spectrum is a valid finite discriminator for the P2774 pair, but current artifacts do not export a strict source law making that spectrum canonical, a global uniqueness theorem over the allowed graph class, or a variational coupling to K/L_total.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["full_laplacian_spectrum_witness"]
    lines = [
        "# P2775/S1725 full-Laplacian-spectrum pair discriminator",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Pair-local spectral result",
        f"- entropy_matches_alpha_geo={witness['entropy_matches_alpha_geo']}",
        f"- same_laplacian_trace={witness['same_laplacian_trace']}",
        f"- distinct_full_laplacian_spectrum_count={witness['distinct_full_laplacian_spectrum_count']}",
        f"- full_spectrum_breaks_p2774_pair_degeneracy={witness['full_spectrum_breaks_p2774_pair_degeneracy']}",
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
    p2774 = read_json(P2774)
    witness = spectral_witness()
    acceptance = acceptance_matrix(witness, p2774)
    payload = {
        "status": "P2775_FULL_LAPLACIAN_SPECTRUM_PAIR_DISCRIMINATOR_NO_CLOSURE",
        "input_hashes": {"P2774": sha(P2774)},
        "input_statuses": {"P2774": p2774.get("status")},
        "audited_question": "Does the full Laplacian spectrum break the P2774 entropy-plus-trace geometry degeneracy?",
        "full_laplacian_spectrum_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2775 only as a pair-local positive discriminator, not as canonical geometry closure.  The next honest move is exactly one sourced spectral principle: either an explicit nadsoliton spectral action/variational law with a target spectrum and bounded uniqueness test over a declared finite graph class, or a cospectral-degeneracy audit showing whether full spectrum still fails on a broader class.  Without that, preserve the P2697-P2775 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2775/S1725 full-Laplacian-spectrum pair discriminator", "## P2775/S1725 full-Laplacian-spectrum pair discriminator\n\n`P2775/S1725` tests one stronger metric principle after P2774: the full graph-Laplacian spectrum.  On the exact P2774 pair, `torus_4x4` and `circulant_pm1_pm2`, the full spectrum distinguishes the two geometries despite shared `H=ln(16)=4 ln 2`, degree sequence, edge count, and Laplacian trace.  This is only a pair-local discriminator: no strict nadsoliton spectral source law, graph-class uniqueness theorem, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2775/S1725 spectral discriminator Ltotal guard", "## P2775/S1725 spectral discriminator Ltotal guard\n\n`P2775/S1725` adds no variational source term.  The full Laplacian spectrum separates the P2774 pair, but current artifacts do not supply a sourced spectral action or nonproxy coupling to `K`/`L_total`; therefore the result remains a bounded pair-local geometry discriminator, not role-bearing `L_total` promotion.\n")
    append_once(AGENTS, "Current full-Laplacian-spectrum pair discriminator guardrail (P2775/S1725, 2026-06-15)", "## Current full-Laplacian-spectrum pair discriminator guardrail (P2775/S1725, 2026-06-15)\n\n- P2775 tests the P2774-recommended full graph-Laplacian-spectrum principle on the exact P2774 pair.\n- The full spectrum distinguishes `torus_4x4` from `circulant_pm1_pm2`, so it is a bounded pair-local degeneracy-breaking witness, but no strict nadsoliton spectral source law, finite graph-class uniqueness theorem, or `K`/`L_total` variational coupling is exported.\n- Do not promote this pair-local spectral discriminator to canonical nadsoliton geometry, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move must provide one sourced spectral action/variational law with bounded uniqueness, or run a broader cospectral-degeneracy audit.\n")
    return payload


if __name__ == "__main__":
    main()
