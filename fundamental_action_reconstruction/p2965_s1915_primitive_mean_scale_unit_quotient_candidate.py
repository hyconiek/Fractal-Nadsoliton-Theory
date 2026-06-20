#!/usr/bin/env python3
"""P2965/S1915: primitive-mean scale/unit quotient candidate gate.

P2964 exposed one missing strict object: a canonical scale/unit quotient with a
coefficient source law for the aggregate-reception coupling.  This step attacks
that exact object without replaying beta normalization or scalar Euler insertion.
It constructs a ray quotient for positive rational multiples of the P2938/P2961
aggregate vector and computes the primitive integer representative plus its mean
coefficient 9/5.

The construction is a proof-grade developmental candidate: it removes arbitrary
positive rescaling inside the aggregate ray and recovers eta=9/5 from the
primitive representative.  It is not strict closure, because current artifacts do
not export this primitive-mean quotient as a physical unit law or as a sourced
coefficient in a nonproxy action density.
"""
from __future__ import annotations

import hashlib, json, math
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2964_s1914_unit_bearing_nonproxy_coupling_reception_no_go import OUT as P2964

OUT = GEN / "p2965_s1915_primitive_mean_scale_unit_quotient_candidate.json"
MD = GEN / "p2965_s1915_primitive_mean_scale_unit_quotient_candidate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

BASE_VECTOR = [1, 2, 2, 2, 2]


def primitive_integer_representative(values: list[Fraction]) -> list[int]:
    den_lcm = 1
    for x in values:
        den_lcm = math.lcm(den_lcm, x.denominator)
    ints = [int(x * den_lcm) for x in values]
    g = 0
    for x in ints:
        g = math.gcd(g, abs(x))
    if g == 0:
        raise ValueError("zero vector has no positive primitive representative")
    ints = [x // g for x in ints]
    if min(ints) < 0:
        ints = [-x for x in ints]
    return ints


def ray_samples() -> list[dict[str, Any]]:
    scalars = [Fraction(1, 4), Fraction(1, 3), Fraction(1, 2), Fraction(1, 1), Fraction(3, 2), Fraction(2, 1), Fraction(5, 1)]
    rows = []
    for q in scalars:
        scaled = [q * x for x in BASE_VECTOR]
        primitive = primitive_integer_representative(scaled)
        rows.append({
            "scale": f"{q.numerator}/{q.denominator}",
            "scaled_vector": [f"{x.numerator}/{x.denominator}" for x in scaled],
            "primitive_representative": primitive,
            "primitive_sum": sum(primitive),
            "primitive_mean": f"{sum(primitive)}/{len(primitive)}",
            "eta": sum(primitive) / len(primitive),
            "ray_quotient_invariant": primitive == BASE_VECTOR,
        })
    return rows


def candidate_quotient_rows() -> list[dict[str, Any]]:
    rows = [
        {
            "candidate": "primitive_ray_mean_quotient",
            "definition": "positive rational ray -> primitive integer representative -> arithmetic mean",
            "scale_orbit_removed": True,
            "returns_eta_9_over_5": True,
            "target_sum_cut": False,
            "beta_normalization_replay": False,
            "physical_unit_exported": False,
            "strict_coefficient_source_law": False,
            "nonproxy_coupling_installed": False,
        },
        {
            "candidate": "l1_normalized_probability_mean",
            "definition": "V -> V/sum(V)",
            "scale_orbit_removed": True,
            "returns_eta_9_over_5": False,
            "target_sum_cut": False,
            "beta_normalization_replay": False,
            "physical_unit_exported": False,
            "strict_coefficient_source_law": False,
            "nonproxy_coupling_installed": False,
        },
        {
            "candidate": "forced_sum_9_section",
            "definition": "choose scale so sum(V)=9",
            "scale_orbit_removed": True,
            "returns_eta_9_over_5": True,
            "target_sum_cut": True,
            "beta_normalization_replay": False,
            "physical_unit_exported": False,
            "strict_coefficient_source_law": False,
            "nonproxy_coupling_installed": False,
        },
        {
            "candidate": "beta_length_normalization_section",
            "definition": "choose length ell so beta ell^eta=1",
            "scale_orbit_removed": True,
            "returns_eta_9_over_5": False,
            "target_sum_cut": False,
            "beta_normalization_replay": True,
            "physical_unit_exported": False,
            "strict_coefficient_source_law": False,
            "nonproxy_coupling_installed": False,
        },
    ]
    for row in rows:
        row["developmental_candidate_accepted"] = row["scale_orbit_removed"] and row["returns_eta_9_over_5"] and not row["target_sum_cut"] and not row["beta_normalization_replay"]
        row["strict_coupling_source_accepted"] = row["developmental_candidate_accepted"] and row["physical_unit_exported"] and row["strict_coefficient_source_law"] and row["nonproxy_coupling_installed"]
    return rows


def obligation_rows(rows: list[dict[str, Any]], samples: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {"obligation": "positive_ray_quotient_constructed", "satisfied": True, "evidence": "primitive representative is computed by lcm denominators and gcd reduction"},
        {"obligation": "finite_scale_samples_invariant", "satisfied": all(r["ray_quotient_invariant"] for r in samples), "evidence": f"{len(samples)} positive rational scale samples collapse to [1,2,2,2,2]"},
        {"obligation": "coefficient_9_over_5_computed_without_sum_cut", "satisfied": any(r["candidate"] == "primitive_ray_mean_quotient" and r["developmental_candidate_accepted"] for r in rows), "evidence": "mean of primitive representative is 9/5"},
        {"obligation": "physical_unit_law_exported", "satisfied": False, "evidence": "primitive ray quotient is dimensionless arithmetic unless a strict physical unit law is added"},
        {"obligation": "strict_coefficient_source_law_exported", "satisfied": False, "evidence": "no current artifact promotes the primitive mean to a sourced action coefficient"},
        {"obligation": "nonproxy_coupling_installed", "satisfied": False, "evidence": "P2964 coupling installation remains missing"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["ray_quotient", "scale_sample_invariance", "eta_9_over_5_without_sum_cut", "no_beta_replay", "physical_unit_law", "strict_coefficient_source", "nonproxy_coupling"]
    full = (1 << len(names)) - 1
    dev_mask = sum(1 << i for i in range(4))
    return [{"mask": m, "present": {n: bool(m & (1 << i)) for i, n in enumerate(names)}, "accepts_developmental_quotient_candidate": (m & dev_mask) == dev_mask, "accepts_strict_coupling_source": m == full} for m in range(1 << len(names))]


def build_payload(p2964: dict[str, Any]) -> dict[str, Any]:
    samples = ray_samples()
    rows = candidate_quotient_rows()
    obligations = obligation_rows(rows, samples)
    matrix = acceptance_matrix()
    return {
        "status": "P2965_PRIMITIVE_MEAN_SCALE_UNIT_QUOTIENT_CANDIDATE_NO_STRICT_EXPORT",
        "input_hashes": {"P2964": hashlib.sha256(P2964.read_bytes()).hexdigest() if P2964.exists() else None},
        "unbounded_lemma": "For any positive rational q, primitive_integer_representative(q*[1,2,2,2,2])=[1,2,2,2,2]; hence the primitive mean is always 9/5 on this ray.",
        "constructed_theoretical_objects": {"candidate_object": "PrimitiveMean_ScaleUnitQuotient_CandidateGate", "ray_sample_rows": samples, "candidate_quotient_rows": rows, "proof_obligation_rows": obligations, "finite_acceptance_matrix": matrix},
        "quotient_certificate": {"sample_count": len(samples), "all_samples_invariant": all(r["ray_quotient_invariant"] for r in samples), "developmental_candidates": [r["candidate"] for r in rows if r["developmental_candidate_accepted"]], "strict_coupling_sources": [r["candidate"] for r in rows if r["strict_coupling_source_accepted"]], "acceptance_matrix_rows": len(matrix), "developmental_rows": sum(1 for r in matrix if r["accepts_developmental_quotient_candidate"]), "strict_rows": sum(1 for r in matrix if r["accepts_strict_coupling_source"])},
        "decision": {
            "positive_progress": "P2965 constructs the missing quotient object as a developmental theorem candidate: positive rescalings of the aggregate ray collapse to the primitive representative [1,2,2,2,2], whose mean is 9/5 without a target_sum cut.",
            "breakthrough": "No strict L_total coupling follows: the quotient is still dimensionless arithmetic until a physical unit law, strict coefficient source law, and nonproxy coupling installation are exported.",
            "negative_export_flags": {k: False for k in ["strict_scale_unit_quotient_exported", "strict_coefficient_source_law_exported", "strict_unit_bearing_nonproxy_coupling_exported", "strict_ratio_package_source_exported", "damping_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay primitive-ray quotient arithmetic, target_sum=9 sections, beta length normalization, scalar Euler insertion, K/C exchange, or typed scalar mediators.  The next proof-grade move must either prove a strict physical unit law that makes the P2965 primitive-mean quotient unit-bearing, prove a strict coefficient source law coupling that quotient into the P2964 action-density schema, or pivot to a genuinely new typed structural object outside the current ratio-package lane while preserving the P2929-P2965 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["quotient_certificate"]
    lines = ["# P2965/S1915 primitive-mean scale/unit quotient candidate", "", f"Status: `{payload['status']}`", "", "## Quotient certificate", f"- sample count: `{cert['sample_count']}`", f"- all scale samples invariant: `{cert['all_samples_invariant']}`", f"- developmental candidates: `{cert['developmental_candidates']}`", f"- strict coupling sources: `{cert['strict_coupling_sources']}`", f"- acceptance matrix rows/developmental/strict: `{cert['acceptance_matrix_rows']}/{cert['developmental_rows']}/{cert['strict_rows']}`", "", "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2964))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2965/S1915 primitive-mean scale/unit quotient candidate", "## P2965/S1915 primitive-mean scale/unit quotient candidate\n\n`P2965/S1915` attacks exactly one P2964 missing object: a canonical scale/unit quotient with coefficient source law.  It constructs a positive rational ray quotient for the P2938/P2961 aggregate: every positive rational rescaling of `[1,2,2,2,2]` reduces by lcm/gcd to the same primitive representative, whose mean is `9/5` without a target_sum cut.  This is proof-grade developmental progress, but it remains dimensionless arithmetic until a strict physical unit law, coefficient source law, and nonproxy coupling installation are exported; no strict ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2965/S1915 primitive-mean quotient `L_total` guard", "## P2965/S1915 primitive-mean quotient `L_total` guard\n\n`P2965/S1915` supplies a clean primitive-ray quotient candidate and computes the coefficient `9/5` from the primitive aggregate mean.  Because no strict physical unit law or coefficient/source map installs that coefficient into the P2964 action-density schema, no sourced damping coefficient enters `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE from P2965.\n")
    append_once(AGENTS, "Current primitive-mean scale/unit quotient candidate guardrail (P2965/S1915, 2026-06-20)", "## Current primitive-mean scale/unit quotient candidate guardrail (P2965/S1915, 2026-06-20)\n\n- P2965 constructs the P2964-requested quotient object as a developmental theorem candidate: positive rational rescalings of the P2938/P2961 aggregate ray reduce to primitive representative `[1,2,2,2,2]`, whose mean is `9/5`, without a target_sum cut.\n- This does not export strict closure: the quotient is still dimensionless arithmetic until a strict physical unit law, strict coefficient source law, and nonproxy coupling installation are proved.\n- Do not promote primitive-ray quotient arithmetic, target_sum sections, beta normalization, scalar Euler insertion, K/C exchange, or typed scalar mediators to strict ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must prove exactly one missing strict object: a physical unit law for the P2965 quotient, a coefficient source law installing it into the P2964 coupling schema, or a genuinely new typed structural object outside the current ratio-package lane; otherwise preserve the P2929-P2965 boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
