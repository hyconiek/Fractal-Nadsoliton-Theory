#!/usr/bin/env python3
"""Recompute the four fixed native interval stresses; preserve AD unchanged."""
import argparse
import ast
import hashlib
import json
import re
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = [
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "fin_global_activity_defect_research_20261020/CONTRACT.md",
    "fin_global_activity_defect_research_20261020/GLOBAL_APPROXIMATE_SCALING.md",
    "fin_global_activity_defect_research_20261020/LOCALIZED_LOG_FLOW_THEOREM.md",
    "fin_global_activity_defect_research_20261020/ACTUAL_DRIFT_NOISE_PAYMENTS.md",
    "fin_global_activity_defect_research_20261020/ORIGINAL_SOURCE_PATH_THEOREM.md",
    "fin_global_activity_defect_research_20261020/NEXT_PROOF_AUDIT.md",
    "fin_global_activity_defect_research_20261020/MANIFEST.sha256",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
    "fin_intrinsic_scaling_research_20261020/GLOBAL_SCALING_THEOREM.md",
    "fin_intrinsic_scaling_research_20261020/ACTUAL_LOCAL_WINDOW.md",
    "fin_intrinsic_scaling_research_20261020/DISCRETE_CRITICAL_GATE.md",
    "fin_sqrt_critical_process_research_20261020/STATIONARY_PATH_THEOREM.md",
    "fin_hyp_history_entropy_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_hyp_history_entropy_review_20261020/PROOF_AUDIT_AND_DUAL.md",
    "fin_reference_transport_research_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
    "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py",
]
PROTECTED = [
    "fin_global_activity_defect_research_20261020",
    "fin_intrinsic_scaling_research_20261020",
    "fin_sqrt_critical_process_research_20261020",
    "fin_sqrt_local_regime_research_20261020",
    "fin_spectator_limit_research_20261020",
    "fin_spectator_composition_research_20261020",
    "fin_reference_transport_research_20261020",
    "fin_hyp_history_entropy_20261020",
    "fin_hyp_history_entropy_review_20261020",
    "fin_critical_metastability_review_20261015",
    "fin_critical_metastability_campaign_20261014",
    "fin_uniform_critical_campaign_20261014",
    "fin_uniform_critical_review_20261014",
    "fin_intrinsic_observation_review_20261015",
    "fin_stationary_kinetic_rigidity_review_20261016",
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder):
    seen = set()
    for line in (folder/"MANIFEST.sha256").read_text().splitlines():
        h, name = line.split("  ", 1)
        assert name not in seen and not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder/name) == h, (folder.name, name)
        seen.add(name)
    actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
    actual |= {n for n in seen if "__pycache__" in Path(n).parts or n.endswith(".pyc")}
    assert actual == seen | {"MANIFEST.sha256"}, (folder.name, actual-seen, seen-actual)
    return {"entries": len(seen), "sha256": sha(folder/"MANIFEST.sha256"), "coverage": "complete_noncache"}


def check_links():
    for p in HERE.glob("*.md"):
        prose, fence = [], False
        for line in p.read_text().splitlines():
            assert line == line.rstrip(), p.name
            if line.lstrip().startswith("```"):
                fence = not fence
            elif not fence and not line.startswith(("    ", "\t")):
                prose.append(line)
        text = re.sub(r"`[^`\n]*`", "", "\n".join(prose))
        for name in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", text):
            name = name.strip("<>").split("#")[0]
            if name and not name.startswith(("http:", "https:")):
                assert (p.parent/name).exists(), (p.name, name)


def main():
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    if args.check:
        manifest(HERE)
    before = {d: manifest(ROOT/d) for d in PROTECTED}
    inputs = {n: sha(ROOT/n) for n in INPUTS}
    from interval_stress import I, J, POPULATIONS, native_context, run_case
    checks = []
    def ck(name, ok):
        assert name not in checks and ok, name
        checks.append(name)
    ck("dyadic_addition_outward", (I(F(1, 3))+I(F(2, 7))).contains(F(13, 21)))
    ck("dyadic_multiplication_outward", (I(F(2, 7))*I(3)).contains(F(6, 7)))
    ck("dyadic_reciprocal_outward", I(F(1, 7)).inv().contains(7))
    ck("integer_root_outward", (I(3).sqrt()*I(3).sqrt()).contains(3))
    ck("positive_exp_series_zero_exact", I(0).exp().lo == I(0).exp().hi == I(1).lo)
    x = J(2, 1, 0)
    ck("jet_product_chain_rule", (x*x).v.contains(4) and (x*x).d.contains(4) and (x*x).s.contains(1))
    inv = x.inv()
    ck("jet_reciprocal_chain_rule", inv.v.contains(F(1, 2)) and inv.d.contains(F(-1, 4)) and inv.s.contains(F(1, 8)))
    ex = J(0, 1, 0).exp()
    ck("jet_exp_chain_rule", ex.v.contains(1) and ex.d.contains(1) and ex.s.contains(F(1, 2)))
    sq = J(4, 1, 0).sqrt()
    ck("jet_sqrt_chain_rule", sq.v.contains(2) and sq.d.contains(F(1, 4)) and sq.s.contains(F(-1, 64)))
    h = F(1, 7)
    intervals = [(F(0), h), (h-F(1, 1000), h+F(1, 1000)), (F(2, 13), F(9, 11))]
    ck("grid_bracket_count_pays_extra_h", all((int(t//h)-int(s//h))*h <= t-s+h for s, t in intervals))
    s, t = intervals[1]
    ck("continuous_density_without_grid_h_countercontrol", (int(t//h)-int(s//h))*h > t-s)
    context = native_context()
    cases = [run_case(name, context) for name in POPULATIONS]
    for row in cases:
        box = row["uniform_half_second_derivative_box"]["exact_dyadic"]
        scale = 1 << box["bits"]
        row["public_integer_remainder_bound"] = (max(abs(box["lo"]), abs(box["hi"]))+scale-1)//scale
    certified = {row["point"]: row["uniform_trial_positivity_on_epsilon_interval_certified"] for row in cases}
    finite_positive = sum(row["trial_positive_at_all_legal_neighbors_certified"] and row["analysis_self_probability_positive_certified"] for case in cases for row in case["direct_points"])
    total_checks = len(checks)+sum(len(row["checks"]) for row in cases)
    certificate = {"date": "2026-10-20", "status": "FIXED_POPULATION_TAYLOR_ENCLOSURES_WITH_POSITIVITY_LIMIT_RETAINED", "matrix": "actual A7 via admitted exact spectral boxes and q12 Fourier projectors", "spectral_sha256": context["spectral_sha256"], "gain_box": context["g"].box(), "lambda": 0, "bits": 80, "case_results": cases, "uniform_in_N_scope": "each frozen p and applicable integer N>=264, conditional on parent analytic Phi(0)=1/Phi'(0)=0 identities", "global_simplex_certificate": False, "uniform_positivity_certified_by_case": certified, "finite_direct_positive_count": finite_positive, "finite_direct_evaluations": 8, "independent_audit": False}
    assert before == {d: manifest(ROOT/d) for d in PROTECTED}
    assert inputs == {n: sha(ROOT/n) for n in INPUTS}
    result = {"date": "2026-10-20", "status": "BOUNDED_REMAINDER_STRESS_WITH_VERTEX_POSITIVITY_INCONCLUSIVE", "bounded_checks": total_checks, "arithmetic_controls": checks, "fixed_population_checks": {row["point"]: row["checks"] for row in cases}, "public_integer_remainder_bounds": {row["point"]: row["public_integer_remainder_bound"] for row in cases}, "uniform_positivity_certified_by_case": certified, "direct_positive_count": finite_positive, "source_inputs": inputs, "interval_engine_sha256": sha(HERE/"interval_stress.py"), "protected_packages": before, "protected_unchanged": True, "AGENTS_full_sha256": inputs["AGENTS.md"], "original_AD_modified": False, "same_author": True, "user_declined_delegation": True, "independent_AD_audit_executed": False, "new_scientific_campaign": False, "global_simplex_theorem_certified": False, "strict_F01": "OPEN_SOURCE", "broad_goal": "ACTIVE_NOT_COMPLETE"}
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    if args.record:
        (HERE/"CERTIFICATE.json").write_text(json.dumps(certificate, indent=2)+"\n")
        (HERE/"RESULTS.json").write_text(json.dumps(result, indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text("".join(f"{h}  {n}\n" for n, h in inputs.items()))
        check_links()
        files = sorted(p for p in HERE.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        assert json.loads((HERE/"CERTIFICATE.json").read_text()) == certificate
        assert json.loads((HERE/"RESULTS.json").read_text()) == result
        assert (HERE/"INPUTS.sha256").read_text() == "".join(f"{h}  {n}\n" for n, h in inputs.items())
        check_links()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "bounded_checks": total_checks, "inputs": len(inputs), "positive_whole_epsilon_boxes": sum(certified.values()), "whole_epsilon_boxes": 4, "positive_direct_controls": finite_positive, "direct_controls": 8, "protected_packages": len(before), "protected_entries": sum(v["entries"] for v in before.values()), "AD_or_AGENTS_modified": False, "independent_audit": False}, indent=2))


if __name__ == "__main__":
    main()
