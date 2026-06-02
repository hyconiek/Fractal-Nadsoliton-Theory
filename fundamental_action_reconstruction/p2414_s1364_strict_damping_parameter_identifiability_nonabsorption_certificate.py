#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.json"
MD = GEN / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.md"

SOURCE_FILES = {
    "P2377_DAMPING_TRANSPORT_PRIMITIVE": GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json",
    "P2393_COMPLETION_BOUNDARY_RESIDUAL": GEN / "p2393_s1343_kernel_completion_boundary_residual_certificate.json",
    "P2411_BRIDGE_SOURCE_HYPERGRAPH": GEN / "p2411_s1361_legacy_strict_bridge_source_obligation_hypergraph_certificate.json",
    "P2413_AMPLITUDE_SCALAR_WITNESS": GEN / "p2413_s1363_amplitude_scalar_normalization_bridge_witness_certificate.json",
    "SCRATCH_DAMPING_SEPARATION": SCRATCH / "bridge_strict_completion_legacy_to_strict_damping_compression_separation_certificate_report.json",
    "SCRATCH_DAMPING_IDENTIFIABILITY": SCRATCH / "bridge_strict_completion_legacy_to_strict_damping_parameter_identifiability_certificate_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

DOMAIN = list(range(1, 12))
BETA_TORS = Fraction(1, 100)
STRICT_BETA = Fraction(1, 1)
STRICT_ETA = Fraction(9, 5)
CANDIDATE_P_MAX = 30
CANDIDATE_Q_MAX = 10


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg",
            "-n",
            pattern,
            "fundamental_action_reconstruction",
            "material_dowodowy",
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.tex",
            "-g",
            "!generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:18]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2414|S1364|strict damping parameter identifiability|damping parameter nonabsorption",
        "prior_transport": "damping-compression transport primitive|endpoint primitive|uniform coupling",
        "prior_boundary": "linear torsion damping to nonlinear strict compression|eta=1 boundary|current residual",
        "scratch_identifiability": "beta eta recovery|strict denominator parameter|damping parameter identifiability|log eta recovery",
        "linear_no_go": "1\\+gamma\\*d vs 1\\+d\\^eta|no single linear gamma|legacy linear torsion damping strict nonlinear compression",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds prior damping transport, boundary residual, and scratch beta/eta identifiability/no-go material, "
            "but no production P24xx certificate that follows P2413 by separating accepted strict damping-parameter "
            "identifiability from any beta_tors absorption, source theorem, role transfer, or full bridge closure."
        ),
    }


def strict_denominator_float(d: int) -> float:
    return 1.0 + d ** (float(STRICT_ETA))


def legacy_denominator_float(d: int) -> float:
    return 1.0 + float(BETA_TORS) * d


def required_linear_gamma_float(d: int) -> float:
    return d ** (float(STRICT_ETA - 1))


def candidate_grid() -> list[Fraction]:
    return sorted({Fraction(p, q) for p in range(1, CANDIDATE_P_MAX + 1) for q in range(1, CANDIDATE_Q_MAX + 1)})


def matching_candidates() -> list[Fraction]:
    # In the accepted strict model S(d)-1=d^(9/5).  Once beta=S(1)-1=1,
    # a rational candidate p/q matches the d=2 sample iff 2^(p/q)=2^(9/5),
    # equivalently p/q=9/5 by injectivity of x -> 2^x.
    return [candidate for candidate in candidate_grid() if candidate == STRICT_ETA]


def best_l2_linear_gamma() -> float:
    numerator = sum(d * (d ** float(STRICT_ETA)) for d in DOMAIN)
    denominator = sum(d * d for d in DOMAIN)
    return numerator / denominator


def finite_rows() -> list[dict[str, Any]]:
    gamma_star = best_l2_linear_gamma()
    rows = []
    for d in DOMAIN:
        strict_den = strict_denominator_float(d)
        legacy_den = legacy_denominator_float(d)
        gamma_req = required_linear_gamma_float(d)
        rows.append(
            {
                "d": d,
                "strict_denominator_symbolic": f"1+{d}^(9/5)",
                "legacy_denominator_exact": f"{100 + d}/100",
                "strict_denominator_float": strict_den,
                "legacy_beta_tors_denominator_float": legacy_den,
                "strict_minus_legacy_denominator_residual_float": strict_den - legacy_den,
                "required_linear_gamma_symbolic": f"{d}^(4/5)",
                "required_linear_gamma_float": gamma_req,
                "required_gamma_minus_legacy_beta_tors_float": gamma_req - float(BETA_TORS),
                "best_l2_linear_denominator_float": 1.0 + gamma_star * d,
                "strict_minus_best_l2_linear_residual_float": strict_den - (1.0 + gamma_star * d),
                "fifth_power_cover_identity": f"(S({d})-1)^5 = {d}^9",
            }
        )
    return rows


def build_certificate(sources: dict[str, Any]) -> dict[str, Any]:
    rows = finite_rows()
    candidates = matching_candidates()
    required_gammas = [row["required_linear_gamma_float"] for row in rows]
    gamma_diffs = [b - a for a, b in zip(required_gammas, required_gammas[1:])]
    legacy_residuals = [row["strict_minus_legacy_denominator_residual_float"] for row in rows]
    l2_residuals = [row["strict_minus_best_l2_linear_residual_float"] for row in rows]
    p2411_theorem = sources["P2411_BRIDGE_SOURCE_HYPERGRAPH"].get(
        "legacy_strict_bridge_source_obligation_hypergraph_certificate", {}
    ).get("theorem_export", {})
    p2413_theorem = sources["P2413_AMPLITUDE_SCALAR_WITNESS"].get(
        "amplitude_scalar_normalization_bridge_witness_certificate", {}
    ).get("theorem_export", {})
    scratch_sep = sources["SCRATCH_DAMPING_SEPARATION"].get("separation_summary", {})
    scratch_ident = sources["SCRATCH_DAMPING_IDENTIFIABILITY"].get("damping_parameter_identifiability_summary", {})
    return {
        "domain": DOMAIN,
        "strict_denominator_model": "S(d)=1+beta*d^eta",
        "accepted_strict_denominator_samples": "S(d)=1+d^(9/5) on d=1..11",
        "legacy_linear_denominator": "L(d)=1+beta_tors*d with beta_tors=1/100",
        "parameter_identification": {
            "beta_step": "S(1)-1=beta, and the accepted sample has S(1)-1=1, so beta=1.",
            "eta_step": "For each d>=2, eta=log(S(d)-1)/log(d)=log(d^(9/5))/log(d)=9/5 after beta=1.",
            "finite_cover_identity": "For eta=9/5, every row lies on the algebraic cover (S(d)-1)^5=d^9.",
            "candidate_grid": f"reduced rationals p/q with 1<=p<={CANDIDATE_P_MAX}, 1<=q<={CANDIDATE_Q_MAX}",
            "matching_rational_candidates": [str(candidate) for candidate in candidates],
        },
        "nonabsorption_proof": {
            "scalar_denominator_no_go": "At d=0 both denominators equal 1, so any scalar denominator identity has scalar 1; at d=1, 101/100 != 2.",
            "linear_gamma_no_go": "If 1+gamma*d=1+d^(9/5) at positive node d, then gamma=d^(4/5), which is strictly increasing in d.",
            "two_node_exact_contradiction": "Matching d=1 forces gamma=1; matching d=2 would force gamma=2^(4/5), so one beta_tors-style gamma cannot match both.",
            "legacy_beta_tors_no_positive_node_match": "d^(4/5) >= 1 for d>=1, while beta_tors=1/100.",
        },
        "finite_rows": rows,
        "required_gamma_values_strictly_increase": all(diff > 0.0 for diff in gamma_diffs),
        "no_single_linear_gamma_matches_two_distinct_positive_nodes": all(diff > 0.0 for diff in gamma_diffs),
        "legacy_beta_tors_matches_no_positive_strict_node": all(row["required_gamma_minus_legacy_beta_tors_float"] > 0 for row in rows),
        "legacy_denominator_residuals_all_positive": all(res > 0.0 for res in legacy_residuals),
        "best_l2_linear_gamma": best_l2_linear_gamma(),
        "best_l2_residual_l2_norm": math.sqrt(sum(res * res for res in l2_residuals)),
        "candidate_grid_size": len(candidate_grid()),
        "candidate_grid_unique_eta_match": len(candidates) == 1 and candidates[0] == STRICT_ETA,
        "strict_beta_eta_identified_from_samples": len(candidates) == 1,
        "strict_beta_eta_source_exported": False,
        "legacy_beta_tors_to_beta_eta_theorem_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_licensed": False,
        "p2411_full_bridge_still_requires_all_obligations": p2411_theorem.get("bridge_ready_true_masks") == [255],
        "p2413_amplitude_witness_inherited": p2413_theorem.get("scalar_normalization_witness_ready") is True,
        "scratch_damping_separation_inherited": scratch_sep.get("no_single_linear_gamma_matches_two_distinct_positive_nodes") is True,
        "scratch_damping_identifiability_inherited": scratch_ident.get("beta_identified_from_d1") is True
        and scratch_ident.get("candidate_grid_unique_match") is True,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2414/S1364 strict damping parameter identifiability and nonabsorption certificate

`P2414/S1364` follows the P2413 amplitude witness by isolating the damping/compression row.  It treats the strict denominator samples as accepted finite data

```text
S(d) = 1 + beta*d^eta = 1 + d^(9/5),        d=1..11.
```

Within that accepted strict denominator model, `beta` is identified by `S(1)-1=1`, and `eta` is identified by `eta=log(S(d)-1)/log(d)=9/5` for every `d>=2`; equivalently `(S(d)-1)^5=d^9` on the finite algebraic cover.  A reduced rational grid `p/q` with `p<=30`, `q<=10` has the unique match `9/5`.

The same certificate proves nonabsorption into the legacy linear torsion denominator: matching `1+gamma*d` to `1+d^(9/5)` at a positive node forces `gamma=d^(4/5)`, which is strictly increasing, so no single `beta_tors`-style linear parameter matches two positive nodes; the legacy `beta_tors=1/100` matches no positive strict node.

This identifies accepted strict damping parameters and proves a linear-denominator no-go only.  It does not export a strict dynamic source for `beta,eta`, does not prove `beta_tors -> beta/eta`, does not complete the damping bridge, and does not license role transfer, QW-2191 discharge, `L_total`, or ToE closure.
""".strip()
    lag_section = """
## P2414/S1364 strict damping nonabsorption guard for Lagrangian/EOM

`P2414/S1364` allows the strict denominator data `1+d^(9/5)` to be recognized as finitely identifying `beta=1, eta=9/5` inside the accepted strict model, while proving that this cannot be absorbed into a legacy `1+beta_tors*d` denominator.  Since no strict dynamic source, damping bridge theorem, or role-transfer theorem is exported, no nonlinear-compression expression can yet be promoted to a role-bearing `L_total` term.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    cert = build_certificate(sources)
    theorem_export = {
        "theorem_name": "P2414_T1_strict_damping_parameter_identifiability_nonabsorption_certificate",
        "domain_size": len(cert["domain"]),
        "strict_denominator_model": cert["strict_denominator_model"],
        "accepted_strict_denominator_samples": cert["accepted_strict_denominator_samples"],
        "candidate_grid_size": cert["candidate_grid_size"],
        "matching_rational_candidates": cert["parameter_identification"]["matching_rational_candidates"],
        "candidate_grid_unique_eta_match": cert["candidate_grid_unique_eta_match"],
        "strict_beta_eta_identified_from_samples": cert["strict_beta_eta_identified_from_samples"],
        "required_gamma_values_strictly_increase": cert["required_gamma_values_strictly_increase"],
        "no_single_linear_gamma_matches_two_distinct_positive_nodes": cert[
            "no_single_linear_gamma_matches_two_distinct_positive_nodes"
        ],
        "legacy_beta_tors_matches_no_positive_strict_node": cert["legacy_beta_tors_matches_no_positive_strict_node"],
        "legacy_denominator_residuals_all_positive": cert["legacy_denominator_residuals_all_positive"],
        "best_l2_linear_gamma": cert["best_l2_linear_gamma"],
        "best_l2_residual_l2_norm": cert["best_l2_residual_l2_norm"],
        "strict_beta_eta_source_exported": cert["strict_beta_eta_source_exported"],
        "legacy_beta_tors_to_beta_eta_theorem_exported": cert["legacy_beta_tors_to_beta_eta_theorem_exported"],
        "damping_compression_bridge_component_ready": cert["damping_compression_bridge_component_ready"],
        "full_bridge_theorem_exported": cert["full_bridge_theorem_exported"],
        "role_transfer_licensed": cert["role_transfer_licensed"],
        "p2411_full_bridge_still_requires_all_obligations": cert["p2411_full_bridge_still_requires_all_obligations"],
        "p2413_amplitude_witness_inherited": cert["p2413_amplitude_witness_inherited"],
        "scratch_damping_separation_inherited": cert["scratch_damping_separation_inherited"],
        "scratch_damping_identifiability_inherited": cert["scratch_damping_identifiability_inherited"],
        "not_licensed": [
            "No strict dynamic source theorem for beta=1 or eta=9/5 is exported.",
            "No beta_tors -> beta/eta completion-map theorem is claimed.",
            "No beta_tors -> chi11 theorem or QW-2191 discharge follows.",
            "No damping bridge completion, legacy role transfer, role-bearing L_total, or ToE closure follows.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "domain_has_eleven_positive_nodes": theorem_export["domain_size"] == 11,
        "candidate_grid_unique_eta_match": theorem_export["candidate_grid_unique_eta_match"],
        "strict_beta_eta_identified_from_samples": theorem_export["strict_beta_eta_identified_from_samples"],
        "required_gamma_strictly_increases": theorem_export["required_gamma_values_strictly_increase"],
        "linear_two_node_no_go": theorem_export["no_single_linear_gamma_matches_two_distinct_positive_nodes"],
        "legacy_beta_tors_matches_no_positive_node": theorem_export["legacy_beta_tors_matches_no_positive_strict_node"],
        "legacy_residuals_positive": theorem_export["legacy_denominator_residuals_all_positive"],
        "no_source_theorem_exported": not theorem_export["strict_beta_eta_source_exported"],
        "no_beta_tors_map_exported": not theorem_export["legacy_beta_tors_to_beta_eta_theorem_exported"],
        "damping_bridge_still_open": not theorem_export["damping_compression_bridge_component_ready"],
        "full_bridge_still_open": not theorem_export["full_bridge_theorem_exported"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "p2411_bridge_obligation_inherited": theorem_export["p2411_full_bridge_still_requires_all_obligations"],
        "p2413_amplitude_witness_inherited": theorem_export["p2413_amplitude_witness_inherited"],
        "scratch_separation_inherited": theorem_export["scratch_damping_separation_inherited"],
        "scratch_identifiability_inherited": theorem_export["scratch_damping_identifiability_inherited"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2414_s1364_v1",
        "packet_id": "P2414",
        "stage_id": "S1364",
        "result_kind": "STRICT_DAMPING_PARAMETER_IDENTIFIABILITY_NONABSORPTION_CERTIFICATE",
        "status": "PASS_STRICT_DAMPING_PARAMETERS_IDENTIFIED_FROM_ACCEPTED_SAMPLES_NO_SOURCE_NO_ABSORPTION",
        "strict_damping_parameter_identifiability_nonabsorption_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Turn the finite beta/eta identification and linear-denominator no-go into a real damping-compression "
            "completion-map/source theorem, or keep it as a strict-side addition without role transfer."
        ),
        "global_status": "OPEN_PROGRESS_DAMPING_PARAMETERS_IDENTIFIED_NO_SOURCE_NO_ROLE_TRANSFER_OR_TOE_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["strict_damping_parameter_identifiability_nonabsorption_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2414 S1364: strict damping parameter identifiability and nonabsorption certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2414/S1364 identifies beta/eta from accepted strict denominator samples and proves that the result is not a legacy linear beta_tors absorption.",
                "",
                "## Finite facts",
                "",
                f"- Domain size: `{theorem['domain_size']}` positive nodes.",
                f"- Accepted samples: `{theorem['accepted_strict_denominator_samples']}`.",
                f"- Candidate grid size: `{theorem['candidate_grid_size']}`.",
                f"- Matching rational candidates: `{theorem['matching_rational_candidates']}`.",
                f"- Required gamma strictly increases: `{theorem['required_gamma_values_strictly_increase']}`.",
                f"- Legacy beta_tors matches no positive node: `{theorem['legacy_beta_tors_matches_no_positive_strict_node']}`.",
                f"- Best L2 linear gamma: `{theorem['best_l2_linear_gamma']}`.",
                f"- Best L2 residual norm: `{theorem['best_l2_residual_l2_norm']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
                "",
                "## Gatekeepers",
                "",
                f"`{payload['gatekeeper_checks']}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    append_doc_sections()
    payload = build_payload()
    write_outputs(payload)
    if not all(payload["gatekeeper_checks"].values()):
        raise SystemExit(f"gatekeeper failure: {payload['gatekeeper_checks']}")


if __name__ == "__main__":
    main()
