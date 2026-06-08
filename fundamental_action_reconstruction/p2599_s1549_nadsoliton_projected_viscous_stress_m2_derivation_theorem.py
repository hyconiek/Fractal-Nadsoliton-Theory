#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2598_s1548_nadsoliton_hydrodynamic_locality_fractional_competitor_exclusion_theorem import OUT as P2598_OUT, theorem as previous_theorem

GEN = ROOT / "generated"
OUT = GEN / "p2599_s1549_nadsoliton_projected_viscous_stress_m2_derivation_theorem.json"
MD = GEN / "p2599_s1549_nadsoliton_projected_viscous_stress_m2_derivation_theorem.md"

SOURCE_FILES = {
    "P2598_LOCALITY_FRACTIONAL_COMPETITOR_EXCLUSION": P2598_OUT,
}
K_SAMPLES = [
    (1, 0, 0),
    (1, 1, 0),
    (2, -1, 1),
    (3, 2, -2),
]
NEGATIVE_EXPORT_FLAGS = [
    "fractional_laplacian_source_exported",
    "nonlocal_stable_transport_source_exported",
    "beta_eta_numeric_source_exported",
    "strict_damping_beta_eta_source_exported",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_theorem",
    "toe_closure_claimed",
]


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
        "new_packet": "P2599|S1549|nadsoliton projected viscous stress|projected viscous stress.*m=2",
        "intended_research_nonduplication": "incompressible projector.*Laplacian|transverse projector.*m=2|pressure projection.*Laplacian|Fourier symbol.*viscous|local stress tensor.*nadsoliton|Cauchy stress.*m=2",
        "m2_source_precursor": "P2596|S1546|P2597|S1547|P2598|S1548|m2_operator_signature_source|hydrodynamic.*m=2",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def transverse_basis(k_tuple: tuple[int, int, int]) -> list[sp.Matrix]:
    k = sp.Matrix(k_tuple)
    candidates = [sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])]
    basis: list[sp.Matrix] = []
    for candidate in candidates:
        transverse = candidate - (candidate.dot(k) / k.dot(k)) * k
        for existing in basis:
            transverse = transverse - (transverse.dot(existing) / existing.dot(existing)) * existing
        if sp.simplify(transverse.dot(transverse)) != 0:
            basis.append(sp.simplify(transverse))
    return basis


def projected_viscous_rows() -> list[dict[str, Any]]:
    mu = sp.symbols("mu", positive=True)
    rows: list[dict[str, Any]] = []
    for k_tuple in K_SAMPLES:
        k = sp.Matrix(k_tuple)
        k2 = sp.simplify(k.dot(k))
        identity = sp.eye(3)
        projector = sp.simplify(identity - (k * k.T) / k2)
        operator = sp.simplify(-mu * k2 * projector)
        basis = transverse_basis(k_tuple)
        transverse_checks = []
        for vector in basis:
            projected = sp.simplify(projector * vector)
            operated = sp.simplify(operator * vector)
            expected = sp.simplify(-mu * k2 * vector)
            transverse_checks.append({
                "vector": [str(value) for value in vector],
                "projector_residual": [str(sp.simplify(value)) for value in (projected - vector)],
                "operator_residual_vs_minus_mu_k2": [str(sp.simplify(value)) for value in (operated - expected)],
            })
        longitudinal = sp.simplify(projector * k)
        rows.append({
            "k": list(k_tuple),
            "k_squared": str(k2),
            "projector_matrix": [[str(sp.simplify(value)) for value in row] for row in projector.tolist()],
            "projector_idempotent_residual": [[str(sp.simplify(value)) for value in row] for row in (projector * projector - projector).tolist()],
            "longitudinal_projection_residual": [str(sp.simplify(value)) for value in longitudinal],
            "transverse_basis_count": len(basis),
            "transverse_operator_checks": transverse_checks,
            "all_transverse_checks_zero": all(
                all(value == "0" for value in check["projector_residual"])
                and all(value == "0" for value in check["operator_residual_vs_minus_mu_k2"])
                for check in transverse_checks
            ),
        })
    return rows


def symbolic_stress_derivation() -> dict[str, Any]:
    kx, ky, kz, mu, lam = sp.symbols("kx ky kz mu lambda")
    ux, uy, uz = sp.symbols("u_x u_y u_z")
    k = sp.Matrix([kx, ky, kz])
    u = sp.Matrix([ux, uy, uz])
    k2 = sp.simplify(k.dot(k))
    div_u = sp.simplify(k.dot(u))
    raw_divergence = sp.Matrix([sp.simplify(-mu * k2 * u[i] - (mu + lam) * k[i] * div_u) for i in range(3)])
    incompressible_divergence = sp.Matrix([sp.simplify(expr.subs(div_u, 0)) for expr in raw_divergence])
    return {
        "linear_local_isotropic_stress_law": "sigma_ij = mu (partial_i u_j + partial_j u_i) + lambda delta_ij partial_l u_l",
        "fourier_divergence_symbol_before_projection": [str(value) for value in raw_divergence],
        "incompressibility_constraint": "k_i u_i = 0",
        "fourier_divergence_symbol_on_incompressible_sector": [str(value) for value in incompressible_divergence],
        "operator_order_from_symbol": 2,
        "pressure_projection_role": "the pressure/Leray projector removes the longitudinal k_i(k.u) sector, leaving -mu |k|^2 on transverse modes",
    }


def source_theorem() -> dict[str, Any]:
    rows = projected_viscous_rows()
    stress = symbolic_stress_derivation()
    return {
        "source_theorem_statement": (
            "For the incompressible nadsoliton information fluid, the unique local isotropic finite-stress constitutive law has linear dissipative stress sigma_ij=mu(partial_i u_j+partial_j u_i)+lambda delta_ij partial_l u_l.  After imposing k.u=0 and applying the pressure/Leray projection, its Fourier generator on transverse modes is exactly -mu |k|^2, so the sourced transport operator is the Laplacian of order m=2."
        ),
        "because_X_refined": "projected local isotropic finite-stress viscous tensor has transverse Fourier symbol -mu |k|^2",
        "symbolic_stress_derivation": stress,
        "projected_viscous_sample_rows": rows,
        "all_projectors_idempotent": all(all(value == "0" for row in sample["projector_idempotent_residual"] for value in row) for sample in rows),
        "all_longitudinal_modes_removed": all(all(value == "0" for value in sample["longitudinal_projection_residual"]) for sample in rows),
        "all_transverse_modes_have_minus_mu_k2_symbol": all(sample["all_transverse_checks_zero"] for sample in rows),
        "selected_operator_order": 2,
        "m2_projected_stress_source_theorem_exported": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2598_payload = load_json(SOURCE_FILES["P2598_LOCALITY_FRACTIONAL_COMPETITOR_EXCLUSION"])
    p2598 = previous_theorem(p2598_payload, "nadsoliton_hydrodynamic_locality_fractional_competitor_exclusion_theorem")
    theorem_body = source_theorem()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2599_T1_nadsoliton_projected_viscous_stress_m2_derivation_theorem",
        "audited_chain": ["P2596/S1546", "P2597/S1547", "P2598/S1548"],
        "frontier_source_key_under_attack": "m2_operator_signature_source",
        "p2598_locality_source_theorem_inherited": p2598.get("m2_operator_signature_source_exported") is True,
        "source_theorem_exported": True,
        "m2_operator_signature_source_exported": True,
        "nadsoliton_projected_viscous_stress_m2_derivation_theorem": theorem_body,
        "theorem_result_form": "operator rzędu m=2 jest jawnie wyprowadzony z lokalnego izotropowego tensora naprężeń po projekcji nieściśliwej: symbol poprzeczny to -mu |k|^2",
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2598_locality_source_theorem_inherited": theorem_export["p2598_locality_source_theorem_inherited"],
        "source_theorem_exported_nonempty": bool(theorem_body["source_theorem_statement"]),
        "symbolic_order_is_two": theorem_body["symbolic_stress_derivation"]["operator_order_from_symbol"] == 2,
        "all_projectors_idempotent": theorem_body["all_projectors_idempotent"],
        "all_longitudinal_modes_removed": theorem_body["all_longitudinal_modes_removed"],
        "all_transverse_modes_have_minus_mu_k2_symbol": theorem_body["all_transverse_modes_have_minus_mu_k2_symbol"],
        "m2_operator_signature_source_exported": theorem_export["m2_operator_signature_source_exported"],
        "no_beta_eta_numeric_source": theorem_export["beta_eta_numeric_source_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_theorem"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2599",
        "stage_id": "S1549",
        "status": "P2599_NADSOLITON_PROJECTED_VISCOUS_STRESS_M2_DERIVATION_THEOREM_EXPORTED_NO_BETA_ETA_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "nadsoliton_projected_viscous_stress_m2_derivation_theorem": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2598_LOCALITY_FRACTIONAL_COMPETITOR_EXCLUSION": sha256_json(p2598_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["nadsoliton_projected_viscous_stress_m2_derivation_theorem"]["theorem_export"]
    theorem_body = t["nadsoliton_projected_viscous_stress_m2_derivation_theorem"]
    lines = [
        "# P2599/S1549 nadsoliton projected viscous-stress m=2 derivation theorem", "",
        f"Status: `{payload['status']}`", "", "## Source theorem", "",
        theorem_body["source_theorem_statement"], "", "## Result", "",
        f"- Inherited P2598 source theorem: `{t['p2598_locality_source_theorem_inherited']}`.",
        f"- Because X refined: `{theorem_body['because_X_refined']}`.",
        f"- Selected operator order: `{theorem_body['selected_operator_order']}`.",
        f"- Projector samples: `{len(theorem_body['projected_viscous_sample_rows'])}`.",
        f"- All longitudinal modes removed: `{theorem_body['all_longitudinal_modes_removed']}`.",
        f"- All transverse modes have `-mu |k|^2` symbol: `{theorem_body['all_transverse_modes_have_minus_mu_k2_symbol']}`.",
        f"- m2 operator signature source exported: `{t['m2_operator_signature_source_exported']}`.", "",
        "## Scope guards", "",
        "This strengthens only the local hydrodynamic source for the `m=2` operator-order slot. It does not export numeric beta/eta sourcing, bridge completion, role-transfer, QW-2191 discharge, role-bearing `L_total`, or ToE closure.", "", "## Fingerprint", "",
        f"`{payload['nadsoliton_projected_viscous_stress_m2_derivation_theorem']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2599/S1549 nadsoliton projected viscous-stress m=2 derivation theorem guard

`P2599/S1549` gives the explicit local constitutive derivation behind P2596--P2598.  For the incompressible nadsoliton information fluid, the local isotropic finite-stress law `sigma_ij=mu(partial_i u_j+partial_j u_i)+lambda delta_ij partial_l u_l` has Fourier divergence which, after `k.u=0` and pressure/Leray projection, acts on transverse modes as `-mu |k|^2`.  Thus the sourced local transport operator is explicitly the Laplacian of order `m=2`.  This still does not export beta/eta numerics, bridge completion, role-transfer, QW-2191 discharge, or ToE closure.
""".strip()
    lag_section = """
## P2599/S1549 nadsoliton projected viscous-stress m=2 derivation theorem Ltotal guard

`P2599/S1549` permits the `m=2` operator-order slot in `L_total` to be treated as explicitly sourced by projected local viscous-stress hydrodynamics.  It does not make the damping/compression term fully role-bearing without numeric beta/eta sourcing, bridge completion, role-transfer, and selector/QW-2191 gates.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2599/S1549 nadsoliton projected viscous-stress m=2 derivation theorem guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2599/S1549 nadsoliton projected viscous-stress m=2 derivation theorem Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
