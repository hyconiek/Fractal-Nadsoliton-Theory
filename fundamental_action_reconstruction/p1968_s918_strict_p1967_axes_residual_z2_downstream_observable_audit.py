#!/usr/bin/env python3
"""P1968 S918 strict P1967 axes residual-Z2 downstream observable audit.

This is the next honest step after P1967.  It uses the strict Shannon
axis-only source instantiated as Delta_sel tensors and asks which downstream
claims are already safe under the remaining residual sign u -> -u.

The result is deliberately split:

* projective/quadratic Release-7 OS observables are safe, reusing P709/N706;
* sign-sensitive directed/global selector closure is not safe and remains
  obstructed, reusing P687/N687/N707.

No global QW-2191 discharge, admissible S_sel_int, physical sign datum, or ToE
closure is claimed.
"""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1968_s918_strict_p1967_axes_residual_z2_downstream_observable_audit.json"

SOURCE_FILES = [
    "P709_CURRENT_STRICT_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDIT_PROBE.md",
    "N706_CURRENT_STRICT_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_PACKAGE_THEOREM.md",
    "P687_CURRENT_STRICT_T173_GLOBAL_EDGE_SIGN_COHERENCE_SOLVABILITY_AUDIT_PROBE.md",
    "N687_CURRENT_STRICT_T173_GLOBAL_EDGE_SIGN_COHERENCE_OBSTRUCTION_BOUNDARY_THEOREM.md",
    "N707_CURRENT_STRICT_T173_PREVIOUS_METHODOLOGY_SURVIVAL_AND_GLOBAL_GAP_BOUNDARY_THEOREM.md",
    "N493_CURRENT_FIRST_STRICT_QW2191_RESIDUAL_Z2_SIGN_FLIP_GAUGE_EQUIVALENCE_THEOREM.md",
]

KEYWORDS = [
    "residual",
    "sign",
    "gauge",
    "irrelevance",
    "directed",
    "physical orientation",
    "QW-2191",
    "ToE closure",
    "projective",
    "observable",
    "Z2",
]


def load_generated(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def file_digest(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(obj: object) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode("utf-8")).hexdigest()


def grep_evidence() -> list[dict[str, Any]]:
    evidence: list[dict[str, Any]] = []
    for name in SOURCE_FILES:
        path = ROOT / name
        if not path.exists():
            evidence.append({"path": name, "present": False, "matches": []})
            continue
        matches = []
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            lower = line.lower()
            if any(keyword.lower() in lower for keyword in KEYWORDS):
                matches.append({"line": lineno, "text": line.strip()[:260]})
            if len(matches) >= 10:
                break
        evidence.append({"path": name, "present": True, "sha256": file_digest(path), "matches": matches})
    return evidence


def symbolic_z2_invariance() -> dict[str, Any]:
    s = sp.symbols("s", real=True)
    h11, h12, h22, u1, u2, l1, l2 = sp.symbols("h11 h12 h22 u1 u2 l1 l2", real=True)
    H = sp.Matrix([[h11, h12], [h12, h22]])
    u = sp.Matrix([u1, u2])
    ell = sp.Matrix([l1, l2])

    projector = u * u.T
    signed_projector = (s * u) * (s * u).T
    rayleigh = (u.T * H * u)[0]
    signed_rayleigh = ((s * u).T * H * (s * u))[0]
    linear = (ell.T * u)[0]
    signed_linear = (ell.T * (s * u))[0]

    substitutions = {s**2: 1}
    return {
        "projector_difference_under_s_squared_1": sp.sstr(sp.expand(signed_projector - projector).subs(substitutions)),
        "quadratic_rayleigh_difference_under_s_squared_1": sp.sstr(sp.expand(signed_rayleigh - rayleigh).subs(substitutions)),
        "linear_channel_difference_under_s_squared_1": sp.sstr(sp.expand(signed_linear - linear).subs(substitutions)),
        "classification": {
            "projector_or_ray_observables": "Z2-even / sign-gauge-invariant",
            "quadratic_forms": "Z2-even / sign-gauge-invariant",
            "linear_or_directed_output_channels": "Z2-odd unless an extra sign convention or physical sign datum is supplied",
        },
    }


def numeric_delta_projector_checks(p1967: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in p1967.get("pair_rows", []):
        mat = np.array(row["Delta_sel_pair_m"], dtype=float)
        vals, vecs = np.linalg.eigh(mat)
        u = vecs[:, 1]
        projector = np.outer(u, u)
        signed_projector = np.outer(-u, -u)
        rows.append(
            {
                "pair": row["pair"],
                "delta_sel_gap": row["Delta_sel_gap"],
                "projector_sign_flip_inf_norm": float(np.linalg.norm(projector - signed_projector, ord=np.inf)),
                "quadratic_value_sign_flip_diff": float(abs(((-u).T @ mat @ (-u)) - (u.T @ mat @ u))),
                "linear_sum_channel_sign_flip_diff_example": float(abs(np.sum(-u) - np.sum(u))),
                "projector_and_quadratic_even_pass": bool(
                    np.linalg.norm(projector - signed_projector, ord=np.inf) < 1e-12
                    and abs(((-u).T @ mat @ (-u)) - (u.T @ mat @ u)) < 1e-12
                ),
            }
        )
    return rows


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    p1967 = load_generated("p1967_s917_strict_shannon_axis_source_to_delta_sel_tensor_map.json")
    p709 = load_generated("p709_current_strict_release_7_os_residual_sign_gauge_irrelevance_audit_probe_summary.json")
    p687 = load_generated("p687_current_strict_t173_global_edge_sign_coherence_solvability_audit_probe_summary.json")
    n681 = load_generated("n681_current_strict_t173_directed_output_sign_lift_obstruction_boundary_theorem_summary.json")
    n687 = load_generated("n687_current_strict_t173_global_edge_sign_coherence_obstruction_boundary_theorem_summary.json")

    symbolic = symbolic_z2_invariance()
    numeric_rows = numeric_delta_projector_checks(p1967)

    p1967_axis_pass = p1967.get("machine_checks", {}).get("strict_axis_only_source_pass") is True
    p709_sign_safe = p709.get("sign_ok") is True and p709.get("status") == "PASS_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDITED"
    p687_obstructed = p687.get("sign_system_solvable") is False and p687.get("status") == "PASS_SIGN_COHERENCE_OBSTRUCTED"
    n681_blocks_directed = str(n681.get("status", "")).endswith("NO_FALSE_PASS") or "OBSTRUCTION" in str(n681.get("status", ""))
    n687_blocks_global = str(n687.get("status", "")).endswith("NO_FALSE_PASS") or "OBSTRUCTION" in str(n687.get("status", ""))
    numeric_even_pass = all(row["projector_and_quadratic_even_pass"] for row in numeric_rows)

    out = {
        "packet_id": "P1968",
        "stage_id": "S918",
        "status": "P1967_AXES_DOWNSTREAM_Z2_AUDIT_PASS_FOR_PROJECTIVE_OS__SIGN_SENSITIVE_GLOBAL_SELECTOR_STILL_BLOCKED",
        "route": "strict_only_projective_observable_audit_after_shannon_axis_source",
        "legacy_bridge_used": False,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "sympy": sp.__version__,
        "repo_grep_evidence": grep_evidence(),
        "input_generated_hashes": {
            "p1967_sha256": digest(p1967),
            "p709_summary_sha256": digest(p709),
            "p687_summary_sha256": digest(p687),
            "n681_summary_sha256": digest(n681),
            "n687_summary_sha256": digest(n687),
        },
        "symbolic_z2_parity_lemma": symbolic,
        "numeric_delta_sel_projector_checks": numeric_rows,
        "imported_existing_machine_audits": {
            "p1967_axis_only_source_pass": p1967_axis_pass,
            "p709_release7_os_sign_gauge_irrelevance_pass": p709_sign_safe,
            "p709_patterns_tested": p709.get("n_patterns_tested"),
            "p709_max_abs_diffs": {
                "p694_m2": p709.get("max_abs_diff_p694_m2_under_sign_patterns"),
                "p696_channel_m2": p709.get("max_abs_diff_p696_channel_m2_under_sign_patterns"),
                "p696_offdiag_to_diag_ratio": p709.get("max_abs_diff_p696_offdiag_to_diag_ratio_under_sign_patterns"),
                "p696_offblock_max_fro": p709.get("max_abs_diff_p696_offblock_max_fro_under_sign_patterns"),
            },
            "p687_global_edge_sign_coherence_obstructed": p687_obstructed,
            "n681_directed_output_sign_lift_blocks_physical_sign_claim": n681_blocks_directed,
            "n687_global_edge_sign_coherence_blocks_global_directed_selector": n687_blocks_global,
        },
        "machine_checks": {
            "numeric_p1967_projector_and_quadratic_even_pass": numeric_even_pass,
            "release7_projective_os_safe_under_residual_z2": bool(p1967_axis_pass and p709_sign_safe and numeric_even_pass),
            "global_sign_sensitive_selector_closure_pass": False,
            "global_sign_sensitive_selector_blocked_by_existing_audits": bool(p687_obstructed and n681_blocks_directed and n687_blocks_global),
        },
        "theorem_export": {
            "positive_projective_statement": "For P1967 axes, projectors/rays and quadratic Release-7 OS observables are residual-Z2 even; P709 already audits all 2^5 sign patterns for the downstream OS values.",
            "negative_directed_statement": "Linear/directed output channels remain residual-Z2 sensitive unless an additional sign convention or strict physical sign datum is supplied; P687/N687 keep global edge sign coherence obstructed.",
            "not_promoted_to": [
                "directed/sign-sensitive physical orientation datum in strict core",
                "admissible S_sel_int",
                "kernel-alone/global QW-2191 discharge",
                "global T173 discharge",
                "ToE closure",
            ],
        },
        "false_pass_guard": "P1968 licenses projective/quadratic OS continuation under the P1967 axes, not global sign-sensitive selector closure.",
        "next_honest_step": "If strict ToE closure needs only projective/quadratic OS observables, cite P1968 as the Z2-safety bridge; if a directed observable is required, build a new global sign-provider or prove a no-go theorem for that directed channel.",
        "lay_explanation": "Mamy już osie z P1967. P1968 sprawdza, co zależy od odwrócenia strzałki na osi. Wielkości typu długość, projektor i energia kwadratowa się nie zmieniają, więc operacyjny tor OS jest bezpieczny. Ale wielkości kierunkowe, które widzą zwrot strzałki, nadal wymagają dodatkowego globalnego znaku.",
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
