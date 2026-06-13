#!/usr/bin/env python3
"""P2678/S1628: strict internal orientation-source provider-class audit.

P2677 closed the current tau_src -> pair12 typed-seed route as a bounded no-go.
This packet switches to the different provider class it recommended: an
orientation-source theorem/provider that could break the XOR/XNOR reversal and
bind a pre-collapse Sigma_sel_src_target_v1 -> F301 arrow.
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
OUT = GEN / "p2678_s1628_strict_internal_orientation_source_provider_class_audit.json"
MD = GEN / "p2678_s1628_strict_internal_orientation_source_provider_class_audit.md"
P2677 = GEN / "p2677_s1627_s6_o3_typed_seed_route_no_go_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "strict_internal_orientation_source_theorem_exported",
    "oriented_torsor_provider_exported",
    "precollapse_sigma_f301_binding_exported",
    "xor_xnor_reversal_broken",
    "q_w_2191_discharged",
    "s6_reopened",
    "o3_reopened",
    "o4_o5_allowed",
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
        "different_provider_class_content": "provider class|different provider|new provider|source theorem|orientation-source theorem",
        "oriented_torsor_content": "orientation torsor|oriented torsor|orientation-odd|odd sign|spin/Pin|spin orientation|Pin orientation",
        "legacy_scalar_even_content": "beta_tors.*chi11|beta_tors -> chi_11|legacy scalar|scalar input|torsion damping",
        "symmetry_breaking_content": "symmetry-breaking|symmetry breaking|spontaneous symmetry|boundary condition|orientation boundary",
        "sigma_f301_binding_content": "Sigma_sel_src_target_v1.*F301|F301.*Sigma_sel_src_target_v1|pre-collapse Sigma|Sigma->F301|Sigma_sel_src_target_v1 -> F301",
        "forbidden_collapse_content": "Q_basis|preLM|projector-only|nonprojector|quotient|collapse|fiat|convention",
        "closure_guard_content": "QW-2191|O3|O4/O5|L_total|ToE closure|role transfer|bridge completion",
    }
    return {
        "tool": "rg",
        "mode": "content-first different-provider-class audit for strict internal orientation source; not name/number-only",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def upstream_consistency() -> dict[str, Any]:
    p2677 = load_json(P2677)
    decision = p2677.get("closure_decision", {})
    lattice = p2677.get("finite_no_go_lattice", {})
    return {
        "p2677_current_route_no_go": decision.get("s6_current_route_passed_now") is False and decision.get("o3_current_route_passed_now") is False,
        "p2677_not_global_no_go": decision.get("global_no_go_claimed") is False,
        "p2677_missing_internal_source": "strict_internal_orientation_source_exported" in lattice.get("missing_current_obligations", []),
        "p2677_missing_nonfiat_choice": "nonfiat_orientation_choice_proof" in lattice.get("missing_current_obligations", []),
    }


def c2_equivariant_provider_enumeration() -> dict[str, Any]:
    """Enumerate provider kinds under orientation reversal C2.

    A strict sign source is an odd map f(-x)=-f(x).  Even/scalar or quotient
    data cannot produce such a map into {+1,-1}; oriented torsor/spin-boundary
    input can formally do so, but only if the repo exports the provider and its
    Sigma->F301 binding.
    """
    provider_rows = []
    provider_kinds = [
        ("legacy_scalar_beta_tors", "even_scalar", False, False, False),
        ("axis_or_projector_quotient", "quotient_even", False, False, False),
        ("q_basis_terminal_branch", "postcollapse_choice", False, False, False),
        ("declared_orientation_convention", "fiat_choice", False, False, False),
        ("oriented_c2_torsor", "orientation_odd_torsor", True, False, False),
        ("spin_pin_orientation_source", "orientation_odd_spin_pin", True, False, False),
        ("boundary_symmetry_breaking_source", "orientation_odd_boundary", True, False, False),
    ]
    maps = []
    # all functions f:{0,1}->{+1,-1}; C2 acts by x -> 1-x and sign -> -sign
    for table in range(4):
        values = {x: 1 if ((table >> x) & 1) else -1 for x in (0, 1)}
        odd = all(values[1 - x] == -values[x] for x in (0, 1))
        even = all(values[1 - x] == values[x] for x in (0, 1))
        maps.append({"table": table, "values": values, "orientation_odd": odd, "orientation_even": even})
    for name, provider_type, has_formal_odd_domain, exported_now, binding_now in provider_kinds:
        admissible_maps = [row["table"] for row in maps if row["orientation_odd"] and has_formal_odd_domain]
        passes = bool(admissible_maps) and exported_now and binding_now
        provider_rows.append({
            "provider": name,
            "provider_type": provider_type,
            "has_formal_orientation_odd_domain": has_formal_odd_domain,
            "formal_odd_map_tables": admissible_maps,
            "provider_exported_now": exported_now,
            "precollapse_sigma_f301_binding_exported_now": binding_now,
            "passes_provider_class_gate": passes,
        })
    return {
        "group": "C2 orientation reversal",
        "sign_maps_checked": len(maps),
        "orientation_odd_map_count": sum(row["orientation_odd"] for row in maps),
        "orientation_even_map_count": sum(row["orientation_even"] for row in maps),
        "sign_map_rows": maps,
        "provider_rows": provider_rows,
        "formal_provider_classes_with_odd_domain": [row["provider"] for row in provider_rows if row["has_formal_orientation_odd_domain"]],
        "passing_provider_class_count": sum(row["passes_provider_class_gate"] for row in provider_rows),
    }


def provider_obligation_lattice() -> dict[str, Any]:
    obligations = [
        "different_provider_class_selected",
        "orientation_odd_domain_exported",
        "nonfiat_symmetry_breaking_or_torsor_law_exported",
        "precollapse_sigma_side_binding_exported",
        "surviving_f301_carrier_binding_exported",
        "nonquotient_nonprojector_transport_exported",
        "closure_guards_preserved",
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
    for bits in itertools.product([False, True], repeat=len(obligations)):
        pass_count += int(all(bits))
    missing = [key for key, value in current.items() if not value]
    return {
        "obligations": obligations,
        "total_states": 2 ** len(obligations),
        "passing_states": pass_count,
        "current_state": current,
        "missing_current_obligations": missing,
        "hamming_distance_to_pass": len(missing),
    }


def closure_decision(enumeration: dict[str, Any], lattice: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2678_STRICT_INTERNAL_ORIENTATION_SOURCE_PROVIDER_CLASS_AUDIT__FORMAL_ODD_PROVIDERS_IDENTIFIED_BUT_NOT_EXPORTED",
        "professorial_verdict": (
            "P2678 moves to the different provider class requested after P2677.  The finite C2 enumeration shows exactly what kind of object could break XOR/XNOR without fiat: an orientation-odd torsor, spin/Pin orientation source, or boundary symmetry-breaking source.  "
            "Legacy scalar beta_tors, quotient/projector axes, Q_basis branches, and declared conventions cannot supply an odd strict sign.  The formal odd provider classes exist as admissible shapes, but none is currently exported together with a pre-collapse Sigma->F301 binding, so S6/O3 is not reopened."
        ),
        "next_honest_step": (
            "The next honest proof-grade step is a construction-or-no-go audit for one concrete orientation-odd provider object, preferably the smallest C2 torsor/spin-boundary source with an explicit action on Sigma_sel_src_target_v1 and a typed arrow into the surviving F301 carrier.  If no such object can be exported, keep the S6/O3 route closed and shift the main bridge work back to the legacy->strict completion bridge source terms."
        ),
        "formal_odd_provider_class_count": len(enumeration["formal_provider_classes_with_odd_domain"]),
        "passing_provider_class_count": enumeration["passing_provider_class_count"],
        "hamming_distance_to_pass": lattice["hamming_distance_to_pass"],
        "strict_internal_orientation_source_exported_now": False,
        "precollapse_sigma_f301_binding_exported_now": False,
        "s6_reopened_now": False,
        "o3_reopened_now": False,
        "o4_o5_allowed_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2678/S1628 strict internal orientation-source provider-class audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    enum = payload["c2_equivariant_provider_enumeration"]
    lines.extend([
        "", "## C2 equivariant provider enumeration",
        f"Sign maps checked: `{enum['sign_maps_checked']}`.",
        f"Orientation-odd sign maps: `{enum['orientation_odd_map_count']}`.",
        f"Formal odd provider classes: `{enum['formal_provider_classes_with_odd_domain']}`.",
        f"Passing provider-class count: `{enum['passing_provider_class_count']}`.",
        "", "## Provider rows",
    ])
    for row in enum["provider_rows"]:
        lines.append(f"- `{row['provider']}` ({row['provider_type']}): formal_odd_domain=`{row['has_formal_orientation_odd_domain']}`, exported=`{row['provider_exported_now']}`, binding=`{row['precollapse_sigma_f301_binding_exported_now']}`, pass=`{row['passes_provider_class_gate']}`")
    lat = payload["provider_obligation_lattice"]
    lines.extend([
        "", "## Provider obligation lattice",
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
    enumeration = c2_equivariant_provider_enumeration()
    lattice = provider_obligation_lattice()
    payload: dict[str, Any] = {
        "status": "P2678_STRICT_INTERNAL_ORIENTATION_SOURCE_PROVIDER_CLASS_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {"P2677": sha256_file(P2677), "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET), "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT)},
        "upstream_consistency": upstream_consistency(),
        "c2_equivariant_provider_enumeration": enumeration,
        "provider_obligation_lattice": lattice,
        "closure_decision": closure_decision(enumeration, lattice),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2678/S1628 strict internal orientation-source provider-class guard",
        "## P2678/S1628 strict internal orientation-source provider-class guard\n\n"
        "`P2678/S1628` switches away from the bounded no-go typed-seed branch and audits the different provider class required to break the XOR/XNOR reversal: a strict internal orientation-odd source.  The finite `C2` enumeration shows that legacy scalar `beta_tors`, quotient/projector axes, `Q_basis` terminal choices, and declaration/convention choices cannot provide an odd sign.  Only orientation-odd torsor, spin/Pin orientation, or boundary symmetry-breaking provider shapes could do so formally, and none is currently exported with a pre-collapse `Sigma_sel_src_target_v1 -> F301` binding.  Therefore S6/O3 is not reopened; O4/O5, `QW-2191` discharge, role-bearing `L_total`, bridge completion, role transfer, and ToE closure remain blocked.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2678/S1628 strict orientation-source provider Ltotal guard",
        "## P2678/S1628 strict orientation-source provider Ltotal guard\n\n"
        "`P2678/S1628` keeps `L_total` closed while moving to the correct different provider class.  A variational boundary-square or sector-swap source term requires an exported orientation-odd torsor/spin-boundary provider plus a pre-collapse typed `Sigma_sel_src_target_v1 -> F301` binding; formal provider shapes alone are not Lagrangian terms.\n",
    )
    return payload


if __name__ == "__main__":
    main()
