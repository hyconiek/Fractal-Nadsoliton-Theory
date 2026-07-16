"""P3170/S2120: S_+ source-obligation hitting-set theorem.

P3169 closed the non-binary Lambda_origin profile branch as receiver-only and
recommended pivoting to the other P3168 hard branch: a genuinely scale-charged
strict S_+ value coupled to Omega_M/K_dim.  P3170 does not invent that value.
Instead it constructs the proof object that tells exactly what is missing.

The audit reads the P3162 S_+ candidate gate matrix, builds a finite blocker
hypergraph (candidate routes -> failed source obligations), and exhaustively
computes minimal hitting sets.  This turns the informal request for S_+ into a
sharp theorem: every current route is blocked by the single common atom
`strict_nadsoliton_source_exported`; the least-repair admissible route is the
`new_strict_scale_charged_datum_schema`, which still needs exactly the paired
atoms `nonzero_value_exported` and `strict_nadsoliton_source_exported`.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
GEN.mkdir(exist_ok=True)
OUT = GEN / "p3170_s2120_s_plus_source_obligation_hitting_set.json"
MD = GEN / "p3170_s2120_s_plus_source_obligation_hitting_set.md"
AGENTS = REPO / "AGENTS.md"
SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
INPUTS = {
    "P3161": GEN / "p3161_s2111_omega_scale_positive_torsor_source_law_audit.json",
    "P3162": GEN / "p3162_s2112_s_plus_scale_charged_source_datum_intake_audit.json",
    "P3167": GEN / "p3167_s2117_s_plus_monomial_source_exhaustion.json",
    "P3168": GEN / "p3168_s2118_post_p3167_no_strict_unit_state_map_certificate.json",
    "P3169": GEN / "p3169_s2119_ternary_origin_datum_exhaustive_audit.json",
    "SUMMARY_GROK": REPO / "SUMMARY_GROK.md",
}
ADMISSIBLE_SCHEMA = "new_strict_scale_charged_datum_schema"


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    except Exception:
        return {}


def rg(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        ["rg", "-n", pattern, "AGENTS.md", "SUMMARY_GROK.md", "fundamental_action_reconstruction", "-g", "*.md", "-g", "*.json", "-g", "*.py"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = [line for line in proc.stdout.splitlines() if line]
    return {"count": len(lines), "samples": lines[:20]}


def append_once(path: Path, marker: str, text: str) -> None:
    old = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in old:
        path.write_text(old.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def candidate_failures(p3162: dict[str, Any]) -> dict[str, set[str]]:
    failures: dict[str, set[str]] = {}
    for row in p3162["constructed_theoretical_objects"]["candidate_gate_rows"]:
        failures.setdefault(row["candidate"], set())
        if not row["passed"]:
            failures[row["candidate"]].add(row["gate"])
    return failures


def minimal_hitting_sets(failures: dict[str, set[str]]) -> list[tuple[str, ...]]:
    universe = sorted(set().union(*failures.values()))
    routes = list(failures.values())
    out: list[tuple[str, ...]] = []
    for size in range(1, len(universe) + 1):
        for combo in itertools.combinations(universe, size):
            cset = set(combo)
            if all(cset & route for route in routes):
                if not any(set(prev).issubset(cset) for prev in out):
                    out.append(combo)
        if out:
            return out
    return out


def route_table(failures: dict[str, set[str]]) -> list[dict[str, Any]]:
    rows = []
    for candidate, failed in sorted(failures.items()):
        rows.append(
            {
                "candidate": candidate,
                "failed_obligations": sorted(failed),
                "failed_count": len(failed),
                "is_admissible_schema": candidate == ADMISSIBLE_SCHEMA,
                "accepted_S_plus_source": False,
            }
        )
    return rows


def repair_frontier(failures: dict[str, set[str]]) -> list[dict[str, Any]]:
    rows = []
    for candidate, failed in sorted(failures.items(), key=lambda item: (len(item[1]), item[0])):
        import_or_selector_blocked = bool({"not_external_unit_import", "selector_bridge_ltotal_toe_free"} & failed)
        dimensionless_blocked = bool({"weight_plus_one_under_Rpos", "not_dimensionless_invariant"} & failed)
        rows.append(
            {
                "candidate": candidate,
                "missing_atoms": sorted(failed),
                "missing_count": len(failed),
                "recommended_for_repair": candidate == ADMISSIBLE_SCHEMA,
                "reason": (
                    "least-repair admissible schema; needs actual nonzero strict source export"
                    if candidate == ADMISSIBLE_SCHEMA
                    else "blocked by import/selector lane" if import_or_selector_blocked
                    else "blocked by dimensionless invariant lane" if dimensionless_blocked
                    else "formal or receiver-only route without strict source value"
                ),
            }
        )
    return rows


def payload() -> dict[str, Any]:
    p3162 = load(INPUTS["P3162"])
    failures = candidate_failures(p3162)
    hits = minimal_hitting_sets(failures)
    routes = route_table(failures)
    repairs = repair_frontier(failures)
    common_failures = sorted(set.intersection(*failures.values()))
    admissible_missing = sorted(failures[ADMISSIBLE_SCHEMA])
    return {
        "status": "P3170_S_PLUS_SOURCE_OBLIGATION_HITTING_SET_THEOREM",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: load(path).get("status") or ("markdown_or_missing_json" if path.exists() else "missing") for key, path in INPUTS.items()},
        "repo_grep": {
            "s_plus_frontier": rg(r"S_\+|scale-charged|Omega_M|K_dim|positive torsor|weight-one"),
            "strict_source_atom": rg(r"strict_nadsoliton_source_exported|strict source|source law|nonzero strict"),
            "closed_imports": rg(r"Planck|apparatus|selector replay|bridge completion|role transfer|L_total|ToE"),
        },
        "constructed_theoretical_objects": {
            "S_plus_source_obligation_hypergraph": "candidate S_+ routes as hyperedges of failed strict-source obligations from P3162",
            "minimal_universal_blocker_cut": [list(item) for item in hits],
            "common_failed_atoms": common_failures,
            "route_failure_table": routes,
            "least_repair_frontier": repairs,
        },
        "finite_certificate": {
            "candidate_routes": len(failures),
            "distinct_failed_obligations": len(set().union(*failures.values())),
            "route_failure_rows": len(routes),
            "minimal_hitting_sets": len(hits),
            "minimal_hitting_set_size": len(hits[0]) if hits else 0,
            "common_failed_atoms_count": len(common_failures),
            "admissible_schema_missing_atoms": len(admissible_missing),
            "accepted_S_plus_sources": 0,
        },
        "finite_theorem": {
            "name": "P3170_T1_strict_source_atom_is_the_minimal_S_plus_blocker",
            "statement": "In the P3162 S_+ candidate gate matrix, the singleton atom {strict_nadsoliton_source_exported} is a minimal universal hitting set: every current candidate route fails it.  The least-repair admissible route is not an existing receiver but the typed schema new_strict_scale_charged_datum_schema, which still lacks exactly nonzero_value_exported and strict_nadsoliton_source_exported.  Therefore no current artifact exports S_+, and the next proof move must construct that paired source/value atom rather than scan more dimensionless receivers.",
        },
        "decision": {
            "bounded_result": "P3170 sharpens the S_+ frontier to a finite blocker theorem instead of another inventory: all current routes share the missing strict source atom.",
            "next_honest_step": "Construct exactly one formula-level strict source law Source_S_plus(nadsoliton)->s in V_chi with chi(c)=c, prove s > 0 and nonzero, and prove its coupling to Omega_M/K_dim.  Do not use dimensionless invariants, Planck/apparatus calibration, selector replay, bridge role transfer, or formal Omega_M self-source.  If that law cannot be supplied, stop at the P3168-P3170 no-new-live-frontier certificate or write CA+SA only as non-strict conditioning.",
            "negative_export_flags": {
                "S_plus_source_exported": False,
                "Omega_M_fixed": False,
                "K_dim_functor_exported": False,
                "unit_source_exported": False,
                "selector_closure_exported": False,
                "bridge_completion_exported": False,
                "role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
        },
    }


def write_md(data: dict[str, Any]) -> None:
    cert = data["finite_certificate"]
    lines = [
        "# P3170/S2120 S_plus source-obligation hitting-set theorem",
        "",
        f"Status: `{data['status']}`",
        "",
        "## Constructed objects",
        "- `S_plus_source_obligation_hypergraph`: candidate `S_+` routes as failed-obligation hyperedges from P3162.",
        "- `minimal_universal_blocker_cut`: finite hitting-set computation over all current `S_+` routes.",
        "- `least_repair_frontier`: ordered list of routes by missing atoms; the admissible schema is isolated from receiver/import lanes.",
        "",
        "## Finite certificate",
    ]
    for key, value in cert.items():
        lines.append(f"- `{key}`: `{value}`")
    lines += [
        "",
        "## Theorem",
        f"`{data['finite_theorem']['name']}`: {data['finite_theorem']['statement']}",
        "",
        "## Next honest step",
        data["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    data = payload()
    OUT.write_text(json.dumps(data, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_md(data)
    append_once(
        SHEET,
        "P3170/S2120 S_plus source-obligation hitting-set theorem",
        """## P3170/S2120 S_plus source-obligation hitting-set theorem

`P3170/S2120` converts the remaining `S_+` branch into a finite blocker hypergraph using the `P3162` candidate gate matrix.  Across `12` candidate routes, the singleton atom `{strict_nadsoliton_source_exported}` is a minimal universal hitting set: every route fails it.  The least-repair admissible route is `new_strict_scale_charged_datum_schema`, which still needs exactly `nonzero_value_exported` and `strict_nadsoliton_source_exported`.  No `S_+`, `Omega_M`, `K_dim`, unit-bearing action, selector closure, bridge completion, role transfer, `L_total`, or ToE closure is exported.""",
    )
    append_once(
        DRAFT,
        "P3170/S2120 S_plus blocker theorem leaves action units conditional",
        """## P3170/S2120 S_plus blocker theorem leaves action units conditional

`P3170/S2120` shows that the remaining unit branch is blocked by a precise missing atom: a nonzero strict nadsoliton source value in the weight-one `S_+` representation coupled to `Omega_M/K_dim`.  Until that source law is supplied, all action units, mass units, nonproxy action density, EOM terms, and `L_total` uses remain conditional or imported, not strict exports.""",
    )
    append_once(
        AGENTS,
        "Current S_plus source-obligation hitting-set guardrail (P3170/S2120, 2026-07-16)",
        """## Current S_plus source-obligation hitting-set guardrail (P3170/S2120, 2026-07-16)

- P3170 builds the finite `S_+` source-obligation hypergraph from the P3162 gate matrix: `12` candidate routes and their failed strict-source obligations.
- The singleton atom `{strict_nadsoliton_source_exported}` is a minimal universal blocker: every current `S_+` route fails it.
- The least-repair admissible route is the typed `new_strict_scale_charged_datum_schema`, which still lacks exactly `nonzero_value_exported` and `strict_nadsoliton_source_exported`.
- Do not continue dimensionless receiver scans, Planck/apparatus imports, selector replay, bridge role-transfer, formal `Omega_M` self-sourcing, or profile-origin detours as `S_+`, unit, `L_total`, or ToE closure.
- The next proof-grade move must provide one formula-level source law `Source_S_plus(nadsoliton) -> s in V_chi`, prove `s > 0`, and prove coupling to `Omega_M/K_dim`; otherwise preserve the P3168-P3170 no-new-live-frontier certificate or write CA+SA only as explicit non-strict conditioning.
""",
    )
    return data


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
