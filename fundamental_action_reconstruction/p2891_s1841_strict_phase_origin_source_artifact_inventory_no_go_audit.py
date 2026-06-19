#!/usr/bin/env python3
"""P2891/S1841: strict phase/origin source artifact inventory no-go audit.

P2890 showed that Fourier power is too invariant and Fourier phase needs an
external phase/origin pin.  The next honest move is not another carrier replay;
it is an intake gate for an actually exported strict phase/origin source artifact
with nonconventional sign/phase and a coupling theorem to the 9/5 variational
density.

This packet scans current generated JSON artifacts for that exact missing object.
It is a current-artifact certificate, not a future impossibility theorem.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2890 = GEN / "p2890_s1840_fourier_phase_source_law_9_over_5_no_go_audit.json"
OUT = GEN / "p2891_s1841_strict_phase_origin_source_artifact_inventory_no_go_audit.json"
MD = GEN / "p2891_s1841_strict_phase_origin_source_artifact_inventory_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

PHASE_SOURCE_TOKENS = (
    "phase_source",
    "phase_origin",
    "origin_source",
    "translation_breaking",
    "symmetry_breaking",
    "chiral",
    "pseudoscalar",
    "nonconventional_sign",
    "nonconventional_phase",
    "phase_pin",
)
COUPLING_TOKENS = (
    "9_over_5",
    "9:5",
    "9/5",
    "variational_density",
    "variational_chain_rule",
    "localized_action_density",
    "unit_bearing",
)
POSITIVE_VERBS = ("export", "accepted", "accept", "select", "derive", "coupl", "source")
NEGATIVE_TOKENS = ("negative_export_flags", "without", "missing", "blocked", "no_go", "no-go", "obstruction", "requires", "need")
NEGATIVE_PREFIXES = ("no_", "not_")


def json_files() -> list[Path]:
    return sorted(path for path in GEN.glob("*.json") if path.is_file() and path != OUT)


def walk(value: Any, prefix: str = "") -> list[tuple[str, Any]]:
    rows: list[tuple[str, Any]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            child_prefix = f"{prefix}.{key}" if prefix else str(key)
            rows.extend(walk(child, child_prefix))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            rows.extend(walk(child, f"{prefix}[{index}]"))
    else:
        rows.append((prefix, value))
    return rows


def text_of(path: str, value: Any) -> str:
    return f"{path} {value}".lower()


def has_any(text: str, tokens: tuple[str, ...]) -> bool:
    return any(token in text for token in tokens)


def is_positive_export_path(path: str, value: Any) -> bool:
    if value is not True:
        return False
    lowered = path.lower()
    leaf = lowered.rsplit(".", 1)[-1]
    if leaf.startswith(NEGATIVE_PREFIXES) or "_no_" in leaf or "_not_" in leaf:
        return False
    if any(token in lowered for token in NEGATIVE_TOKENS):
        return False
    return any(verb in lowered for verb in POSITIVE_VERBS)


def inventory_records() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for file_path in json_files():
        try:
            payload = read_json(file_path)
        except Exception as exc:
            records.append({"file": str(file_path.relative_to(ROOT)), "read_error": repr(exc)})
            continue
        phase_hits = []
        coupling_hits = []
        positive_phase_booleans = []
        positive_coupling_booleans = []
        for path, value in walk(payload):
            text = text_of(path, value)
            if has_any(text, PHASE_SOURCE_TOKENS):
                phase_hits.append({"path": path, "value": str(value)[:220]})
                if is_positive_export_path(path, value):
                    positive_phase_booleans.append({"path": path, "value": value})
            if has_any(text, COUPLING_TOKENS):
                coupling_hits.append({"path": path, "value": str(value)[:220]})
                if is_positive_export_path(path, value):
                    positive_coupling_booleans.append({"path": path, "value": value})
        if phase_hits or coupling_hits:
            records.append(
                {
                    "file": str(file_path.relative_to(ROOT)),
                    "status": payload.get("status"),
                    "phase_hit_count": len(phase_hits),
                    "coupling_hit_count": len(coupling_hits),
                    "positive_phase_boolean_count": len(positive_phase_booleans),
                    "positive_coupling_boolean_count": len(positive_coupling_booleans),
                    "candidate_coupled_positive_export": bool(positive_phase_booleans and positive_coupling_booleans),
                    "sample_phase_hits": phase_hits[:8],
                    "sample_coupling_hits": coupling_hits[:8],
                    "positive_phase_booleans": positive_phase_booleans[:8],
                    "positive_coupling_booleans": positive_coupling_booleans[:8],
                }
            )
    return records


def build_payload(p2890: dict[str, Any]) -> dict[str, Any]:
    files = json_files()
    records = inventory_records()
    coupled_positive_records = [record for record in records if record.get("candidate_coupled_positive_export")]
    facts = {
        "p2890_rechecked": p2890.get("status") == "P2890_FOURIER_PHASE_SOURCE_LAW_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE",
        "generated_json_inventory_nonempty": len(files) > 0,
        "phase_or_coupling_terms_found": len(records) > 0,
        "no_coupled_positive_phase_origin_9_over_5_export_found": len(coupled_positive_records) == 0,
        "p2890_obligation_remains_unsatisfied": len(coupled_positive_records) == 0,
    }
    return {
        "status": "P2891_STRICT_PHASE_ORIGIN_SOURCE_ARTIFACT_INVENTORY_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2890": sha(P2890)},
        "strict_phase_origin_source_artifact_inventory_no_go_audit": {
            "input_status_rechecked": p2890.get("status"),
            "candidate_class": "current generated-artifact inventory for a strict phase/origin source with nonconventional sign/phase and coupling theorem to the 9/5 variational density",
            "generated_json_file_count_excluding_self": len(files),
            "relevant_record_count": len(records),
            "coupled_positive_export_record_count": len(coupled_positive_records),
            "sample_relevant_records": records[:24],
            "proof_certificate": {
                "inventory_rule": "Scan generated JSON artifacts for phase/origin/translation-breaking/chiral/pseudoscalar source terms and 9/5 variational-density coupling terms; count a candidate only when positive non-negative booleans for both sides occur in the same artifact.",
                "finite_result": "Relevant source and coupling language exists in current artifacts, mostly as negative guardrails or uncoupled scoped results, but no generated artifact exports the coupled strict phase/origin source required after P2890.",
                "sourcehood_step": "A future move must supply a new explicit source artifact and coupling theorem rather than another Fourier/carrier/orbit replay or a bare inventory hit.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_strict_phase_origin_source_artifact": False,
            "exports_nonconventional_phase_or_sign_source": False,
            "exports_coupling_to_9_over_5_variational_density": False,
            "exports_nonimported_9_over_5_variational_chain_rule": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "strict_phase_origin_source_artifact_exported": False,
                "nonconventional_phase_or_sign_source_exported": False,
                "coupling_to_9_over_5_variational_density_exported": False,
                "nonimported_9_over_5_variational_chain_rule_exported": False,
                "localized_action_density_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2891 performs the exact post-P2890 intake gate across current generated artifacts.  It finds phase/source/coupling language, but no artifact with coupled positive exports for a strict phase/origin source and a 9/5 variational-density theorem.  The P2890 obligation remains unsatisfied on current artifacts.",
            "next_honest_step": "Do not replay Fourier-power signatures, Fourier-character phase pins, P2888/P2890 carrier/orbit representatives, bounded coefficients, C12-neutral unit measures, or inventory hits as strict sourcehood.  A next proof-grade move must introduce one new explicit strict phase/origin source formula with computed nonconventional sign/phase and a coupling theorem to the 9/5 variational density, or pivot to a genuinely different typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["strict_phase_origin_source_artifact_inventory_no_go_audit"]
    lines = [
        "# P2891/S1841 strict phase/origin source artifact inventory no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Phase/origin source artifact inventory",
        f"- generated JSON files scanned, excluding self: `{audit['generated_json_file_count_excluding_self']}`",
        f"- relevant record count: `{audit['relevant_record_count']}`",
        f"- coupled positive export record count: `{audit['coupled_positive_export_record_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2890))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2891/S1841 strict phase/origin source artifact inventory no-go audit",
        "## P2891/S1841 strict phase/origin source artifact inventory no-go audit\n\n"
        "`P2891/S1841` performs the exact post-`P2890` intake gate across generated artifacts for a strict phase/origin source with nonconventional sign/phase and a coupling theorem to the `9/5` variational density.  Relevant phase/source/coupling language exists, but no artifact has coupled positive exports for both the strict source and the `9/5` variational-density theorem.  No strict translation-breaking phase source, nonimported `9/5` variational chain rule, localized action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2891/S1841 strict phase/origin source artifact inventory `L_total` guard",
        "## P2891/S1841 strict phase/origin source artifact inventory `L_total` guard\n\n"
        "`P2891/S1841` adds no strict action term.  Current generated artifacts do not export a coupled strict phase/origin source plus `9/5` variational-density theorem, so no unit-bearing localized action density or variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source is licensed.\n",
    )
    append_once(
        AGENTS,
        "Current strict phase/origin source artifact inventory guardrail (P2891/S1841, 2026-06-19)",
        "## Current strict phase/origin source artifact inventory guardrail (P2891/S1841, 2026-06-19)\n\n"
        "- P2891 scans current generated artifacts for the exact post-P2890 missing object: a strict phase/origin source with nonconventional sign/phase and a coupling theorem to the `9/5` variational density.\n"
        "- Relevant phase/source/coupling language exists, but no artifact exports coupled positive source and `9/5` variational-density booleans satisfying the obligation.\n"
        "- Do not promote Fourier-power signatures, Fourier-character phase pins, P2888/P2890 carrier/orbit representatives, bounded coefficients, `C12`-neutral unit measures, or inventory hits to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must introduce one new explicit strict phase/origin source formula with computed nonconventional sign/phase and coupling theorem to the `9/5` variational density, pivot to a genuinely different typed object, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
