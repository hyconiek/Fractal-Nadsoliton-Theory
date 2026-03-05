#!/usr/bin/env python3
"""
QW-2117: K_total locality/operator-status audit gate.

Purpose:
- classify whether current strict implementation uses K_total as
  index-space mixing (local in spacetime) or as a spacetime-nonlocal kernel,
- provide auditable evidence from strict chain code/docs.
"""

from __future__ import annotations

import ast
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2117_ktotal_locality_operator_audit.json"
OUT_MD = ROOT / "RAPORT_QW2117_KTOTAL_LOCALITY_OPERATOR_AUDIT.md"


STRICT_CODE_FILES = [
    "QW_2063_DERIVATIONAL_RECONSTRUCTION_SHARED_FLAVOR_BASIS.py",
    "QW_2069_FULL_SM_GR_DERIVATION_PACKAGE.py",
    "QW_2071_SM_GR_FULL_PRECISION_CLOSURE_GATE.py",
    "QW_2097_CKM_CP_TARGET_REFINEMENT_GATE.py",
    "QW_2115_GRAVITY_HIERARCHY_STRICT_BRIDGE_GATE.py",
]

CANONICAL_DOC_FILES = [
    "TOE_FINAL_DOCUMENTATION.tex",
    "RELEASE_5.md",
    "RELEASE_5_TEXTBOOK_EN_PL.md",
]

NONLOCAL_TOKEN_PATTERNS = [
    r"\bfftconvolve\b",
    r"\bconvolve\b",
    r"\bintegral\s*d\^?4y\b",
    r"\bnp\.fft\b",
    r"\bscipy\.signal\b",
    r"\bK\(\s*x\s*-\s*y\s*\)",
]


def load_text(name: str) -> str:
    return (ROOT / name).read_text(encoding="utf-8", errors="ignore")


def ast_kernel_signature_index_only(src: str) -> Dict[str, object]:
    tree = ast.parse(src)
    out: Dict[str, object] = {
        "kernel_fn_exists": False,
        "kernel_fn_args": [],
        "kernel_fn_uses_spacetime_symbols": False,
    }
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "kernel_fn":
            out["kernel_fn_exists"] = True
            args = [a.arg for a in node.args.args]
            out["kernel_fn_args"] = args
            names = {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}
            out["kernel_fn_uses_spacetime_symbols"] = bool(names.intersection({"x", "y", "t", "r"}))
            break
    return out


def find_nonlocal_hits(file_names: List[str]) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for name in file_names:
        txt = load_text(name)
        hits: List[str] = []
        for pat in NONLOCAL_TOKEN_PATTERNS:
            if re.search(pat, txt, flags=re.IGNORECASE):
                hits.append(pat)
        if hits:
            out[name] = hits
    return out


def main() -> None:
    src2063 = load_text("QW_2063_DERIVATIONAL_RECONSTRUCTION_SHARED_FLAVOR_BASIS.py")
    ksig = ast_kernel_signature_index_only(src2063)

    uses_cyclic_q_distance = "cyclic_distance_matrix(q, q" in src2063
    uses_spacetime_convolution_form = bool(re.search(r"\bK\(\s*x\s*-\s*y\s*\)", src2063))

    strict_nonlocal_hits = find_nonlocal_hits(STRICT_CODE_FILES)
    canonical_nonlocal_hits = find_nonlocal_hits(CANONICAL_DOC_FILES)

    formal_draft = load_text("FORMALNE_WYPROWADZENIA_BRAKUJACYCH_POZYCJI_FIN_V5_DRAFT.md")
    has_explicit_a_vs_b_locality_clause = (
        "Przypadek A: operator lokalny w czasoprzestrzeni" in formal_draft
        and "Przypadek B: operator nielokalny spacetime" in formal_draft
    )

    flags = {
        "kernel_fn_exists": bool(ksig["kernel_fn_exists"]),
        "kernel_fn_args_distance_only_signature": bool(
            ksig["kernel_fn_args"] == ["d", "omega", "phi", "beta", "eta"]
        ),
        "kernel_fn_no_spacetime_symbols": bool(not ksig["kernel_fn_uses_spacetime_symbols"]),
        "flavor_hamiltonian_uses_cyclic_q_distance": bool(uses_cyclic_q_distance),
        "no_spacetime_convolution_form_in_qw2063": bool(not uses_spacetime_convolution_form),
        "strict_chain_has_no_nonlocal_tokens": bool(len(strict_nonlocal_hits) == 0),
        "formal_a_vs_b_locality_clause_present": bool(has_explicit_a_vs_b_locality_clause),
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    implementation_local = bool(
        flags["kernel_fn_exists"]
        and flags["kernel_fn_args_distance_only_signature"]
        and flags["kernel_fn_no_spacetime_symbols"]
        and flags["flavor_hamiltonian_uses_cyclic_q_distance"]
        and flags["no_spacetime_convolution_form_in_qw2063"]
        and flags["strict_chain_has_no_nonlocal_tokens"]
    )
    operator_status = (
        "INDEX_SPACE_MIXING_OPERATOR_LOCAL_IN_SPACETIME_IMPLEMENTATION"
        if implementation_local
        else "UNRESOLVED_OR_NONLOCAL_IN_IMPLEMENTATION"
    )

    verdict = (
        "KTOTAL_LOCALITY_OPERATOR_AUDIT_PASS_IMPLEMENTATION_LOCAL"
        if implementation_local
        else "KTOTAL_LOCALITY_OPERATOR_AUDIT_FAIL_OR_UNRESOLVED"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "strict_code_files": STRICT_CODE_FILES,
            "canonical_doc_files": CANONICAL_DOC_FILES,
            "formal_draft": "FORMALNE_WYPROWADZENIA_BRAKUJACYCH_POZYCJI_FIN_V5_DRAFT.md",
        },
        "kernel_signature_audit": ksig,
        "implementation_evidence": {
            "uses_cyclic_q_distance": bool(uses_cyclic_q_distance),
            "uses_spacetime_convolution_form": bool(uses_spacetime_convolution_form),
        },
        "nonlocal_token_hits": {
            "strict_code": strict_nonlocal_hits,
            "canonical_docs": canonical_nonlocal_hits,
        },
        "operator_status": operator_status,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "FREEZE_OPERATOR_STATUS_DECLARATION_IN_CANONICAL_STRICT_PROTOCOL"
            if implementation_local
            else "RESOLVE_OPERATOR_STATUS_AND_RERUN_QW2117"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2117: KTOTAL LOCALITY OPERATOR AUDIT",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        f"- operator_status: `{operator_status}`",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Nonlocal token hits",
            f"- strict_code hits: `{len(strict_nonlocal_hits)}`",
            f"- canonical_docs hits: `{len(canonical_nonlocal_hits)}`",
            "",
            "## Artifact",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2117] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2117] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2117] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

