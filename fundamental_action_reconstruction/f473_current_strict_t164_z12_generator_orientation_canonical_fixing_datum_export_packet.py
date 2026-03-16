#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-16"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_Z12 = GENERATED / "z_12_v1_group.json"
IN_AUT = GENERATED / "aut_z12_v1_units_group.json"

OUT_OBJ = GENERATED / "kappa_z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1.json"
OUT_SUMMARY = (
    GENERATED
    / "f473_current_strict_t164_z12_generator_orientation_canonical_fixing_datum_export_packet_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing = [str(p.relative_to(REPO)) for p in (IN_Z12, IN_AUT) if not p.exists()]
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "F473",
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": missing,
                    "expected": [
                        "generated/z_12_v1_group.json (F329/N450)",
                        "generated/aut_z12_v1_units_group.json (F331/N453)",
                    ],
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    z12 = load_json(IN_Z12)
    aut = load_json(IN_AUT)
    if str(z12.get("object") or "") != "Z_12_v1":
        raise SystemExit(
            json.dumps(
                {"stage": "F473", "status": "INVALID_Z12_OBJECT", "got": z12.get("object"), "no_false_pass": True},
                ensure_ascii=True,
            )
        )
    if str(aut.get("object") or "") != "Aut_Z12_v1":
        raise SystemExit(
            json.dumps(
                {"stage": "F473", "status": "INVALID_AUT_OBJECT", "got": aut.get("object"), "no_false_pass": True},
                ensure_ascii=True,
            )
        )

    # Premise-based fixing: pick generator 1, thus successor suc(k) = k+1 mod 12.
    g_fix = 1
    suc_fix_table = [(k + 1) % 12 for k in range(12)]

    # Aut(Z12) elements are units {1,5,7,11} acting by k ↦ (u*k) mod 12.
    # The only unit preserving g_fix=1 is u=1.
    aut_preserving_fix = [1]

    obj = {
        "object": "Kappa_Z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1",
        "status": "actual_exported_strict_provenance_datum_object__premise_based",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export one explicit fixing datum that canonically fixes (in a declared scope) a generator/orientation of the "
            "typed strict group object Z_12_v1, so future directed/sign-sensitive constructions can be made slotless without "
            "smuggling an untracked generator choice. By N462, no Aut(Z_12)-invariant fixing exists; therefore this fixing datum "
            "is exported as an explicit premise-based strict provenance object (no false pass)."
        ),
        "typed_prerequisites": {
            "Z_12_v1": str(IN_Z12.relative_to(REPO)),
            "Aut_Z12_v1": str(IN_AUT.relative_to(REPO)),
        },
        "declared_domain": {
            "group_object": "Z_12_v1",
            "aut_group_object": "Aut_Z12_v1",
            "index_set": "I_12_v1 = {0,1,...,11}",
            "scope_note": (
                "This fixing datum is intended only for declared directed continuations (e.g. T171). "
                "It is not an Aut(Z_12)-invariant construction; it explicitly restricts the admissible symmetry to the subgroup "
                "preserving the fix."
            ),
        },
        "definition": {
            "D_fix_v1": {
                "generator_fixed": g_fix,
                "successor_map": "suc_fix(k) := (k+1) mod 12",
                "successor_table": suc_fix_table,
                "meaning": "Fix an oriented 12-cycle structure induced by generator 1.",
            }
        },
        "reduced_symmetry": {
            "Aut_Z12_v1_full": (aut.get("carrier") or {}).get("elements"),
            "Aut_Z12_v1_preserving_D_fix_v1": aut_preserving_fix,
            "note": (
                "The datum reduces the admissible Aut(Z_12) ambiguity to the subgroup preserving the fixed generator/orientation. "
                "Here the stabilizer of generator 1 inside Aut(Z_12) is the identity only."
            ),
        },
        "compatibility_with_T163_T164": {
            "T164_target": "Kappa_Z12_generator_orientation_canonical_fixing_datum_target_v1",
            "T163_branch": "T163(4a) canonical-fixing-datum route (nontrivial phase embeddings require tracked symmetry breaking).",
            "note": (
                "This object is a premise-based strict provenance fixing datum intended to satisfy the T164 acceptance pattern: "
                "it exports an explicit D_fix_v1 and makes the reduced admissible symmetry explicit. It does not claim Aut-invariant canonicity."
            ),
        },
        "provenance": {
            "classification": "strict_source_upgraded",
            "method": "explicit_strict_side_premise_generator_orientation_fix",
            "premise_id": "P_z12_generator_orientation_fix_v1",
            "boundary_reference": "N462 (no Aut(Z12)-invariant canonical generator choice from typed structure alone)",
            "note": (
                "Premise-based strict provenance fixing datum: it explicitly breaks/restricts Aut(Z_12) and therefore must be cited "
                "by any directed/sign-sensitive downstream object. This does not derive any theta angles and does not claim selector closure."
            ),
        },
        "hard_limits": [
            "Does not claim Aut(Z_12)-invariant canonicity (ruled out by N462).",
            "Does not discharge T163 nor export any canonical phase embedding.",
            "Does not export theta_1/theta_2 and does not imply any sigma-int -> theta upgrade.",
            "Does not claim strict-core selector closure nor global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    OUT_OBJ.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    summary = {
        "stage": "F473",
        "status": "PASS_EXPORTED_T164_FIXING_DATUM_OBJECT_PREMISE_BASED",
        "exported": str(OUT_OBJ.relative_to(REPO)),
        "t164_discharged": True,
        "premise_based": True,
        "no_false_pass": True,
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

