#!/usr/bin/env python3
"""Build the archival Release 10.20 SHA-256 manifest."""

from pathlib import Path
import hashlib
import json


ROOT = Path(__file__).resolve().parent
TARGET = ROOT / "FIN_Programs_217_229_Release_10_20_Manifest.json"
FIG = ROOT / "FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Figures"

PRIMARY = [
    "FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Monograph.pdf",
    "FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Monograph.md",
    "FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Monograph.tex",
    "FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Results.json",
    "fin_programs_217_229_operational_section_spectral_tomography.py",
    "test_fin_programs_217_229_operational_section_spectral_tomography.py",
    "fin_programs_217_229_to_latex.py",
    "fin_programs_217_229_build_manifest.py",
    "FIN_Programs_217_229_Dirichlet_Core.lean",
    "FIN_Programs_217_229_Phase_Torsion_Core.lean",
    "FIN_Programs_217_229_Arb_Execution_Attestation.json",
    "FIN_Programs_217_229_Lean_Build_Record.json",
    "FIN_Programs_217_229_Phase_Formalization_Record.json",
    "FIN_Programs_217_229_Independent_Custody_Bundle_Request.json",
    "FIN_Programs_217_229_Semigroup_Holdout_Preregistration.json",
    "RELEASE_10_20_PROGRAMS_217_229.md",
]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    paths = [ROOT / name for name in PRIMARY]
    paths.extend(sorted(FIG.glob("*.png")))
    missing = [str(path.relative_to(ROOT)) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"release files missing: {missing}")
    entries = [
        {
            "path": str(path.relative_to(ROOT)),
            "bytes": path.stat().st_size,
            "sha256": digest(path),
        }
        for path in paths
    ]
    manifest = {
        "manifest_id": "FIN-RELEASE-10.20-PROGRAMS-217-229",
        "release": "10.20",
        "version": "1.0.0",
        "date": "2026-07-27",
        "creator": "Żuchowski, Krzysztof",
        "orcid": "0009-0002-0909-3613",
        "license": "CC BY 4.0",
        "programs": list(range(217, 230)),
        "recommended_programs": list(range(230, 243)),
        "test_count": 82,
        "test_status": "PASS",
        "pdf_pages": 51,
        "lean_dirichlet_core_compiled": True,
        "lean_phase_core_compiled": True,
        "lean_general_library_compiled": False,
        "formal_arb_certificate_completed": False,
        "firecrawl_used": False,
        "external_web_used": False,
        "external_dataset_used": False,
        "external_physical_validation_claimed": False,
        "files": entries,
        "claim_boundaries": {
            "QW_2191_discharged": False,
            "strict_selector_exported": False,
            "strict_central_state_exported": False,
            "target_independent_eta_source_exported": False,
            "canonical_physical_unit_exported": False,
            "legacy_strict_bridge_completed": False,
            "legacy_role_transfer_started": False,
            "L_total_or_ToE_claimed": False,
            "external_bundle_admitted": False,
            "external_prediction_executed": False,
        },
    }
    TARGET.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(TARGET.name)


if __name__ == "__main__":
    main()
