from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p1950_s900_strict_renormalization_exact_integration.py",
    ROOT / "p2094_s1044_strict_b1_quotient_renormalization_rank_repair.py",
    ROOT / "p2095_s1045_strict_b1_gb_derived_channel_certificate.py",
    ROOT / "p2096_s1046_strict_b1_renormalization_closure_contract.py",
    ROOT / "p2097_s1047_strict_b1_quotient_closure_stability_minigrid.py",
    ROOT / "p2098_s1048_strict_precutkosky_readiness_contract.py",
    ROOT / "p2099_s1049_strict_u1_same_scheme_lock_witness.py",
    ROOT / "p2100_s1050_strict_u2_phase_space_quadrature_witness.py",
    ROOT / "p2101_s1051_strict_u3_residue_positivity_uncertainty_witness.py",
    ROOT / "p2102_s1052_strict_task2_entry_gate_summary.py",
    ROOT / "p2103_s1053_strict_dressed_discontinuity_backend_import_contract.py",
    ROOT / "p2104_s1054_strict_d1_dressed_pole_residue_source_object.py",
]
OUT = ROOT / "generated" / "p2104_s1054_strict_d1_dressed_pole_residue_source_object.json"


class TestP2104StrictD1DressedPoleResidueSourceObject(unittest.TestCase):
    def test_p2104_exports_d1_source_object(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2104_s1054_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_D1_DRESSED_POLE_RESIDUE_SOURCE_OBJECT_WITH_TRACE__SEED_LEVEL",
        )

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["d1_source_object_present"])
        self.assertTrue(checks["seed_level_only"])
        self.assertFalse(checks["full_dressed_backend_import_executed"])
        self.assertFalse(checks["full_cutkosky_closure_proven"])
        self.assertFalse(checks["c3_theorem_proven"])


if __name__ == "__main__":
    unittest.main()
