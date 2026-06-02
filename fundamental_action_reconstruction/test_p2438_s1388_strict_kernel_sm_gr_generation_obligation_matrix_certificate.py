import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2438_s1388_strict_kernel_sm_gr_generation_obligation_matrix_certificate.py"
OUT = ROOT / "generated" / "p2438_s1388_strict_kernel_sm_gr_generation_obligation_matrix_certificate.json"
MD = ROOT / "generated" / "p2438_s1388_strict_kernel_sm_gr_generation_obligation_matrix_certificate.md"
P2437 = ROOT / "generated" / "p2437_s1387_legacy_kernel_physical_value_derivation_methodology_audit_certificate.json"


class P2438StrictKernelSmGrGenerationObligationMatrixCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2437.exists():
            subprocess.run([sys.executable, str(ROOT / "p2437_s1387_legacy_kernel_physical_value_derivation_methodology_audit_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_kernel_sm_gr_generation_obligation_matrix_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_matrix_size(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2438")
        self.assertEqual(self.payload["stage_id"], "S1388")
        self.assertEqual(
            self.payload["status"],
            "PASS_STRICT_KERNEL_SM_GR_GENERATION_OBLIGATION_MATRIX_ALL_TARGETS_BLOCKED_NO_CLOSURE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["obligation_count"], 8)
        self.assertEqual(self.theorem["target_count"], 6)
        self.assertEqual(len(self.cert["finite_obligation_matrix"]["target_rows"]), 6)
        self.assertEqual(len(self.cert["finite_obligation_matrix"]["obligation_rows"]), 8)

    def test_targets_blocked_and_masks(self) -> None:
        self.assertEqual(self.theorem["ready_target_count_now"], 0)
        self.assertEqual(self.theorem["ready_targets_now"], [])
        self.assertTrue(self.theorem["all_targets_blocked_now"])
        self.assertEqual(self.theorem["current_discharge_mask"], 0)
        self.assertEqual(self.theorem["full_discharge_mask"], 255)
        self.assertGreaterEqual(self.theorem["minimum_missing_obligations_for_any_target"], 5)
        self.assertEqual(self.theorem["maximum_missing_obligations_for_target"], 8)

    def test_inherited_blockers_and_hard_limits(self) -> None:
        self.assertTrue(self.theorem["legacy_value_route_retired_as_generator"])
        self.assertTrue(self.theorem["p1421_qw2191_missing_inherited"])
        self.assertTrue(self.theorem["p1646_scaffold_open_inherited"])
        self.assertTrue(self.theorem["p1705_qg_required_inherited"])
        self.assertTrue(self.theorem["p1981_background_open_inherited"])
        self.assertTrue(self.theorem["p2437_no_strict_value_generator_inherited"])
        self.assertFalse(self.theorem["strict_sm_gr_generation_theorem_exported"])
        self.assertFalse(self.theorem["strict_observable_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("strict-kernel SM/GR generation", MD.read_text(encoding="utf-8"))
        self.assertIn("P2438/S1388", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2438/S1388", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
