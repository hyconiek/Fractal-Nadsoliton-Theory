import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2440_s1390_current_strict_tuple_coefficient_replay_rank_certificate.py"
OUT = ROOT / "generated" / "p2440_s1390_current_strict_tuple_coefficient_replay_rank_certificate.json"
MD = ROOT / "generated" / "p2440_s1390_current_strict_tuple_coefficient_replay_rank_certificate.md"
P2439 = ROOT / "generated" / "p2439_s1389_strict_coefficient_source_consistency_audit_certificate.json"


class P2440CurrentStrictTupleCoefficientReplayRankCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2439.exists():
            subprocess.run([sys.executable, str(ROOT / "p2439_s1389_strict_coefficient_source_consistency_audit_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["current_strict_tuple_coefficient_replay_rank_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_replay_size(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2440")
        self.assertEqual(self.payload["stage_id"], "S1390")
        self.assertEqual(
            self.payload["status"],
            "PASS_CURRENT_STRICT_TUPLE_COEFFICIENT_REPLAY_PHASE_NULL_NO_GENERATOR",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["coefficient_count"], 14)
        self.assertEqual(len(self.theorem["current_tuple_replayed_coefficients"]), 14)
        self.assertEqual(len(self.theorem["jacobian_numeric"]), 14)

    def test_rank_phase_null_and_inverse(self) -> None:
        self.assertEqual(self.theorem["jacobian_parameter_order"], ["omega", "phi", "beta", "eta", "A"])
        self.assertEqual(self.theorem["jacobian_real_rank"], 4)
        self.assertTrue(self.theorem["phi_column_zero"])
        self.assertFalse(self.theorem["phase_parameter_recovered_by_replay"])
        self.assertTrue(self.theorem["local_inverse"]["local_inverse_pass_on_recovered_four_parameters"])
        self.assertEqual(self.theorem["local_inverse"]["unrecovered_parameters"], ["phi"])

    def test_quarantine_and_hard_limits(self) -> None:
        self.assertEqual(self.theorem["p1664_current_tuple_coefficient_comparison"]["row_count"], 14)
        self.assertGreater(self.theorem["p1664_current_tuple_coefficient_comparison"]["changed_count"], 0)
        self.assertTrue(self.theorem["p1562_closure_flags_conflict_with_current_no_closure_state"])
        self.assertTrue(self.theorem["p2439_no_acceptable_generator_inherited"])
        self.assertFalse(self.theorem["strict_kernel_to_coefficient_map_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("coefficient replay rank certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2440/S1390", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2440/S1390", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
