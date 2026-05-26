from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.py",
    ROOT / "p2133_s1083_strict_cmp2_real_extension_merge_contract.py",
    ROOT / "p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.py",
    ROOT / "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.py",
]
OUT = ROOT / "generated" / "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.json"


class TestP2135StrictCmp2MergedRealBlockVariantStabilityAudit(unittest.TestCase):
    def test_p2135_exports_variant_stability_audit(self) -> None:
        for s in SCRIPTS:
            subprocess.run([sys.executable, str(s)], check=True)

        d = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(d["schema_version"], "p2135_s1085_v1")
        self.assertEqual(d["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertIn(d["result_kind"], {
            "PASS_STRICT_CMP2_MERGED_REAL_BLOCK_VARIANT_STABILITY_AUDIT_WITH_TRACE",
            "OPEN_STRICT_CMP2_MERGED_REAL_BLOCK_VARIANT_STABILITY_AUDIT_BLOCKED",
        })

        checks = d["gatekeeper_checks"]
        self.assertFalse(checks["full_d3_covariance_transport_proven"])
        self.assertFalse(checks["full_cutkosky_closure_proven"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
