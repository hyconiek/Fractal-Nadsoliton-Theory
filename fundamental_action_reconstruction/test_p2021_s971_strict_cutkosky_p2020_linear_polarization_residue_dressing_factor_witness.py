import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2021_s971_strict_cutkosky_p2020_linear_polarization_residue_dressing_factor_witness.py"
OUT = ROOT / "generated" / "p2021_s971_strict_cutkosky_p2020_linear_polarization_residue_dressing_factor_witness.json"


def test_p2021_proxy_residue_transport_is_rejected_as_dressed_cutkosky_input():
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["result_kind"] == "OPEN_PROXY_RESIDUE_TRANSPORT_SANITY_ONLY_NOT_P1953_ADMISSIBLE"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert data["legacy_bridge_used"] is False
    assert data["admissible_for_p1953_dressed_interface"] is False
    assert all(data["local_algebra_sanity_checks"].values())

    admissibility = data["admissibility_checks_for_dressed_cutkosky"]
    assert admissibility == {
        "loop_derived_from_L_total": False,
        "same_scheme_as_P1953_MSbar_B1_seed_locked": False,
        "DiscM_common_basis_exported": False,
        "DiscM_minus_CutSum_simplified_evaluated": False,
        "BRST_physical_state_projector_exported": False,
        "proxy_not_promoted_to_dressed_residue": True,
    }

    factor = data["proxy_residue_factor_sandbox"]
    assert factor["Z_proxy(s)"] == "(51*s + 50)/(50*(s + 1))"
    assert factor["R_proxy(s)=Z_proxy(s)^2"] == "(51*s + 50)**2/(2500*(s + 1)**2)"
    assert factor["domain"] == "s > 0 real"
    assert "not accepted as P1953" in factor["admissibility_warning"]

    transported = data["proxy_transport_linear_polarization_matrices_over_kappa2_Zgauge2"]
    assert transported["basis_order"] == ["plus", "cross"]
    assert transported["accepted_as_dressed_residue"] is False
    assert transported["no_identical_symmetry"] == [
        ["(51*s + 50)**2/(2500*pi*(s + 1)**2)", "0"],
        ["0", "(51*s + 50)**2/(2500*pi*(s + 1)**2)"],
    ]
    assert transported["identical_final_state"] == [
        ["(51*s + 50)**2/(5000*pi*(s + 1)**2)", "0"],
        ["0", "(51*s + 50)**2/(5000*pi*(s + 1)**2)"],
    ]

    non_update = data["p1953_contract_non_update"]
    assert non_update["M_dressed_common_basis"] == "OPEN_NOT_EXPORTED__P2021_PROXY_FACTOR_REJECTED_AS_DRESSED_INPUT"
    assert non_update["DiscM_common_basis"] == "OPEN_NOT_EVALUATED"
    assert non_update["DiscM_minus_CutSum_simplified"] == "OPEN_NOT_EVALUATED"
    assert "explicitly rejects that proxy" in data["false_pass_guard"]

    grid = data["numeric_grid_certificate"]
    assert grid["max_exact_numeric_eigenvalue_l2_error"] < 1e-14
    assert all(row["positive_semidefinite_numeric"] for row in grid["rows"])
