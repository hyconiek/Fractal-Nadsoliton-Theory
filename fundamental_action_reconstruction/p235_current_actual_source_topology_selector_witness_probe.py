#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F147 = ROOT / "generated" / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
OUT = ROOT / "generated" / "p235_current_actual_source_topology_selector_witness_probe_summary.json"


def main() -> None:
    f147 = json.loads(IN_F147.read_text(encoding="utf-8"))

    passed = (
        f147["input_packet"] == "tau_src_candidate_v1"
        and f147["witness"] == "Pi_sel_src_actual_witness_v1"
        and f147["codomain_packet"] == "Sigma_sel_src_target_v1"
        and f147["actual_selector_witness_exported"] is True
        and f147["chart_bound_preobserver_realization"] is True
        and f147["tau_src_identified_with_s_prelm"] is False
        and f147["observer_role"] == "downstream_only"
        and f147["basis_independence_discharged"] is False
        and f147["qw2191_quotient_safe_discharged"] is False
        and f147["current_selector_closure"] is False
        and f147["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P235",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_BELOW_BASIS_INDEPENDENCE_AND_QW2191_AFTER_P235"
            if passed
            else "P235_FAIL"
        ),
        "input_packet": f147["input_packet"],
        "witness": f147["witness"],
        "codomain_packet": f147["codomain_packet"],
        "support_packet_id": f147["support_packet_id"],
        "observer_role": f147["observer_role"],
        "chart_bound_preobserver_realization": f147["chart_bound_preobserver_realization"],
        "tau_src_identified_with_s_prelm": f147["tau_src_identified_with_s_prelm"],
        "actual_full_source_topology_nontriviality_witness_exported": f147["actual_full_source_topology_nontriviality_witness_exported"],
        "full_source_topology_nontriviality_discharged": f147["full_source_topology_nontriviality_discharged"],
        "actual_selector_witness_exported": f147["actual_selector_witness_exported"],
        "basis_independence_discharged": f147["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f147["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f147["current_selector_closure"],
        "current_global_qw2191_discharge": f147["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
