#!/usr/bin/env python3
"""P3104/S2054: spectral-triple/geometry-interface obstruction audit.

P3103 found a finite Hilbert-like Z12 proxy but no non-imported physical
Hilbert/state-vector source.  P3104 attacks exactly one adjacent standard-
physics interface atom: whether the Z12 Dirichlet/Laplacian branch plus the
legacy entropy number alpha_geo=4 ln 2 can internally source a spectral triple,
metric distance, physical length/action unit, and empirical geometry/readout
interface without importing noncommutative-geometry axioms, continuum manifolds,
measurement apparatus, observed light, selector closure, L_total,
bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit import Z12_SIZE, eigenvalue
from p3103_s2053_hilbert_state_vector_completion_obstruction_audit import OUT as P3103

OUT = GEN / "p3104_s2054_spectral_triple_geometry_interface_obstruction_audit.json"
MD = GEN / "p3104_s2054_spectral_triple_geometry_interface_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = Z12_SIZE
ALPHA_GEO = 4.0 * math.log(2.0)
CONTENT_PATTERNS = {
    "spectral_triple_atom": r"spectral-triple|Dirac|distance formula|geometry-interface",
    "alpha_entropy_unit_atom": r"alpha_geo|4 ln\(2\)|entropy|unit|action",
    "blocked_promotions": r"observed light|L_total|ToE|selector closure|bridge/role-transfer|continuum manifold",
}
GATES = (
    "algebra_representation_exported", "dirac_operator_exported", "compact_resolvent_finite_proxy",
    "bounded_commutators_exported", "distance_formula_exported", "physical_length_unit_exported",
    "alpha_geo_unit_map_exported", "geometry_readout_interface_exported", "nonimported_physical_geometry_source_exported",
)
CANDIDATES = (
    {"id": "z12_diagonal_algebra_representation", "algebra_representation_exported": True, "dirac_operator_exported": False, "compact_resolvent_finite_proxy": True, "bounded_commutators_exported": True, "distance_formula_exported": False, "physical_length_unit_exported": False, "alpha_geo_unit_map_exported": False, "geometry_readout_interface_exported": False, "nonimported_physical_geometry_source_exported": False, "blocker": "finite diagonal C^12 representation is formal and has no sourced Dirac/metric unit"},
    {"id": "laplacian_sqrt_dirac_proxy", "algebra_representation_exported": True, "dirac_operator_exported": True, "compact_resolvent_finite_proxy": True, "bounded_commutators_exported": True, "distance_formula_exported": False, "physical_length_unit_exported": False, "alpha_geo_unit_map_exported": False, "geometry_readout_interface_exported": False, "nonimported_physical_geometry_source_exported": False, "blocker": "sqrt-Laplacian is a computable operator proxy but not an internally sourced physical Dirac operator"},
    {"id": "cycle_graph_connes_distance_proxy", "algebra_representation_exported": True, "dirac_operator_exported": True, "compact_resolvent_finite_proxy": True, "bounded_commutators_exported": True, "distance_formula_exported": True, "physical_length_unit_exported": False, "alpha_geo_unit_map_exported": False, "geometry_readout_interface_exported": False, "nonimported_physical_geometry_source_exported": False, "blocker": "graph distances are dimensionless Z12 combinatorics without physical length calibration"},
    {"id": "alpha_geo_entropy_to_length_unit_candidate", "algebra_representation_exported": True, "dirac_operator_exported": True, "compact_resolvent_finite_proxy": True, "bounded_commutators_exported": True, "distance_formula_exported": True, "physical_length_unit_exported": False, "alpha_geo_unit_map_exported": False, "geometry_readout_interface_exported": False, "nonimported_physical_geometry_source_exported": False, "blocker": "alpha_geo=4 ln 2 is a real scalar/entropy clue but no bit-to-length/action conversion law is exported"},
    {"id": "imported_ncg_spectral_triple_template", "algebra_representation_exported": True, "dirac_operator_exported": True, "compact_resolvent_finite_proxy": True, "bounded_commutators_exported": True, "distance_formula_exported": True, "physical_length_unit_exported": True, "alpha_geo_unit_map_exported": False, "geometry_readout_interface_exported": False, "nonimported_physical_geometry_source_exported": False, "blocker": "passes only by importing noncommutative-geometry/manifold semantics and still lacks alpha_geo unit sourcing"},
    {"id": "imported_empirical_geometry_apparatus_template", "algebra_representation_exported": True, "dirac_operator_exported": True, "compact_resolvent_finite_proxy": True, "bounded_commutators_exported": True, "distance_formula_exported": True, "physical_length_unit_exported": True, "alpha_geo_unit_map_exported": True, "geometry_readout_interface_exported": True, "nonimported_physical_geometry_source_exported": False, "blocker": "empirical geometry/readout passes only by importing rods/clocks/light/apparatus"},
)

def content_grep() -> list[dict[str, Any]]:
    rows=[]
    for lane, pattern in CONTENT_PATTERNS.items():
        proc=subprocess.run(["rg","-n",pattern,"fundamental_action_reconstruction","AGENTS.md","-S"],cwd=REPO,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
        hits=[line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane":lane,"pattern":pattern,"hit_count":len(hits),"sample_hits":hits[:8]})
    return rows

def algebra_rows() -> list[dict[str, Any]]:
    return [{"basis_id":f"e_{i}","diagonal_idempotent":True,"star_closed":True,"represented_on_C12":True,"physical_coordinate_chart_attached":False} for i in range(N)]

def dirac_rows() -> list[dict[str, Any]]:
    rows=[]
    for k in range(N):
        lam=eigenvalue(k); sval=math.sqrt(max(lam,0.0))
        rows.append({"mode":k,"laplacian_eigenvalue":round(lam,12),"sqrt_laplacian_dirac_abs_eigenvalue":round(sval,12),"finite_compact_resolvent_proxy":True,"physical_dirac_unit_attached":False})
    return rows

def distance_rows() -> list[dict[str, Any]]:
    rows=[]
    for i in range(N):
        for j in range(i+1,N):
            graph_dist=min((j-i)%N,(i-j)%N)
            alpha_scaled=ALPHA_GEO*graph_dist
            rows.append({"left":i,"right":j,"z12_graph_distance":graph_dist,"alpha_geo_scaled_distance":round(alpha_scaled,12),"distance_formula_proxy":True,"physical_length_unit_attached":False,"alpha_geo_converted_to_unit":False})
    return rows

def commutator_rows() -> list[dict[str, Any]]:
    rows=[]
    lipschitz_samples=("constant","delta0","ramp","alternating")
    values={"constant":[1.0]*N,"delta0":[1.0 if i==0 else 0.0 for i in range(N)],"ramp":[float(i) for i in range(N)],"alternating":[1.0 if i%2==0 else -1.0 for i in range(N)]}
    for name in lipschitz_samples:
        vals=values[name]
        edge_norm=max(abs(vals[(i+1)%N]-vals[i]) for i in range(N))
        rows.append({"sample_function":name,"finite_commutator_lipschitz_proxy":round(edge_norm,12),"bounded_commutator_proxy":True,"physical_gradient_unit_attached":False})
    return rows

def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate":c["id"],"required_gate":g,"gate_passed":bool(c[g]),"blocking_if_failed":not bool(c[g]),"detail":"passed" if c[g] else c["blocker"]} for c in CANDIDATES for g in GATES]

def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows=[]
    for c in CANDIDATES:
        subset=[r for r in gates if r["candidate"]==c["id"]]
        rows.append({"candidate":c["id"],"passed_gates":sum(r["gate_passed"] for r in subset),"failed_gates":sum(not r["gate_passed"] for r in subset),"accepted_nonimported_geometry_source":all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_geometry_source_exported"])})
    return rows

def build_payload() -> dict[str, Any]:
    p3103=read_json(P3103)
    greps=content_grep(); alg=algebra_rows(); dr=dirac_rows(); dist=distance_rows(); comm=commutator_rows(); gates=gate_rows(); aggs=aggregate(gates)
    accepted=[r for r in aggs if r["accepted_nonimported_geometry_source"]]
    obligations=[
        {"obligation":"read_p3103_next_atom","satisfied":True,"detail":"P3103 selected spectral-triple/geometry-interface as the next atom"},
        {"obligation":"construct_algebra_representation_rows","satisfied":len(alg)==N and all(r["represented_on_C12"] for r in alg),"detail":"finite diagonal algebra representation is explicit"},
        {"obligation":"construct_dirac_spectrum_rows","satisfied":len(dr)==N and all(not r["physical_dirac_unit_attached"] for r in dr),"detail":"sqrt-Laplacian Dirac proxy is explicit but unitless"},
        {"obligation":"construct_distance_formula_rows","satisfied":len(dist)==66 and all(not r["physical_length_unit_attached"] for r in dist),"detail":"Z12 distance and alpha_geo-scaled distances are explicit but dimensionless"},
        {"obligation":"test_alpha_geo_unit_conversion","satisfied":False,"detail":"alpha_geo=4 ln 2 supplies a scalar entropy clue, not an exported bit-to-length/action/unit map"},
        {"obligation":"export_nonimported_physical_geometry_source","satisfied":False,"detail":"0 candidates pass as internal geometry/readout sources"},
    ]
    return {"status":"P3104_SPECTRAL_TRIPLE_GEOMETRY_INTERFACE_OBSTRUCTION_BOUNDED_NO_GO","input_hashes":{"P3103":hashlib.sha256(P3103.read_bytes()).hexdigest() if P3103.exists() else None},"constructed_theoretical_objects":{"content_first_repo_grep":greps,"spectral_triple_geometry_audit_object":{"object":"Z12DirichletSpectralTripleGeometryInterfaceObstructionAudit","source_reused":"P3103 recommendation plus alpha_geo=4 ln 2 entropy/unit clue","alpha_geo":ALPHA_GEO,"required_gates":list(GATES),"candidate_geometry_sources":[c["id"] for c in CANDIDATES],"acceptance_predicate":"algebra representation, Dirac operator, compact-resolvent finite proxy, bounded commutators, distance formula, physical length unit, alpha_geo unit map, geometry/readout interface, and non-imported physical geometry source"},"algebra_representation_rows":alg,"dirac_spectrum_rows":dr,"distance_formula_rows":dist,"commutator_bound_rows":comm,"candidate_gate_rows":gates,"candidate_aggregate_certificate":aggs},"finite_certificate":{"content_grep_lanes":len(greps),"content_grep_hits":sum(r["hit_count"] for r in greps),"p3103_accepted_nonimported_hilbert_sources":p3103.get("finite_certificate",{}).get("accepted_nonimported_hilbert_sources"),"algebra_representation_rows":len(alg),"dirac_spectrum_rows":len(dr),"dirac_rows_with_physical_units":sum(r["physical_dirac_unit_attached"] for r in dr),"distance_formula_rows":len(dist),"distance_rows_with_physical_length_units":sum(r["physical_length_unit_attached"] for r in dist),"distance_rows_with_alpha_geo_unit_conversion":sum(r["alpha_geo_converted_to_unit"] for r in dist),"commutator_bound_rows":len(comm),"commutator_rows_with_physical_gradient_units":sum(r["physical_gradient_unit_attached"] for r in comm),"geometry_candidates":len(CANDIDATES),"required_gates":len(GATES),"candidate_gate_rows":len(gates),"accepted_nonimported_geometry_sources":len(accepted),"proof_obligations":len(obligations),"satisfied_proof_obligations":sum(r["satisfied"] for r in obligations)},"proof_obligations":obligations,"decision":{"bounded_result":"P3104 constructs the requested spectral-triple/geometry-interface obstruction audit.  The Z12 branch supports a finite diagonal algebra on C^12, a sqrt-Laplacian Dirac-like spectrum, bounded finite commutator proxies, graph-distance rows, and alpha_geo-scaled distance rows using alpha_geo=4 ln 2.  These are real mathematical witnesses, but no current artifact exports a physical Dirac unit, a physical length/action unit, an alpha_geo bit-to-length/action conversion law, or an empirical geometry/readout interface.  Therefore no non-imported physical spectral-triple/geometry source is exported.","negative_export_flags":{key:False for key in ["physical_dirac_unit_exported","physical_length_unit_exported","alpha_geo_bit_to_length_action_map_exported","nonimported_spectral_triple_source_exported","geometry_readout_interface_exported","observed_light_interface_exported","physical_hamiltonian_exported","spacetime_eom_exported","QW_2191_discharged","strict_selector_closure_exported","ltotal_exported","bridge_role_transfer_exported","toe_closure_exported"]},"positive_scoped_flags":{"algebra_representation_rows_computed":True,"dirac_spectrum_rows_computed":True,"distance_formula_rows_computed":True,"commutator_bound_rows_computed":True,"alpha_geo_scaled_distance_rows_computed":True},"next_honest_step":"Attack exactly the missing alpha_geo unit-map atom: construct a bounded entropy-to-action/length conversion obstruction audit for alpha_geo=4 ln 2, testing candidate maps from the 16-bit entropy scalar to dimensionful action, length, Hamiltonian time, and detector calibration units.  Do not import Planck units, observed light, quantum measurement postulates, continuum spacetime, selector closure, L_total, bridge/role-transfer, or ToE; accept only an internally sourced bit-to-unit conversion law."}}

def write_markdown(payload: dict[str, Any]) -> None:
    c=payload["finite_certificate"]
    lines=["# P3104/S2054 spectral-triple/geometry-interface obstruction audit","",f"Status: `{payload['status']}`","","## Finite certificate",f"- P3103 accepted non-imported Hilbert sources: `{c['p3103_accepted_nonimported_hilbert_sources']}`",f"- algebra representation rows: `{c['algebra_representation_rows']}`",f"- Dirac spectrum rows: `{c['dirac_spectrum_rows']}`",f"- Dirac rows with physical units: `{c['dirac_rows_with_physical_units']}`",f"- distance formula rows: `{c['distance_formula_rows']}`",f"- distance rows with physical length units: `{c['distance_rows_with_physical_length_units']}`",f"- alpha_geo unit conversions: `{c['distance_rows_with_alpha_geo_unit_conversion']}`",f"- commutator bound rows: `{c['commutator_bound_rows']}`",f"- geometry candidates: `{c['geometry_candidates']}`",f"- candidate gate rows: `{c['candidate_gate_rows']}`",f"- accepted non-imported geometry sources: `{c['accepted_nonimported_geometry_sources']}`",f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`","","## Decision",payload["decision"]["bounded_result"],"","## Recommendation",payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines)+"\n",encoding="utf-8")

def main() -> dict[str, Any]:
    payload=build_payload()
    OUT.write_text(json.dumps(payload,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET,"P3104/S2054 spectral-triple/geometry-interface obstruction audit","## P3104/S2054 spectral-triple/geometry-interface obstruction audit\n\n`P3104/S2054` attacks the P3103-recommended spectral-triple/geometry-interface atom while explicitly testing the only available unit clue, `alpha_geo=4 ln 2`.  It constructs `12` algebra-representation rows, `12` sqrt-Laplacian Dirac-spectrum rows, `66` graph-distance/alpha_geo-scaled distance rows, `4` commutator-bound rows, and a `6 x 9 = 54` candidate gate matrix.  The finite spectral geometry remains formal: no physical Dirac unit, physical length/action unit, alpha_geo bit-to-unit conversion law, empirical geometry/readout interface, observed-light interface, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT,"P3104/S2054 spectral geometry and alpha_geo unit map remain unsourced","## P3104/S2054 spectral geometry and alpha_geo unit map remain unsourced\n\n`P3104/S2054` confirms that the Z12 branch supports a finite algebra representation, a sqrt-Laplacian Dirac-like proxy, bounded commutators, graph-distance witnesses, and alpha_geo-scaled metric rows.  A Lagrangian/EOM, physical geometry, observed-light, or empirical-readout reading still needs an internally sourced bit-to-length/action conversion law for `alpha_geo=4 ln 2` plus a physical geometry/readout interface; imported NCG axioms, Planck units, rods/clocks, and apparatus templates do not supply those strict sources.\n")
    append_once(AGENTS,"Current spectral-triple/geometry-interface obstruction guardrail (P3104/S2054, 2026-06-25)","## Current spectral-triple/geometry-interface obstruction guardrail (P3104/S2054, 2026-06-25)\n\n- P3104 follows the P3103 recommendation and audits one standard-physics interface atom: a spectral-triple/geometry-interface source for the Z12 Dirichlet/Laplacian branch, including the `alpha_geo=4 ln 2` entropy/unit clue.\n- The finite audit computes `12` algebra-representation rows, `12` Dirac-spectrum rows, `66` distance/alpha_geo-scaled rows, `4` commutator-bound rows, and `54` candidate gate rows; `0` candidates export an internal non-imported physical geometry/source law.\n- Do not promote finite diagonal algebra, sqrt-Laplacian Dirac proxies, bounded finite commutators, graph distances, alpha_geo-scaled distances, imported NCG templates, Planck-unit templates, or empirical apparatus templates to physical geometry, action/length units, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one alpha_geo entropy-to-action/length unit-map obstruction audit, unless a genuinely new spectral-geometry source theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
