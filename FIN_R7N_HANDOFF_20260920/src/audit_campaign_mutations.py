from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import json, copy
ROOT=Path(__file__).resolve().parents[1]

def load(p): return json.load(open(ROOT/p))

def main():
    phase=load('results/R7N-046_phase_exhaustion_final.json')
    paudit=load('results/R7N-046_phase_exhaustion_audit.json')
    pract=load('results/R7N-019_024_practical_lane_final.json')
    target=load('results/R7N-017_target_p_implication.json')
    union=load('results/safe_union_v2_audited.json')
    claims=load('CLAIMS.json')
    checks={}
    # Phase mutations were executed by independent phase audit.
    mut=paudit['checks']['mutations']['detail']
    checks['phase_delete_leaf_rejected']=bool(mut['delete_leaf_rejected'])
    checks['phase_shift_boundary_rejected']=bool(mut['shift_boundary_rejected'])
    checks['phase_gradient_sign_rejected']=bool(mut['damaged_gradient_sign_rejected'])
    checks['phase_enlarged_collar_rejected']=bool(mut['enlarged_collar_rejected'])
    checks['phase_fixture_mutation_rejected']=bool(mut['source_amplitude_mutation_rejected'])
    checks['phase_bad_symmetry_rejected']=bool(mut['odd_translation_symmetry_rejected'])
    # Formula replay must be complete.
    checks['K16_formula_replay_complete']=phase['formula_level_replay']['K16']['complete'] and not phase['formula_level_replay']['K16']['failed']
    checks['K20_formula_replay_complete']=phase['formula_level_replay']['K20']['complete'] and not phase['formula_level_replay']['K20']['failed']
    # Target P must reject global PASS while residual remains.
    residual=pract['final_certified_partial']['residual_volume_fraction_of_compact_hull']
    checks['target_P_global_pass_rejected_with_residual']=residual>0 and pract['global_pass'] is False and claims['TARGET_P']['status'].startswith('OPEN')
    # A tau0 result must not be treated as sigma: thresholds are distinct and sigma<tau0.
    checks['sigma_tau0_swap_rejected']=target['sigma_lt_tau0'] and target['tau0']=='67/250' and claims['TARGET_S']['status'].startswith('OPEN')
    # Navigation buffer must not be proof input.
    checks['navigation_buffer_rejected']=union.get('navigation_buffer_used') is False
    # Canonical safe union count.
    direct=len(union['direct_domains']); repaired=len(union['repaired_domains']); leaves=sum(len(x['leaves']) for x in union['repaired_domains'])
    checks['canonical_union_shape']=direct==99 and repaired==7 and leaves==18
    # Full exact phase claims only when all audit/replay conditions pass.
    checks['full_phase_exact_claim_guarded']=(claims['FULL_PHASE_EXHAUSTION']['count']==60 and phase['geometry_and_mutation_audit']['all_pass'] and phase['formula_level_replay']['all_pass'] and phase['full_root_collar_separation']['all_60_collars_pairwise_disjoint'])
    out={'task':'R7N-058','checks':checks,'all_pass':all(checks.values()),'target_P_residual_fraction':residual,'note':'Independent campaign-level guard checks. Phase mutation mechanics are delegated to the independent phase audit; this layer tests cross-claim promotion/firewall rules.'}
    (ROOT/'results/R7N-058_mutation_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    if not out['all_pass']: raise SystemExit(1)
if __name__=='__main__': main()
