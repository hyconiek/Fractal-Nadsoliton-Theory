from __future__ import annotations
from fractions import Fraction as F
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
# accepted historical exact rational certificates from the supplied discord report
DELTA_L=F(234218204114629,2000000000000000)
DELTA_SMALL=F(27234855667,12500000000000)
D0=F(911271283377980905464649,26400000000000000000000000000000)
# current strict spectral interval endpoints embedded in campaign certificates
L5_HI=F(22986062720791,10000000000000)
L6_LO=F(234218204114629,100000000000000)
CURRENT_SMALL=(L6_LO-L5_HI)/20

def run():
    factor=F(49,50)**2
    return {
      'R7P-113':{
        'scope_matrix':{
          'canonical_V':'fixed canonical two-body interaction; stationary-state statements with declared marginals',
          'mixed_family':'uniform discord floor for stationary states whose first marginal is the supplied C; existence for every family parameter is not asserted',
          'pure_flow_equivalent_family':'tradeoff uses average labelled marginals after swap twirl; individual marginal equality is stronger and not silently assumed',
          'separable_stationarity':'22-product rank-132 explicit separable stationary construction belongs to the separable-stationarity report, not a complete summary of later discord robustness'
        }},
      'R7P-114':{
        'density_gap_exact_expressions':['lambda6/20','(lambda6-lambda5)/20'],
        'historical_certificates':{'Delta_L':str(DELTA_L),'delta_L':str(DELTA_SMALL),'d0':str(D0)},
        'current_provider_gap_lower':str(CURRENT_SMALL),
        'historical_minus_current_gap_lower':str(DELTA_SMALL-CURRENT_SMALL),
        'replay_status':'REPLAY_PROVENANCE_GAP: current rounded strict interval is lower by exactly 5e-16; historical theorem remains an accepted source-package result but is not independently re-proved by this rounded provider alone.'},
      'R7P-115':{
        'shared':['strict spectral eigenvalues/gaps as numerical-algebraic inputs'],
        'missing_causal_maps':['classical p/theta landscape -> bipartite density operator','rank-seven gain law -> quantum stationary generator','discord-preparation resources -> localization selector','quantum access model -> active-gain source'],
        'verdict':'No causal discord-to-localization theorem follows from a shared gap factor.'},
      'R7P-116':{
        'status':'EXISTING_THEOREM_REPLAY_SCOPED',
        'facts':['explicit separable stationary state with both marginals C','conditional local preparation uses 22 oriented product terms and success probability 1/2','joint/correlation information is not fixed by local marginals alone'],
        'nonconclusion':'No claim that local FIN geometry identifies the joint preparation law.'},
      'R7P-117':{
        'loading_family':'C_gamma=I/12+gamma W, gamma in [49/1000,1/20]',
        'positivity_reason':'C_gamma=(1-alpha)I/12+alpha C_0 with alpha=20 gamma in [49/50,1]; convexity preserves positivity.',
        'gap_scaling':'both relevant simple spectral gaps scale linearly with gamma for gamma>0',
        'conditional_floor_factor':str(factor),
        'conditional_floor_lower':str(factor*D0),
        'status':'CONDITIONAL_ROBUSTNESS_THEOREM, conditional on the accepted historical delta_L/Delta_L certificates and the same canonical interaction proof hypotheses.'},
      'R7P-118':{
        'resources':['program copies C-prime','shared classical cut choice','local projective measurements','classical comparison/heralding','joint Swap or other joint measurement when joint distinguishability is discussed'],
        'separability_note':'Separable output does not by itself imply a universal unknown-input LOCC channel.'},
      'R7P-119':{'status':'OUT_OF_SCOPE_NO_NEW_LEGACY_COMPARISON_NEEDED','reason':'No new theorem in O depends on a legacy signed-weight program; strict constants are not transferred to legacy roles.'},
      'R7P-120':{'classification':'ONE_NEW_SCOPED_ROBUSTNESS_RESULT_PLUS_OPERATIONAL_NONIDENTIFIABILITY; NO_CAUSAL_SPECTRAL_BRIDGE','new_result':'R7P-117 loading interval','replayed_results':['discord floor and tradeoff','separable stationary construction'],'nonclosure':['no active gain source','no selector','no quantum-to-localization causal map']}
    }

def main():
    out=run(); (ROOT/'results/R7P-113_120_quantum_bridge.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps({'gap':out['R7P-114']['historical_minus_current_gap_lower'],'factor':out['R7P-117']['conditional_floor_factor']},indent=2))
if __name__=='__main__': main()
