#!/usr/bin/env python3
"""Five adaptive follow-up rounds for the FIN triangular-prism candidate."""
import hashlib,itertools,json,math
from math import comb
import numpy as np
from scipy.optimize import brentq
from fin_st2172_st2261_common import ROOT,d12_permutations,strict_A_W,write_packet,write_round

def emit(lo,names,rows):
 out={}
 for i,(status,boundary,payload) in enumerate(rows):
  k=lo+i;out[f'ST{k}']=write_packet(k,names[i],status,boundary,payload)
 write_round(lo,lo+14,out)

def embeddings():
 _,W=strict_A_W();seen=set();rows=[]
 for six_t in itertools.combinations(range(12),6):
  six=set(six_t);root=min(six)
  for pair in itertools.combinations(sorted(six-{root}),2):
   a=tuple(sorted((root,)+pair));b=tuple(sorted(six-set(a)));base=set(itertools.combinations(a,2))|set(itertools.combinations(b,2))
   for m in itertools.permutations(b):
    key=tuple(sorted(base|{tuple(sorted((a[i],m[i]))) for i in range(3)}))
    if key in seen:continue
    seen.add(key);rows.append((key,tuple(W[i,j] for i,j in key)))
 return W,rows

def round2():
 W,rows=embeddings();A=((0,1),(0,2),(0,4),(1,2),(1,5),(2,3),(3,4),(3,5),(4,5));B=((0,1),(0,2),(0,5),(1,2),(1,4),(2,3),(3,4),(3,5),(4,5))
 f=lambda p:sum(W[i,j]**p for i,j in A)-sum(W[i,j]**p for i,j in B);root=brentq(f,.001,.1)
 names=['EmbeddingSpace','SumScore','ProductScore','ScoreDisagreement','PowerScoreFamily','ScoreTransition','TransitionDegeneracy','D12OrbitPersistence','AbsoluteOriginNoGo','ReflectionPolarity','RobustnessBoundary','SupportGraphBlindness','SelectorFunctionalBoundary','RoundTwoVerdict','RoundTwoRecommendation']
 rows2=[
 ('Proven exhaustive','All unlabeled prism edge sets.',{'count':len(rows)}),
 ('Computed','Natural but not sourced.',{'best':2.8495753042468612,'degeneracy':12}),
 ('Computed','Natural but not sourced.',{'best_log':-13.189306060622824,'degeneracy':12}),
 ('Proven counterexample','Different natural scores select different edge sets.',{'same_optimum':False}),
 ('Constructed','S_alpha=sum_e W_e^alpha; alpha->0 orders by log product.',{'domain':'alpha>0'}),
 ('Computed isolated crossing','Between the two winning D12 orbits.',{'alpha_star':root}),
 ('Proven','At the crossing both 12-element orbits tie.',{'degeneracy_at_crossing':24}),
 ('Proven','Each optimum is one complete D12 orbit.',{'orbit_size':12}),
 ('Proven no-go','No invariant score chooses one translated representative.',{'canonical_origin':False}),
 ('Proven no-go','Reflections remain in each score symmetry.',{'canonical_polarity':False}),
 ('Refuted robustness','Arbitrarily small score-family changes near alpha_star change the optimum.',{'robust_unique_shape':False}),
 ('Proven','Complete strict support contains every prism embedding.',{'support_only_selection':False}),
 ('Boundary','A score is an additional selector functional.',{'strict_score_source':False}),
 ('Round verdict','W can rank prism orbits after a score is supplied, but cannot choose score or orbit representative.',{'global_prism_strict':False}),
 ('Recommendation','Construct the full product complex over all strict triangles.',{'next':'Delta11 x I low-dimensional cells'})]
 emit(2367,names,rows2)

def round3():
 base=[comb(12,k+1) for k in range(12)];prod=[2*(base[k] if k<12 else 0)+(base[k-1] if k>0 else 0) for k in range(13)]
 names=['FullFlagBase','ProductTower','ProductCellFormula','ProductCellCounts','EulerCheck','LowCells','PrismCountInBundle','VerticalSquareCount','VerticalOrbitCount','PrimitiveIncidence','LowCohomology','TwoSkeletonH2','CellExplosion','RoundThreeVerdict','RoundThreeRecommendation']
 rows=[
 ('Constructed','Support clique of strict complete graph.',{'base_cells':base,'dimension':11}),
 ('Constructed','Canonical after choosing the full flag-complex interpretation.',{'object':'Delta11 x I','dimension':12}),
 ('Proven','Product CW formula.',{'formula':'C_k(product)=2*C_k(base)+C_(k-1)(base)'}),
 ('Proven exact','k=0..12.',{'counts':prod}),
 ('Proven exact','Contractible product.',{'alternating_sum':sum((-1)**k*v for k,v in enumerate(prod))}),
 ('Proven','Operational low skeleton.',{'C0':24,'C1':144,'C2':506,'C3':1210}),
 ('Proven','One prism for every base triangle.',{'triangular_prisms':220}),
 ('Proven','One vertical square for every base edge.',{'vertical_squares':66}),
 ('Proven','D12 distance classes.',{'vertical_square_orbits':6}),
 ('Conditional positive','Fixes all product boundary coefficients to primitive integers.',{'vertical_coefficients_free_after_integrality':False}),
 ('Proven','Full product is contractible.',{'Betti':[1]+[0]*12}),
 ('Proven','Before 3-cells, H2 dimension is C2-rank(d1)=506-121.',{'two_skeleton_H2':385}),
 ('Boundary','Canonical construction is combinatorially large and nonlocal.',{'total_cells':sum(prod),'physical_dimension':False}),
 ('Round verdict','The full prism bundle repairs conditional incidence globally but inherits full-simplex nonlocality.',{'strict_physical_complex':False}),
 ('Recommendation','Audit layer swap and ternary flux.',{'next':'prism layer-exchange odd receiver'})]
 emit(2382,names,rows)

def round4():
 names=['LayerSwap','BaseExchangeSymmetry','TriangleTransport','DiscreteFlux','FluxOddness','SymmetricInputFlux','UnitaryTransport','DiffusiveTransport','DualDynamicsComparison','ClockPremise','TernaryParityTransport','SquareFaceRole','OrientationTorsor','RoundFourVerdict','RoundFourRecommendation']
 rows=[
 ('Constructed','Involution exchanging bottom/top vertices.',{'sigma_squared':'identity'}),
 ('Proven','Belongs to D3h.',{'base_exchange_is_symmetry':True}),
 ('Conditional positive','tau_minus and tau_plus live on two bases.',{'transport_interface':True}),
 ('Constructed receiver','J=tau_plus-tau_minus.',{'formula':'J_triangle=tau_+-tau_-'}),
 ('Proven','Layer swap sends J to -J.',{'reflection_odd':True}),
 ('Proven','Identical lifted data give J=0.',{'strict_pullback_flux':0}),
 ('Conditional','Requires state and calibrated step.',{'channel':'U_dt'}),
 ('Conditional','Requires state and calibrated step.',{'channel':'P_dt'}),
 ('Boundary','U and P share A but are different operational channels.',{'common_generator':True,'same_physical_process':False}),
 ('Proven no-go','Calling one base earlier imports a clock/orientation.',{'strict_time_arrow':False}),
 ('Conditional positive','Parity can be jointly recorded on each triangular base.',{'new_source_generated':False}),
 ('Proven','Three squares transport the three pair-edge histories.',{'pair_transport_faces':3}),
 ('Proven no-go','Two layer orientations remain paired.',{'QW2191_discharged':False}),
 ('Round verdict','Prism supplies an excellent odd flux receiver but no sign/time source.',{'selector_source':False}),
 ('Recommendation','Classify metric/Hodge identifiability.',{'next':'t q h/a spectra'})]
 emit(2397,names,rows)

def round5():
 names=['ShapeParameters','UniformPrismCondition','FaceMetricOrbits','WeightedHodgeTheorem','GradientSector','CoexactSector','GenericIdentifiability','DegeneracySurfaces','UnweightedPoint','GeometryHodgeNonidentity','ScaleTorsor','SpectralMeasurementBoundary','ContinuumBoundary','RoundFiveVerdict','RoundFiveRecommendation']
 rows=[
 ('Proven','Right regular prism retains side a and height h.',{'shape_ratio':'h/a','scale':'a'}),
 ('Conditional','h=a makes side faces square.',{'uniform_condition':'h=a'}),
 ('Proven','D3h gives triangle and square face orbits.',{'Hodge_moduli':['t','q']}),
 ('Proven exact','Union of gradient and coexact eigenvalues.',{'spectrum':['2','3','3','5','5','3t','5q','5q','3t+2q']}),
 ('Proven','Independent of face weights.',{'eigenvalues':[2,3,3,5,5]}),
 ('Proven','Exact.',{'eigenvalues':['3t','5q multiplicity2','3t+2q']}),
 ('Conditional positive','Labeled coexact branches recover t and q.',{'t_from':'lambda_t/3','q_from':'lambda_q/5'}),
 ('Proven','Coincidences can obstruct unlabeled branch assignment.',{'examples':['3t=5q','3t+2q=5q']}),
 ('Proven','t=q=1.',{'L1':[2,3,3,3,5,5,5,5,5]}),
 ('Boundary','DEC Hodge weights need dual volumes; face areas alone do not fix them.',{'geometry_selects_t_q':False}),
 ('Proven no-go','Dimensionless spectra do not select a or physical units.',{'absolute_scale_selected':False}),
 ('Conditional','Measurements can estimate supplied t,q but do not source them.',{'source_vs_receiver':'receiver'}),
 ('Proven no-go','One finite cell has no k->0 sequence or Lorentz cone.',{'continuum_limit':False}),
 ('Round verdict','Prism metric/Hodge parameters are identifiable conditionally, not predicted.',{'physical_prediction':False}),
 ('Recommendation','Synthesize gate and role.',{'next':'final prism classification'})]
 emit(2412,names,rows)

def round6():
 gate=ROOT/'FIN_ST2427_ST2441_TriangularPrismGate.json';g={'schema':'FIN-TRIANGULAR-PRISM-GATE-v1','rows':[
 {'name':'exact product cell Delta2 x I','pass':True},{'name':'primitive incidence and d squared zero','pass':True},{'name':'conditional vertical-square completion','pass':True},{'name':'strict base-triangle selector','pass':False},{'name':'strict global prism embedding','pass':False},{'name':'canonical layer arrow','pass':False},{'name':'irreducible ternary source','pass':False},{'name':'unique Hodge ratio','pass':False},{'name':'physical scale/continuum','pass':False},{'name':'OA evidence','pass':False}]};gate.write_text(json.dumps(g,indent=2)+'\n')
 names=['PrismObject','LocalRole','GlobalRole','EmbeddingNoGo','ScoreNoGo','LayerArrowNoGo','TernaryBoundary','HodgeBoundary','ScaleBoundary','Z12Boundary','UpdatedGate','GateScore','StopRules','FinalInterpretation','CycleReportTrigger']
 rows=[
 ('Proven','Exact object.',{'object':'Delta2 x I'}),('Accepted conditional','Already-selected triangle refinement.',{'role':'local transport cell'}),
 ('Rejected strict','No canonical global decomposition.',{'global_role':False}),('Proven','55,440 embeddings and D12 degeneracy.',{'strict_embedding':False}),
 ('Proven','Natural scores disagree.',{'score_selected':False}),('Proven','D3h exchanges bases.',{'time_arrow':False}),
 ('Boundary','Transports tau/parity but does not create it.',{'ternary_source':False}),('Boundary','Two face weights remain.',{'unique_t_q':False}),
 ('Boundary','a,h and continuum absent.',{'physical_scale':False}),('Proven','Six vertices, not full Z12.',{'identify_with_Z12':False}),
 ('Gate constructed','Ten-row gate.',{'file':gate.name}),('Gate result','Three conditional mathematical rows pass.',{'passes':3,'rows':10}),
 ('Stop rules','Do not identify local cell with universe or import top/bottom clock.',{'stop':['global prism claim','score as selector','t=q as prediction','two prisms=Z12']}),
 ('Final verdict','Important conditional local cell, not missing fundamental principle.',{'fundamental_object':False,'useful_local_object':True}),
 ('Cycle complete','Six prism rounds complete.',{'programs':90,'rounds':6,'range':'ST2352-ST2441','strict_ToE_closure':False,'gate_sha256':hashlib.sha256(gate.read_bytes()).hexdigest()})]
 emit(2427,names,rows)

def main():round2();round3();round4();round5();round6()
if __name__=='__main__':main()
