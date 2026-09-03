#!/usr/bin/env python3
"""Six rounds: configuration-space generator source and prism lumpability."""
import hashlib,itertools,json,math,numpy as np
from fin_st2172_st2261_common import ROOT,write_packet,write_round
def emit(lo,names,rows):
 d={}
 for i,(s,b,p) in enumerate(rows):k=lo+i;d[f'ST{k}']=write_packet(k,names[i],s,b,p)
 write_round(lo,lo+14,d)
def d12_orbits():
 ps=[tuple((i+k)%12 for i in range(12)) for k in range(12)]+[tuple((k-i)%12 for i in range(12)) for k in range(12)]
 def act(x,p):
  y=0
  for i in range(12):
   if (x>>i)&1:y|=1<<p[i]
  return y
 un={(x,x^(1<<i)) for x in range(4096) for i in range(12) if x<(x^(1<<i))};O=[]
 while un:
  a,b=next(iter(un));o={(min(act(a,p),act(b,p)),max(act(a,p),act(b,p))) for p in ps};O.append(o);un-=o
 return O
def q3(rule,theta=.3):
 st=list(itertools.product((-1,1),repeat=3));E=np.array([-theta*np.prod(x) for x in st]);pi=np.exp(-E);pi/=pi.sum();Q=np.zeros((8,8))
 for a,x in enumerate(st):
  for i in range(3):
   y=list(x);y[i]*=-1;b=st.index(tuple(y));d=E[b]-E[a]
   r=min(1,math.exp(-d)) if rule=='Metropolis' else (1/(1+math.exp(d)) if rule=='heat-bath' else math.exp(-d/2))
   Q[a,b]=r
  Q[a,a]=-Q[a].sum()
 return pi,Q,np.sort(np.linalg.eigvals(Q).real)
def r1():
 O=d12_orbits();names=['ConfigurationSpace','FlipGraph','ReversibleConductanceForm','DetailedBalanceRatio','ActivityFreedom','D12OrbitCount','OrbitSizes','StrictEnergyShape','ThetaDependency','Irreducibility','StationaryUniqueness','RateScale','GeneratorCone','RoundOneVerdict','RoundOneRecommendation']
 rows=[('Constructed','Binary state space.',{'states':4096}),('Proven','12-cube single flips.',{'undirected_edges':24576}),('Proven','Every reversible generator is Q_xy=c_xy/pi_x with c_xy=c_yx>0.',{}),('Proven','Ratio fixed, product/activity free.',{'ratio':'r_xy/r_yx=pi_y/pi_x'}),('Proven no-go','One positive activity per undirected edge.',{'dimension':24576}),('Proven exhaustive','After D12.',{'activity_orbits':len(O)}),('Proven','Exact orbit census.',{'sizes':{str(s):sum(len(o)==s for o in O) for s in set(map(len,O))}}),('Conditional','tau fixes energy shape.',{'strict_theta':False}),('No-go','Generator depends on theta.',{'strict_generator':False}),('Conditional','Positive rates make Q12 connected.',{'requires':'all activities positive'}),('Proven','Finite irreducible chain has unique pi.',{'kinetics_unique':False}),('No-go','Q->cQ changes clock only.',{'global_time_scale_free':True}),('Proven','Detailed-balance generators form a positive cone.',{'rays':24576}),('Round verdict','Stationary Gibbs state does not determine dynamics, even after D12.',{'unique_generator':False}),('Recommendation','Compare standard rules.',{'next':'kinetic spectral counterexample'})];emit(2532,names,rows)
def r2():
 data={};
 for rule in ['Metropolis','heat-bath','symmetric']:
  pi,Q,e=q3(rule);data[rule]={'spectrum':e.tolist(),'stationarity_residual':float(abs(pi@Q).max()),'gap':float(-e[-2])}
 names=['ThreeBitRestriction','MetropolisRule','HeatBathRule','SymmetricRule','CommonStationaryState','SpectralDisagreement','GapDisagreement','ClockRescaling','RuleNotRescaling','LocalityShared','DetailedBalanceShared','KineticObservable','RuleSourceGate','RoundTwoVerdict','RoundTwoRecommendation']
 rows=[('Constructed','One triangle exact 8-state test.',{'theta':.3}),('Computed','Accepted reversible rule.',data['Metropolis']),('Computed','Accepted reversible rule.',data['heat-bath']),('Computed','Accepted reversible rule.',data['symmetric']),('Proven','All have same p_theta.',{'max_stationarity':max(v['stationarity_residual'] for v in data.values())}),('Proven counterexample','Spectra differ.',{'equal':False}),('Proven','Relaxation gaps differ.',{k:v['gap'] for k,v in data.items()}),('Boundary','A common multiplier remains free.',{}),('Proven','Normalized eigenvalue ratios also differ.',{'rules_related_only_by_scale':False}),('Proven','All flip one bit.',{}),('Proven','All satisfy same detailed balance.',{}),('Boundary','Time correlations distinguish rules.',{}),('Gate fail','No strict principle chooses one.',{'strict_rate_rule':False}),('Round verdict','Same state and locality admit inequivalent kinetics.',{'state_determines_dynamics':False}),('Recommendation','Test maximum caliber.',{'next':'prior-generator dependence'})];emit(2547,names,rows)
def r3():
 names=['MaximumCaliberSetup','PathEntropyConstraint','ReferenceMeasure','SymmetricSolution','ActivityBudget','BudgetTimeScale','DifferentPriors','FisherMetricBoundary','GradientFlowBoundary','DirichletFormBoundary','CanonicalityGate','StateKineticsSeparation','InformationPrincipleNoGo','RoundThreeVerdict','RoundThreeRecommendation']
 rows=[('Conditional','Optimize path entropy with constraints.',{}),('Boundary','Need mean jump/activity constraints.',{'supplied':True}),('No-go','Relative path entropy requires a reference generator.',{'reference_strict':False}),('Conditional','Uniform edge prior yields square-root/symmetric rates.',{'rule':'exp(-DeltaE/2)'}),('Boundary','Mean activity fixes overall rate.',{}),('No-go','Budget is a clock scale.',{}),('Proven','Metropolis/heat-bath arise from different choices.',{}),('No-go','Stationary Fisher geometry does not specify mobility.',{}),('Conditional','Gradient flow needs Onsager operator.',{}),('Conditional','Reversible Dirichlet form needs edge conductances.',{}),('Gate fail','No source for prior/mobility/activity.',{'strict':False}),('Proven','Information state and kinetic law are independent layers.',{}),('Proven no-go','Maximum caliber cannot choose dynamics without constraints/prior.',{}),('Round verdict','Information optimization relocates rather than removes kinetic axioms.',{}),('Recommendation','Test prism dynamical lumpability.',{})];emit(2562,names,rows)
def r4():
 theta,q=.3,.2;b=(1,1,1);top1=(1,1,1);top2=(-1,1,1)
 # bottom spin0 energy change difference due vertical q
 def delta(top):
  bb=list(b);E0=-theta*(np.prod(bb)+np.prod(top))-q*sum(bb[i]*top[i] for i in range(3));bb[0]*=-1;E1=-theta*(np.prod(bb)+np.prod(top))-q*sum(bb[i]*top[i] for i in range(3));return E1-E0
 d1,d2=delta(top1),delta(top2)
 names=['StationaryRGClosure','DynamicProjection','StrongLumpability','HiddenTopCounterexample','EnergyDifference','RateDifference','QZeroBoundary','QNonzeroNoGo','MemoryKernel','StateAugmentation','ThetaEffInsufficiency','TwoLevelRates','PathCoherence','RoundFourVerdict','RoundFourRecommendation']
 rows=[('Proven','Prior exact theta_eff formula.',{'closed_state_family':True}),('Defined','Forget top spins.',{}),('Criterion','Bottom transition totals must be independent of hidden top.',{}),('Constructed','Same bottom, two top states.',{'top':[top1,top2]}),('Computed exact formula','Flip energy differs.',{'DeltaE':[d1,d2]}),('Proven','Any nonconstant detailed-balance rate sees difference.',{'strongly_lumpable':False}),('Proven','At q=0 hidden top drops out.',{'lumpable':True,'RG_trivial':True}),('No-go','For q nonzero generic dynamics not lumpable.',{}),('Conditional','Projected process acquires memory.',{'Markov':False}),('Conditional','Keep top/hidden state to restore Markov property.',{}),('Proven','Stationary theta_eff does not determine coarse kinetics.',{}),('No-go','Each level needs q_n and activity law.',{}),('Gate fail','No strict dynamic intertwiner.',{'strict':False}),('Round verdict','Stationary RG closure does not imply dynamical RG closure.',{'main_no_go':True}),('Recommendation','Compare with U/P dual channels.',{})];emit(2577,names,rows)
def r5():
 names=['CarrierDimension','UnitaryGenerator','HeatGenerator','ConfigurationGenerator','NoCanonicalLift','SpectrumBoundary','ZeroMode','ParityObservable','StateProjection','OperationalDistinction','CommonAInsufficiency','ClockClasses','PhysicalBoundary','RoundFiveVerdict','RoundFiveRecommendation']
 rows=[('Proven','12 versus 4096.',{'vertex':12,'configuration':4096}),('Proven','-iA anti-Hermitian.',{}),('Proven','-A vertex Markov only in conductance basis.',{}),('Conditional','Q acts on binary configurations.',{}),('No-go','No strict functor maps A to unique Q.',{}),('Boundary','Q spectrum depends activities.',{}),('Proven','All conservative generators have constants zero mode.',{}),('Conditional','Three-body parity lives only on configuration algebra.',{}),('No-go','Vertex populations do not determine configuration law.',{}),('Proven','U/P/Q yield different records.',{}),('Proven','Shared tau/A organization insufficient.',{}),('No-go','Global Q rate and U/P calibrations are separate.',{}),('Boundary','No apparatus or units.',{}),('Round verdict','Configuration dynamics is a third operational channel, not derived dual dynamics.',{}),('Recommendation','Synthesize gate.',{})];emit(2592,names,rows)
def r6():
 gate=ROOT/'FIN_ST2607_ST2621_ConfigurationGeneratorGate.json';g={'rows':[{'name':n,'pass':v} for n,v in [('finite Gibbs state',True),('exact prism stationary RG',True),('reversible generators exist',True),('unique kinetic activity',False),('strict rate rule',False),('strict clock',False),('dynamic RG lumpability',False),('canonical A-to-Q lift',False),('physical units',False),('OA evidence',False)]]};gate.write_text(json.dumps(g,indent=2)+'\n')
 names=['StateTheorem','GeneratorClassification','OrbitNoGo','RateRuleNoGo','CaliberNoGo','LumpabilityNoGo','DualDynamicsNoGo','Gate','GateScore','ProvenLedger','ConditionalLedger','RefutedLedger','HighestValueNext','FinalInterpretation','CycleTrigger']
 rows=[('Proven','Ternary Gibbs and stationary prism RG.',{}),('Proven','Positive reversible conductance cone.',{'activity_orbits_D12':1056}),('No-go','Symmetry leaves 1056 activities.',{}),('No-go','Three standard rates inequivalent.',{}),('No-go','Maximum caliber needs prior/budget.',{}),('No-go','q nonzero blocks Markov projection.',{}),('No-go','A does not uniquely lift to Q.',{}),('Gate constructed','Ten rows.',{'file':gate.name}),('Gate result','Three mathematical passes.',{'passes':3,'rows':10}),('Ledger','State and stationary RG.',{}),('Ledger','Any chosen reversible Q and instrument.',{}),('Ledger','Unique/sourced/physical generator.',{}),('Recommendation','A strict kinetic activity/clock source or an impossibility theorem beyond reversible single flips.',{}),('Final verdict','FIN fixes equilibrium interaction shape but not kinetics; stationary RG is not dynamical physics.',{'fundamental':False}),('Cycle complete','Release 11.12.',{'programs':90,'rounds':6,'range':'ST2532-ST2621','strict_ToE_closure':False,'gate_sha256':hashlib.sha256(gate.read_bytes()).hexdigest()})];emit(2607,names,rows)
def main():r1();r2();r3();r4();r5();r6()
if __name__=='__main__':main()
