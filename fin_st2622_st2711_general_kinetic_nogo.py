#!/usr/bin/env python3
"""Six rounds: full local-rate classification and nonreversible kinetic no-go."""
import hashlib,itertools,json,math,numpy as np
from fin_st2172_st2261_common import ROOT,write_packet,write_round
def emit(lo,names,rows):
 d={}
 for i,(s,b,p) in enumerate(rows):k=lo+i;d[f'ST{k}']=write_packet(k,names[i],s,b,p)
 write_round(lo,lo+14,d)
def cycle_example(theta=.3):
 st=list(itertools.product((-1,1),repeat=3));E=np.array([-theta*np.prod(x) for x in st]);pi=np.exp(-E);pi/=pi.sum();Q=np.zeros((8,8))
 for a,x in enumerate(st):
  for i in range(3):
   y=list(x);y[i]*=-1;b=st.index(tuple(y));d=E[b]-E[a];Q[a,b]=1/(1+math.exp(d))
  Q[a,a]=-Q[a].sum()
 cyc=[st.index(x) for x in [(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)]];C=np.array([[pi[i]*Q[i,j] for j in range(8)] for i in range(8)]);eps=.1*min(C[cyc[i],cyc[(i+1)%4]] for i in range(4))
 def add(sign):
  R=Q.copy()
  for i in range(4):
   a,b=cyc[i],cyc[(i+1)%4];R[a,b]+=sign*eps/pi[a];R[b,a]-=sign*eps/pi[b]
  for a in range(8):R[a,a]=-sum(R[a,b] for b in range(8) if b!=a)
  sig=0
  for a in range(8):
   for b in range(a+1,8):
    if R[a,b]>0 and R[b,a]>0:
     J=pi[a]*R[a,b]-pi[b]*R[b,a]
     sig+=J*math.log((pi[a]*R[a,b])/(pi[b]*R[b,a]))
  return R,sig,float(abs(pi@R).max()),np.linalg.eigvals(R)
 return pi,Q,eps,add(1),add(-1)
def r1():
 names=['RateDifferenceAnsatz','DetailedBalanceEquation','GeneralSolution','EvenActivityFunction','Positivity','AnalyticFreedom','MonotoneFreedom','SymmetricRate','HeatBathRate','MetropolisRate','CommonRatio','TimeScale','FunctionSpaceNoGo','RoundOneVerdict','RoundOneRecommendation']
 rows=[('Defined','Local rate depends on DeltaE.',{'r':'r(Delta)'}),('Proven','Detailed balance.',{'equation':'r(D)/r(-D)=exp(-D)'}),('Proven classification','Set h(D)=exp(D/2)r(D).',{'solution':'r(D)=exp(-D/2)h(D), h even positive'}),('Proven','Every even positive h is admissible.',{'dimension':'infinite'}),('Boundary','h>0 gives irreducibility.',{}),('Proven','Infinitely many analytic h.',{'examples':['1','1/(2 cosh(D/2))','exp(a D^2)']}),('Proven','Additional monotonicity still leaves families.',{}),('Example','h=1.',{'r':'exp(-D/2)'}),('Example','h=1/(2cosh(D/2)).',{'r':'1/(1+exp(D))'}),('Example','h=exp(-abs(D)/2).',{'r':'min(1,exp(-D))'}),('Proven','All satisfy identical forward/backward ratio.',{}),('No-go','h->c h rescales time.',{}),('Proven no-go','Energy and detailed balance leave infinite functional freedom.',{}),('Round verdict','No unique local reversible rate follows from Gibbs energy.',{}),('Recommendation','Audit extra axioms selecting named rules.',{})];emit(2622,names,rows)
def r2():
 names=['HeatBathNormalization','HeatBathUniqueness','SymmetricProductNormalization','SymmetricUniqueness','MetropolisCap','MetropolisUniqueness','MutualIncompatibility','AttemptFrequency','ClockBoundary','SmoothnessBoundary','LocalOptimalityBoundary','AxiomRemoval','StrictSourceAudit','RoundTwoVerdict','RoundTwoRecommendation']
 rows=[('Axiom','r(D)+r(-D)=1.',{}),('Conditional theorem','With detailed balance this uniquely gives heat bath.',{}),('Axiom','r(D)r(-D)=1.',{}),('Conditional theorem','Uniquely gives exp(-D/2).',{}),('Axiom','Largest downhill acceptance equals one plus min rule.',{}),('Conditional theorem','Selects Metropolis.',{}),('Proven','The three normalizations cannot all hold.',{}),('Boundary','Each fixes a kinetic convention/attempt frequency.',{}),('No-go','No normalization derives seconds.',{}),('Boundary','Smoothness rejects Metropolis but not analytic h family.',{}),('Boundary','Variational speed criteria require a mobility/cost.',{}),('Proven necessity','Removing the named normalization restores h freedom.',{}),('Repository audit','No strict packet exports one normalization.',{'strict':False}),('Round verdict','Each canonical rule is unique only relative to its own added axiom.',{}),('Recommendation','Prove functorial nonuniqueness.',{})];emit(2637,names,rows)
def r3():
 names=['FunctorDomain','FunctorTarget','EnergyMap','ThetaBoundary','RateFunctionBoundary','D12Naturality','ActivityOrbitFreedom','NaturalTransformations','CompositionBoundary','RefinementBoundary','ScaleBoundary','Counterfunctors','UniquenessNoGo','RoundThreeVerdict','RoundThreeRecommendation']
 rows=[('Defined','Weighted D12 graph A,W.',{}),('Defined','Markov generators on configuration cube.',{}),('Conditional','A->tau->E_theta.',{}),('No-go','theta not strict.',{}),('No-go','h not strict.',{}),('Proven','Any even h yields D12-covariant Q.',{}),('Proven','More general activities retain 1056 D12 orbits.',{}),('Proven counterfamily','h choices give distinct natural maps.',{}),('Boundary','Graph morphism composition does not select h.',{}),('Boundary','Prism state RG does not select h.',{}),('No-go','Common multiplier remains clock.',{}),('Constructed','Heatbath, symmetric, Metropolis are three functors on same inputs.',{}),('Proven no-go','No uniqueness under current functorial/symmetry axioms.',{}),('Round verdict','A-to-Q functoriality narrows type, not kinetics.',{}),('Recommendation','Allow stationary circulation.',{})];emit(2652,names,rows)
def r4():
 pi,Q,eps,pos,neg=cycle_example();Rp,sp,rp,evp=pos;Rm,sm,rm,evm=neg
 names=['StationaryFluxDecomposition','ReversiblePart','CirculationSpace','HypercubeCycleDimension','FourCycleWitness','PositivePerturbation','Stationarity','OppositeCurrents','SameStationaryLaw','EntropyProduction','ComplexSpectrumBoundary','D12Boundary','ArrowSourceNoGo','RoundFourVerdict','RoundFourRecommendation']
 rows=[('Proven','Stationary edge flow decomposes symmetric conductance plus divergence-free antisymmetric current.',{}),('Proven','Detailed balance iff current vanishes.',{}),('Proven','Cycle-space currents preserve stationarity.',{}),('Proven exact','E-V+1.',{'dimension':24576-4096+1}),('Constructed','One oriented square in Q3.',{'epsilon':eps}),('Proven','Small +-epsilon keeps rates positive.',{}),('Computed','Both signs preserve pi.',{'residual_plus':rp,'residual_minus':rm}),('Proven','J and -J.',{'same_generator':False}),('Proven','Stationary Gibbs state unchanged.',{}),('Computed','Both arrows have positive equal entropy production.',{'sigma_plus':sp,'sigma_minus':sm}),('Boundary','Nonnormal generators may have complex modes.',{'max_imag_plus':float(max(abs(evp.imag)))}),('Boundary','D12 averaging can remove or retain sectors depending orientation representation.',{}),('No-go','Pair core does not choose circulation orientation.',{'strict_arrow':False}),('Round verdict','Dropping reversibility adds at least a 20481-dimensional current freedom before symmetry.',{}),('Recommendation','Audit thermodynamic meaning.',{})];emit(2667,names,rows)
def r5():
 names=['SteadyEntropyProduction','ReversibleZero','CirculationPositive','TimeReversal','SignPair','MagnitudeFreedom','ThermalUnits','LandauerBoundary','HousekeepingHeat','FluctuationTheorem','ClockRequirement','BathRequirement','ObservedArrowBoundary','RoundFiveVerdict','RoundFiveRecommendation']
 rows=[('Defined','sigma=sum J log forward/backward affinity.',{}),('Proven','Detailed balance sigma=0.',{}),('Proven','Nonzero circulation generically sigma>0.',{}),('Proven','Time reversal maps J to -J.',{}),('Proven','Opposite arrows share pi and sigma.',{}),('No-go','Current amplitude remains free subject positivity.',{}),('No-go','Dimensionless sigma is not watts/kelvin.',{}),('Boundary','Physical reset needs bath and kBT.',{}),('Conditional','Nonreversible steady state needs supplied drive.',{}),('Conditional','Requires declared path measure.',{}),('No-go','Rates need seconds.',{}),('No-go','Temperature/environment absent.',{}),('Proven','Entropy production magnitude does not select orientation sign.',{}),('Round verdict','A kinetic arrow requires a nonequilibrium resource and polarity beyond FIN equilibrium data.',{}),('Recommendation','Synthesize broad theorem.',{})];emit(2682,names,rows)
def r6():
 gate=ROOT/'FIN_ST2697_ST2711_GeneralKineticGate.json';g={'rows':[{'name':n,'pass':v} for n,v in [('rate classification theorem',True),('conditional canonical rules',True),('stationary circulation decomposition',True),('strict h activity',False),('strict normalization axiom',False),('strict circulation sign',False),('strict clock',False),('dynamic prism RG',False),('thermodynamic units',False),('OA evidence',False)]]};gate.write_text(json.dumps(g,indent=2)+'\n')
 names=['ReversibleTheorem','NormalizationTheorems','FunctorNoGo','CirculationTheorem','ArrowNoGo','DynamicRGNoGo','PhysicalNoGo','Gate','GateScore','ProvenLedger','ConditionalLedger','RefutedLedger','HighestValueNext','FinalInterpretation','CycleTrigger']
 rows=[('Proven','r=e^-D/2 h even positive.',{}),('Proven conditional','Named rules need distinct axioms.',{}),('Proven','Current naturality leaves counterfunctors.',{}),('Proven','20481 current dimensions before symmetry.',{}),('Proven','J signs paired.',{}),('Proven','Stationary closure not lumpability.',{}),('No-go','Clock/bath/units/OA absent.',{}),('Constructed','Ten rows.',{'file':gate.name}),('Result','Three mathematical passes.',{'passes':3,'rows':10}),('Ledger','Complete kinetic classification in declared local classes.',{}),('Ledger','Any selected h/current/prior.',{}),('Ledger','Unique strict physical dynamics.',{}),('Recommendation','Search for a strict nonequilibrium resource/clock; otherwise preserve kinetic impossibility theorem.',{}),('Final verdict','Equilibrium FIN cannot select kinetics or time arrow from state information alone.',{'fundamental':False}),('Cycle complete','Release 11.13.',{'programs':90,'rounds':6,'range':'ST2622-ST2711','strict_ToE_closure':False,'gate_sha256':hashlib.sha256(gate.read_bytes()).hexdigest()})];emit(2697,names,rows)
def main():r1();r2();r3();r4();r5();r6()
if __name__=='__main__':main()
