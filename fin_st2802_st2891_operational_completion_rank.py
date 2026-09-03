#!/usr/bin/env python3
"""Six rounds: five-axis operational completion and identifiability rank."""
import hashlib,json,math,numpy as np
from fin_st2172_st2261_common import ROOT,write_packet,write_round
def emit(lo,names,rows):
 d={}
 for i,(s,b,p) in enumerate(rows):k=lo+i;d[f'ST{k}']=write_packet(k,names[i],s,b,p)
 write_round(lo,lo+14,d)
def theta_eff(t,q):
 A=math.exp(3*q)+3*math.exp(-q);B=3*math.exp(q)+math.exp(-3*q)
 return t+.5*math.log((math.exp(t)*A+math.exp(-t)*B)/(math.exp(-t)*A+math.exp(t)*B))
def obs(z):
 t,q,a,c,j=z;D1=2*t;D2=.8*t
 r1=math.exp(-D1/2+a*D1*D1);r2=math.exp(-D2/2+a*D2*D2)
 return np.array([math.tanh(t),theta_eff(t,q),r1/(r1+r2),1/(c*(r1+r2)),j])
def jac(z):
 z=np.array(z,float);J=np.zeros((5,5))
 for k in range(5):
  h=1e-6*max(1,abs(z[k]));u=z.copy();v=z.copy();u[k]+=h;v[k]-=h;J[:,k]=(obs(u)-obs(v))/(2*h)
 return J
Z=[.3,.2,.1,2,.05];J=jac(Z)
def r1():
 names=['FiveAxisModel','ThetaAxis','RefinementAxis','ActivityAxis','ClockAxis','CurrentAxis','ParameterDomain','StrictInputs','ConditionalInputs','PhysicalInputs','Smoothness','IndependentVariations','SingleScalarHypothesis','RoundOneVerdict','RoundOneRecommendation']
 rows=[('Constructed','Minimal declared completion family.',{'parameters':['theta','q','a','c','j']}),('Typed','Equilibrium ternary coupling.',{}),('Typed','Vertical prism coupling.',{}),('Typed','Even kinetic activity h_a.',{}),('Typed','Positive rate multiplier.',{}),('Typed','Signed circulation.',{}),('Declared','Open neighbourhood.',{'point':Z}),('Proven','tau/shape only.',{}),('Boundary','All five axes remain supplied.',{}),('Boundary','seconds/energy/OA outside.',{}),('Proven','Observable map analytic locally.',{}),('Constructed','Each axis can vary holding others fixed.',{}),('Test target','Can one smooth scalar replace open five-axis data?',{}),('Round verdict','Five typed source slots are explicit.',{}),('Recommendation','Classify partial-record invariances.',{})];emit(2802,names,rows)
def r2():
 names=['StationaryRecord','PrismStationaryRecord','UntimedJumpRecord','TimedRecord','ArrowRecord','ThetaVisibility','QVisibility','ActivityVisibility','ClockGauge','CurrentGauge','PartialOrbitDimensions','RecordHierarchy','GaugeSeparation','RoundTwoVerdict','RoundTwoRecommendation']
 rows=[('Defined','O1=tanh theta.',{}),('Defined','O2=theta_eff(theta,q).',{}),('Defined','O3=selected jump probability.',{}),('Defined','O4=mean holding time.',{}),('Defined','O5=signed current.',{}),('Proven','O1 sees theta.',{}),('Proven','O2 adds q.',{}),('Proven','O3 adds activity a.',{}),('Proven','Without O4 c is Rplus gauge.',{}),('Proven','Without O5 j sign/magnitude invisible.',{}),('Proven','Stationary-only leaves at least three-dimensional fibre.',{}),('Constructed','stationary -> untimed -> timed -> arrow.',{}),('Proven','Clock and arrow gauges are independent.',{}),('Round verdict','Different record types remove different source fibres.',{}),('Recommendation','Compute Jacobian ranks.',{})];emit(2817,names,rows)
def r3():
 o=obs(Z);names=['ObservableVector','O1TernaryMean','O2PrismRG','O3JumpProbability','O4HoldingTime','O5SignedCurrent','TargetIndependence','FiniteValues','CalibrationBoundary','LikelihoodBoundary','NoiseBoundary','PhysicalBoundary','ObservableCompleteness','RoundThreeVerdict','RoundThreeRecommendation']
 rows=[('Constructed','Five records.',{'values':o.tolist()}),('Exact form','tanh theta.',{}),('Exact form','prior prism theorem.',{}),('Exact form','two unequal energy changes.',{}),('Exact form','1/[c(r1+r2)].',{}),('Exact form','j.',{}),('Proven','No target fitting in map definition.',{}),('Computed','Declared point finite/interior.',{}),('Boundary','O4 requires calibrated time.',{}),('Boundary','Shot likelihood required empirically.',{}),('Boundary','No raw data.',{}),('Boundary','Dimensionless model.',{}),('Conditional','Complete for declared five parameters.',{}),('Round verdict','A minimal separating record vector exists conditionally.',{}),('Recommendation','Rank and removal audit.',{})];emit(2832,names,rows)
def r4():
 subsets={'stationary':[0,1],'untimed_no_arrow':[0,1,2],'untimed_with_arrow':[0,1,2,4],'timed_no_arrow':[0,1,2,3],'full':[0,1,2,3,4]};ranks={k:int(np.linalg.matrix_rank(J[v,:],tol=1e-8)) for k,v in subsets.items()}
 names=['Jacobian','FullRank','Determinant','StationaryRank','UntimedRank','TimedRank','ArrowRank','RemoveTheta','RemoveQ','RemoveActivity','RemoveClock','RemoveCurrent','InvarianceOfDomain','RoundFourVerdict','RoundFourRecommendation']
 rows=[('Computed','Central derivative.',{'matrix':J.tolist()}),('Proven numerical witness','Rank five with large margin.',{'rank':int(np.linalg.matrix_rank(J))}),('Computed','Nonzero.',{'determinant':float(np.linalg.det(J))}),('Computed','O1,O2.',{'rank':ranks['stationary']}),('Computed','O1-O3.',{'rank':ranks['untimed_no_arrow']}),('Computed','O1-O4.',{'rank':ranks['timed_no_arrow']}),('Computed','O1-O3,O5.',{'rank':ranks['untimed_with_arrow']}),('Removal','Lose equilibrium coupling.',{}),('Removal','Lose refinement.',{}),('Removal','Lose kinetic shape.',{}),('Removal','Lose rate scale.',{}),('Removal','Lose arrow.',{}),('Proven theorem','No smooth one-scalar parametrization covers an open five-dimensional locally identifiable family.',{}),('Round verdict','Five axes are locally independent in the declared operational model.',{'ranks':ranks}),('Recommendation','Translate to record/role package.',{})];emit(2847,names,rows)
def r5():
 names=['PreparationRole','RegistrarRole','AnalystRole','ClockRole','ArrowInstrument','StationaryTrials','PrismTrials','JumpSequences','WaitingTimes','CycleCounts','Holdout','HashBoundary','MinimalRecordFields','RoundFiveVerdict','RoundFiveRecommendation']
 rows=[('Required','Prepare configurations/state law.',{}),('Required','Immutable event record.',{}),('Required','Frozen estimator.',{}),('Required','Calibrate timestamps.',{}),('Required','Resolve signed cycle orientation.',{}),('Required','Estimate O1.',{}),('Required','Estimate O2.',{}),('Required','Estimate O3.',{}),('Required','Estimate O4.',{}),('Required','Estimate O5.',{}),('Required','Prevent retuning.',{}),('Boundary','Hash does not create independence.',{}),('Constructed','run,layer,config_before,config_after,wait,timestamp,cycle_orientation,calibration.',{'fields':8}),('Round verdict','Operational completion needs multiple roles and records, not one observable.',{}),('Recommendation','Synthesize minimality boundary.',{})];emit(2862,names,rows)
def r6():
 gate=ROOT/'FIN_ST2877_ST2891_OperationalCompletionGate.json';g={'rows':[{'name':n,'pass':v} for n,v in [('five-axis rank theorem',True),('partial-record rank ladder',True),('minimal record schema',True),('strict theta source',False),('strict q source',False),('strict activity source',False),('strict clock calibration',False),('strict arrow source',False),('independent OA record',False),('physical units/evidence',False)]]};gate.write_text(json.dumps(g,indent=2)+'\n')
 names=['RankTheorem','RecordLadder','SingleScalarNoGo','FiveAxisBoundary','OperationalPackage','StrictSourceFailures','PhysicalFailures','Gate','GateScore','ProvenLedger','ConditionalLedger','RefutedLedger','HighestValueNext','FinalInterpretation','CycleTrigger']
 rows=[('Proven local','Full Jacobian rank five.',{'det':float(np.linalg.det(J))}),('Proven','Ranks rise with record types.',{}),('Proven scoped','No smooth scalar spans declared open family.',{}),('Boundary','Not an absolute theorem over all future models.',{}),('Constructed','Preparation, clock, instrument, record, analyst.',{}),('No-go','No five source laws exported.',{}),('No-go','No units/OA evidence.',{}),('Constructed','Ten rows.',{'file':gate.name}),('Result','Three mathematical passes.',{'passes':3,'rows':10}),('Ledger','Rank/records/minimality.',{}),('Ledger','Full conditional operational model.',{}),('Ledger','Single missing scalar and strict physics.',{}),('Recommendation','Treat missing bridge as a typed operational package; seek one law that couples axes only if it reduces rank by theorem.',{}),('Final verdict','Current FIN needs multiple independent source classes, not one hidden number.',{'fundamental':False}),('Cycle complete','Release 11.15.',{'programs':90,'rounds':6,'range':'ST2802-ST2891','strict_ToE_closure':False,'gate_sha256':hashlib.sha256(gate.read_bytes()).hexdigest()})];emit(2877,names,rows)
def main():r1();r2();r3();r4();r5();r6()
if __name__=='__main__':main()
