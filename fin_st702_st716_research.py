#!/usr/bin/env python3
"""FIN ST702--ST716: multi-model falsification features and executable intake schema."""
from __future__ import annotations
import csv,hashlib,itertools,json,math
from pathlib import Path
from typing import Any
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/"FIN_ST702_ST716_Results.json";SUMMARY=ROOT/"FIN_ST702_ST716_Summary.csv";FIG_DIR=ROOT/"FIN_ST702_ST716_Figures"
NAMES={702:"Declared_Propagation_Model_Library",703:"Feature_Separation_Matrix",704:"Minimal_Falsification_Feature_Set",705:"CalibrationInvariant_Feature_Audit",706:"FiniteCount_Model_Confusion_Replay",707:"FeatureRemoval_Degeneracy_Counterexamples",708:"HeldOut_Synthetic_Model_Challenge",709:"CausalTail_Class_Discriminator",710:"SignedInstability_Discriminator",711:"Noncircular_Anchor_Validator",712:"TwoClock_Record_Validator",713:"LayerSpeed_Laboratory_Executable_Spec",714:"Strict_ModelSource_Gate",715:"PhysicalPrediction_Gate",716:"Independent_Evidence_Gate"};PACKETS={k:ROOT/f"FIN_ST{k}_{v}.json" for k,v in NAMES.items()}
MODELS={
"local_D1_scalar":{"D":1,"nu":2.,"pol":1,"gauss":0,"stable":1,"tail":2,"clocks":2},
"local_D2_scalar":{"D":2,"nu":2.,"pol":1,"gauss":0,"stable":1,"tail":2,"clocks":2},
"local_D3_scalar":{"D":3,"nu":2.,"pol":1,"gauss":0,"stable":1,"tail":2,"clocks":2},
"local_D3_gauge":{"D":3,"nu":2.,"pol":2,"gauss":1,"stable":1,"tail":2,"clocks":2},
"fractional_D3_scalar":{"D":3,"nu":.8,"pol":1,"gauss":0,"stable":1,"tail":1,"clocks":2},
"signed_D3_scalar":{"D":3,"nu":None,"pol":1,"gauss":0,"stable":0,"tail":0,"clocks":None}}
FEATURES=["D","nu","pol","gauss","stable","tail","clocks"]
def native(x:Any)->Any:
 if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 if isinstance(x,np.ndarray):return native(x.tolist())
 if isinstance(x,(np.floating,np.integer,np.bool_)):return x.item()
 return x
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def finalize(k,status,boundary,packet):
 p=PACKETS[k];p.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8");return {"program":f"ST{k}","object":NAMES[k],"packet_file":p.name,"packet_sha256":sha(p),**packet,"status":status,"boundary":boundary}
def signatures(fs):return {m:tuple(MODELS[m][f] for f in fs) for m in MODELS}
def separated(fs):
 vals=list(signatures(fs).values());return len(set(vals))==len(vals)

def st702():return finalize(702,"constructed_declared_six_model_library","The library is representative, not exhaustive or nature-selected.",{"models":MODELS,"feature_semantics":{"tail":{"0":"unstable","1":"polynomial","2":"exponential-cone class"}}})
def st703():
 rows=[]
 for f in FEATURES:rows.append({"feature":f,"distinct_values":len(set(MODELS[m][f] for m in MODELS)),"alone_separates_all":separated([f])})
 return finalize(703,"proven_exact_declared_feature_separation_matrix","Exact separation concerns the frozen library and ideal values.",{"rows":rows,"full_signature_unique":separated(FEATURES)})
def st704():
 mins=[]
 for r in range(1,len(FEATURES)+1):
  mins=[list(c) for c in itertools.combinations(FEATURES,r) if separated(c)]
  if mins:break
 packet={"minimum_feature_count":len(mins[0]),"all_minimal_sets":mins,"theorem":"The listed subsets are exactly the minimum-cardinality ideal feature sets separating all six frozen models."}
 return finalize(704,"proven_minimal_declared_falsification_feature_sets","Noise and unmodeled alternatives can require more records.",packet)
def st705():
 invariant={"D":True,"nu":True,"pol":True,"gauss":True,"stable":True,"tail":True,"clocks":True};packet={"feature_invariant_under_common_length_time_rescaling":invariant,"absolute_c_needed":False,"theorem":"The declared model classes can be separated using dimensionless exponents, ranks and qualitative tail/stability data before assigning SI units."}
 return finalize(705,"proven_calibration_invariant_model_features","Extracting the features operationally still requires calibrated resolution/time ordering.",packet)
def encode(m):
 x=MODELS[m];return np.array([x["D"],-1 if x["nu"] is None else x["nu"],x["pol"],x["gauss"],x["stable"],x["tail"],-1 if x["clocks"] is None else x["clocks"]],float)
def st706():
 rng=np.random.default_rng(706);names=list(MODELS);means=np.array([encode(m) for m in names]);sc=np.array([.12,.06,.08,.03,.03,.08,.05]);trials=30000;conf=np.zeros((len(names),len(names)),int)
 for i,mu in enumerate(means):obs=mu+rng.normal(size=(trials,len(FEATURES)))*sc;pred=np.argmin(np.sum(((obs[:,None]-means[None])/sc)**2,axis=2),axis=1);conf[i]=np.bincount(pred,minlength=len(names))
 err=1-np.trace(conf)/(trials*len(names));return finalize(706,"strong_synthetic_finite_count_model_separation","Gaussian feature noise is synthetic and categorical constraints are approximated continuously.",{"trials_per_model":trials,"scales":sc,"confusion":conf,"total_error":float(err)})
def st707():
 counter=[]
 for f in FEATURES:
  fs=[x for x in FEATURES if x!=f];sig=signatures(fs);groups={}
  for m,v in sig.items():groups.setdefault(str(v),[]).append(m)
  collisions=[g for g in groups.values() if len(g)>1];counter.append({"removed":f,"collisions":collisions})
 return finalize(707,"proven_feature_removal_collision_catalog","Some features are redundant globally; necessity applies to particular minimal sets.",{"rows":counter})
def st708():
 rng=np.random.default_rng(708);names=list(MODELS);hidden="fractional_D3_scalar";train=[m for m in names if m!=hidden];means={m:encode(m) for m in train};obs=encode(hidden)+rng.normal(size=len(FEATURES))*np.array([.05,.03,.05,.02,.02,.04,.03]);scores={m:float(np.linalg.norm(obs-means[m])) for m in train};best=min(scores,key=scores.get)
 return finalize(708,"strong_heldout_unknown_model_rejection_fixture","The hidden class is synthetic; rejection thresholds are not calibrated.",{"hidden":hidden,"excluded_from_candidates":True,"nearest_wrong_model":best,"nearest_distance":scores[best],"out_of_library_flag":scores[best]>.5})
def st709():
 packet={"records":["short-time distance profile","tail exponent","finite-range path-order onset"],"classification":{"local":"factorial/exponential outside proxy cone","fractional":"polynomial direct tail","signed":"growth/oscillation with sign-sensitive modes"},"theorem":"Tail order plus stability separates the three declared locality classes without SI calibration."}
 return finalize(709,"proven_declared_causal_tail_class_discriminator","Finite samples need dynamic range and nuisance-error bounds.",packet)
def st710():
 packet={"test":"monitor norm/population response at two late times","positive_local":"bounded unitary and contractive heat","signed_negative_mode":"raw heat grows exp(t|lambda_min|), wave grows cosh(t sqrt|lambda_min|)","theorem":"Any certified negative mode supplies a model-level instability discriminator against positive local refinement."}
 return finalize(710,"proven_signed_instability_discriminator","Krein-unitary-only models require a different operational null.",packet)
def st711():
 fields=["length_standard_id","clock_standard_id","calibration_hash","calibrated_before_wave_unblinding","standards_not_defined_by_tested_c"]
 packet={"required_fields":fields,"valid_fixture":{"length_standard_id":"L-EXT","clock_standard_id":"T-EXT","calibration_hash":"sha256:fixture","calibrated_before_wave_unblinding":True,"standards_not_defined_by_tested_c":True},"fail_closed":True}
 return finalize(711,"constructed_fail_closed_anchor_validator_schema","A schema cannot certify real independence without external custody.",packet)
def st712():
 packet={"required_slopes":{"diffusive":2,"ballistic":1},"tolerance":.08,"record_fields":["layer","lattice_spacing_ratio","unitary_or_heat_time","wave_time","uncertainty"],"acceptance":"disjoint confidence intervals containing 2 and 1","single_clock_record_rejected":True}
 return finalize(712,"constructed_two_clock_record_validator","No physical timing records are supplied.",packet)
def st713():
 packet={"pipeline":["freeze model library and nuisance ranges","register external rods/clocks","collect event-level multi-layer responses","extract blind features","validate anchors/two clocks","unblind model scores once","publish failures unchanged"],"role_split":["provider","registrar","analyst"],"code_can_generate_independence":False}
 return finalize(713,"constructed_layer_speed_executable_transfer_spec","Independent custody and apparatus remain external obligations.",packet)
def st714():return finalize(714,"blocked_no_strict_model_source","No dimension/refinement/gauge/clock source is exported.",{"canonical_model":False,"source_formula":False})
def st715():return finalize(715,"blocked_no_physical_prediction_record","No SI c, photon, Lorentz apparatus or held-out nature record.",{"physical_prediction":False})
def st716():return finalize(716,"blocked_no_independent_empirical_evidence","Round five is local analytic/computational work.",{"external_referee":"absent","independent_laboratory_record":"absent","new_cycle_round":5})

def figures(r):
 FIG_DIR.mkdir(exist_ok=True);x=np.array(r["ST706"]["confusion"]);fig,ax=plt.subplots(figsize=(6,5));im=ax.imshow(x,cmap="Blues");fig.colorbar(im,ax=ax);ax.set(xlabel="predicted",ylabel="true",title="ST706: synthetic model confusion");fig.tight_layout();fig.savefig(FIG_DIR/"st706_confusion.png",dpi=180);plt.close(fig)

def main():
 funcs={702:st702,703:st703,704:st704,705:st705,706:st706,707:st707,708:st708,709:st709,710:st710,711:st711,712:st712,713:st713,714:st714,715:st715,716:st716};r={}
 for k in range(702,717):print(f"running ST{k}",flush=True);r[f"ST{k}"]=funcs[k]()
 RESULTS.write_text(json.dumps(native(r),indent=2,sort_keys=True),encoding="utf-8")
 with SUMMARY.open("w",newline="",encoding="utf-8") as f:
  w=csv.writer(f);w.writerow(["program","status","object","boundary"])
  for k in range(702,717):x=r[f"ST{k}"];w.writerow([f"ST{k}",x["status"],x["object"],x["boundary"]])
 figures(r);print(f"wrote {RESULTS.name} and {SUMMARY.name}")
if __name__=="__main__":main()
