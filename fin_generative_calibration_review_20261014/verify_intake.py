#!/usr/bin/env python3
"""Immutable F intake: independent process/statistics checks and input probes.

Uses versioned safe_seed and the admitted E review, not the F producer for
scientific reconstruction. F helpers are imported only for interface probes.
No physical data are generated. Interval statements assume mpmath.iv.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import html
import importlib.util
import io
import json
import math
import re
import sys
import zipfile
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import expm
from scipy.special import softmax

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"fin_generative_calibration_campaign_20261014"
EREVIEW=ROOT/"fin_generative_operational_review_20261014"
SCALE=10**15


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def load(p,name):
    spec=importlib.util.spec_from_file_location(name,p)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def sample_certificate(e,rho,n):
    I=mp.iv;I.dps=65
    eps=I.mpf(rho)/3
    b=lambda k:8190*I.exp(-2*k*eps**2)
    assert e.endpoints(b(n))[1]<F(1,40)<e.endpoints(b(n-1))[0]
    radius=I.sqrt(I.log(I.mpf(8190)/I.mpf('.025'))/(2*n))
    return {"rho":rho,"n_min_for_this_union_bound":n,
            "tail_bound_at_n":e.box(b(n)),"tail_bound_at_n_minus_one":e.box(b(n-1)),
            "empirical_radius_interval":e.box(radius),"eps_choice":"rho/3",
            "assumptions":"13 outcomes, iid, actual law in independently frozen hypothesis class",
            "arithmetic":"mpmath.iv; exact decimal parameters; outward rational endpoint serialization"}


def binary_certificate(gain,static):
    probes=[]
    for name in ("FIN_A7","FLAT_P7_TRACE"):
        cells=next(x['boxes'] for x in gain['probability_boxes'] if x['model']==name)
        for g,row,ref in (("2.97",cells[0],"gain cover first cell"),
                          ("3",static['individual_law_enclosures'][name],"static point certificate"),
                          ("3.03",cells[-1],"gain cover last cell")):
            probes.append({"model":name,"gain_point":g,"probability_boxes":row,"origin":ref})
    scaled=[]
    for p in probes:
        row=[]
        for lo,hi in p['probability_boxes']:
            a,b=F(lo)*SCALE,F(hi)*SCALE
            assert a.denominator==b.denominator==1
            row.append((int(a),int(b)))
        scaled.append(row)
    def event(row,mask):
        inside=[x for d,x in enumerate(row) if mask>>d&1]
        outside=[x for d,x in enumerate(row) if not mask>>d&1]
        return max(sum(a for a,b in inside),SCALE-sum(b for a,b in outside)),min(sum(b for a,b in inside),SCALE-sum(a for a,b in outside))
    choices=[]
    for mask in range(4096):
        bounds=[event(row,mask) for row in scaled]
        possible=[]
        for i in range(3):
            for j in range(3,6):
                a,b=bounds[i],bounds[j]
                upper=max(abs(a[0]-b[1]),abs(a[1]-b[0]))
                possible.append((upper,i,j))
        upper,i,j=min(possible)
        assert F(upper,SCALE)<F(1,100)
        choices.append([mask,i,j,upper])
    return {"claim":"No single fixed binary event uniformly distinguishes the two hypothetical static law classes",
            "class":{"gains":"independent [2.97,3.03]","output_TV_ball_each":"0.005","nominal_invalid":0},
            "proof":"For each event choose admissible gain points with event-probability difference < .01. The projected TV-ball Bernoulli intervals overlap; each common probability is attainable by mass transfer between event and complement.",
            "probe_laws":probes,"integer_scale":SCALE,"event_choice_rows":choices,
            "valid_label_masks_covered":4096,"proper_13_output_events_covered":8190,
            "max_selected_event_difference_upper":str(F(max(x[3] for x in choices),SCALE)),
            "limits":["per-law output-TV balls, not an imposed shared calibrated channel","one fixed event only; not multivalued statistics, several events or tighter calibration","input enclosures inherited from audited E mpmath.iv certificates"]}


def interface_probes(scratch):
    parser=load(scratch/'validate_records.py','f_parser')
    def fixture(overrides=None,extra=None,empty=False):
        row={k:'synthetic' for k in parser.FIELDS}
        row.update(record_id='r1',dataset_role='calibration',trajectory_id='t',replicate_id='rep',reset_id='reset',elapsed_s='0',latency_s='0',timestamp_utc='2026-10-14T00:00:00Z',outcome='0',validity_code='VALID',phase_a_prepared='0',phase_b_prepared='0',absolute_a_read='0',absolute_b_read='0')
        row.update(overrides or {})
        out=io.StringIO();writer=csv.DictWriter(out,fieldnames=parser.FIELDS);writer.writeheader()
        if not empty:writer.writerow(row)
        if extra:
            other=row.copy();other.update(extra);writer.writerow(other)
        return out.getvalue()
    cases={"valid_synthetic":fixture(),"empty_header":fixture(empty=True),
           "nan_elapsed":fixture({'elapsed_s':'nan'}),"infinite_latency":fixture({'latency_s':'inf'}),
           "blank_record_id":fixture({'record_id':''}),"malformed_timestamp":fixture({'timestamp_utc':'not-a-date'}),
           "same_reset_across_roles":fixture(extra={'record_id':'r2','dataset_role':'validation'}),
           "negative_elapsed":fixture({'elapsed_s':'-1'}),
           "invalid_without_reason":fixture({'outcome':'INVALID','validity_code':''}),
           "noisy_difference":fixture({'absolute_b_read':'1','outcome':'2'})}
    output={}
    for name,content in cases.items():
        parser.open=lambda *a,**k:io.StringIO(content)
        try:
            result=parser.parse('synthetic_fixture.csv')
            output[name]={"accepted":True,"rows":result['rows'],"warnings":result['warnings'],"iid_certified":result['independent_iid_shots_certified']}
        except ValueError as ex:
            output[name]={"accepted":False,"error":str(ex)}
    assert output['nan_elapsed']['accepted'] and not output['negative_elapsed']['accepted']
    assert not output['invalid_without_reason']['accepted'] and output['noisy_difference']['warnings']
    # Import prototype only to check API gates; not an independent proof route.
    sys.path.insert(0,str(scratch))
    decision=load(scratch/'decision_rule.py','f_decision')
    emp=np.r_[softmax(1.5*decision.all_models()['FIN_A7'][0]),0]
    nprobes={}
    for name,n in (("nan",float('nan')),("infinity",float('inf')),("fractional",10.5)):
        try:nprobes[name]={"accepted":True,"outcome":decision.classify(emp,n)['outcome']}
        except (ValueError,AssertionError,OverflowError) as ex:nprobes[name]={"accepted":False,"error":str(ex)}
    return {"parser":output,"sample_size_API":nprobes,"scope":"synthetic in-memory interface probes, zero physical records"}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--replay',type=Path,required=True)
    args=ap.parse_args()
    entries=[]
    for line in (SOURCE/'MANIFEST.sha256').read_text().splitlines():
        h,p=line.split(None,1);assert sha(SOURCE/p)==h,p;entries.append(p)
    uncovered=[str(p.relative_to(SOURCE)) for p in SOURCE.rglob('*') if p.is_file() and p.name!='MANIFEST.sha256' and str(p.relative_to(SOURCE)) not in entries]
    assert not uncovered,uncovered
    stage=0
    for path in SOURCE.glob('RG-F0[1-5]/INPUTS.sha256'):
        for line in path.read_text().splitlines():
            h,p=line.split(None,1);assert sha(SOURCE/p)==h;stage+=1
    with zipfile.ZipFile(SOURCE/'source_context/FIN_RG_E01_E05_FULL_HANDOFF_20261014.zip') as z:
        members=[n for n in z.namelist() if not n.endswith('/')]
        assert all((ROOT/n).is_file() and hashlib.sha256(z.read(n)).hexdigest()==sha(ROOT/n) for n in members)
    assert sha(SOURCE/'source_model.py')==sha(ROOT/'fin_generative_operational_campaign_20261014/source_model.py')
    assert sha(SOURCE/'source_context/STATIC_INTERVAL_CERTIFICATE.json')==sha(ROOT/'fin_generative_operational_campaign_20261014/RG-E03/STATIC_INTERVAL_CERTIFICATE.json')
    e=load(EREVIEW/'verify_intake.py','admitted_e')
    gain_checker=load(EREVIEW/'check_gain_certificate.py','gain_checker')
    gain=read(EREVIEW/'GAIN_BOX_CERTIFICATE.json');assert gain_checker.check(gain)==3600
    core=ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py'
    module=load(core,'versioned_core');A=module.build_rank7()[-1];models=module.build_countermodels(A)
    process=read(SOURCE/'RG-F01/PROCESS_CONTRACT.json')['positive_quotients']
    mu=np.arange(1,13,dtype=float);mu/=mu.sum()
    K=np.zeros((13,12));K[np.arange(12),np.arange(12)]=.95;K[(np.arange(12)+1)%12,np.arange(12)]=.02;K[12]=.03
    numerical=[]
    for name in ('FIN_A7','FULL_POTTS_TRACE','FLAT_P7_TRACE'):
        for kind in ('heat_bath','metropolis','barker'):
            B=e.quotient(models[name][0],kind)
            saved=process[name][kind]
            assert np.max(np.abs(B-saved['generator_12x12']))<1e-13
            p1=mu@expm(.21*B);P=expm(.57*B)
            joint=np.array([(p1*K[o])@P@K.T for o in range(13)])
            delta=float(np.max(np.abs(joint-saved['synthetic_noisy_joint_13x13'])))
            assert delta<1e-13 and abs(joint.sum()-1)<1e-13
            p=mu@expm(.12*B);P1=expm(.16*B);P2=expm(.32*B)
            tri=np.array([[((p*K[i])@P1*K[j])@P2@K.T for j in range(13)] for i in range(13)])
            assert tri.min()>=0 and abs(tri.sum()-1)<1e-12
            shorter=np.array([(p*K[i])@expm(.48*B)@K.T for i in range(13)])
            assert np.max(np.abs(tri.sum(axis=1)-shorter))<1e-13
            numerical.append({"model":name,"kinetic":kind,"noisy_joint_max_discrepancy":delta,"three_time_mass":float(tri.sum()),"three_time_min":float(tri.min())})
    # Independent field witness, without copying the producer's saved PASS.
    L=np.zeros((144,144));ds=np.array([(b-a)%12 for a,b in itertools_product()])
    for a,b in itertools_product():
        x=12*a+b
        r1=softmax(1.5*A[:,b]+.01*(np.arange(12)==0));r2=softmax(1.5*A[:,a])
        for j in range(12):
            if j!=a:L[x,12*j+b]+=r1[j]
            if j!=b:L[x,12*a+j]+=r2[j]
    np.fill_diagonal(L,-L.sum(axis=1));blocks=L@np.eye(12)[ds]
    defect=max(float(np.max(np.abs(blocks[ds==d]-blocks[np.flatnonzero(ds==d)[0]]))) for d in range(12))
    assert abs(defect-.002711545661289394)<1e-12
    samples=[sample_certificate(e,'0.00119114',40278676),sample_certificate(e,'0.00039114',373539331)]
    binary=binary_certificate(gain,read(SOURCE/'source_context/STATIC_INTERVAL_CERTIFICATE.json'))
    binary['source_hashes']={str((EREVIEW/'GAIN_BOX_CERTIFICATE.json').relative_to(ROOT)):sha(EREVIEW/'GAIN_BOX_CERTIFICATE.json'),str((SOURCE/'source_context/STATIC_INTERVAL_CERTIFICATE.json').relative_to(ROOT)):sha(SOURCE/'source_context/STATIC_INTERVAL_CERTIFICATE.json')}
    args.output.with_name('BINARY_EVENT_NO_GO_CERTIFICATE.json').write_text(json.dumps(binary,indent=2)+'\n')
    probes=interface_probes(args.replay)
    normalize=lambda s:re.sub(r'\s+','',html.unescape(s).replace('\\','').replace('*',''))
    prior_dirty='08b6d235 plus preceding uncommitted RG-E intake/AGENTS additions'
    out={"date":"2026-10-14","baseline":prior_dirty,
         "versions":{"python":sys.version.split()[0],"numpy":np.__version__,"scipy":__import__('scipy').__version__,"mpmath":mp.__version__},
         "manifest":{"files":len(entries),"mismatches":0,"stage_entries":stage,"uncovered":uncovered},
         "predecessor_archive_members_byte_matched":len(members),"inherited_gain_cover_bound_locally_bound":True,
         "plan_byte_identical":sha(ROOT/'FIN_GENERATIVE_CALIBRATION_NEXT_CHAIN_20261014.md')==sha(SOURCE/'source_context/FIN_N2_CALIBRATION_NEXT_CHAIN_USER.md'),
         "plan_format_normalized_identical":normalize((ROOT/'FIN_GENERATIVE_CALIBRATION_NEXT_CHAIN_20261014.md').read_text())==normalize((SOURCE/'source_context/FIN_N2_CALIBRATION_NEXT_CHAIN_USER.md').read_text()),
         "input_hashes":{"checker":sha(Path(__file__)),str(core.relative_to(ROOT)):sha(core),'AGENTS.md':sha(ROOT/'AGENTS.md'),'source_manifest':sha(SOURCE/'MANIFEST.sha256')},
         "independent_process_numeric":numerical,"phase_field_witness":defect,
         "statistical_integer_certificates":samples,"interface_probes":probes,
         "binary_event_no_go_summary":{k:v for k,v in binary.items() if k not in ('event_choice_rows','probe_laws')},
         "source_replay":{"directory":str(args.replay),"verifier":read(args.replay/'INDEPENDENT_VERIFICATION.json')},
         "limits":["zero physical records", "numerical likelihood helper not a globally outward distance oracle", "probability/threshold intervals assume mpmath.iv", "synthetic .05 emission error is not target .005 calibration", "overall calibration confidence coverage missing"]}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({"manifest":out['manifest'],"process_cases":len(numerical),"sample_counts":[s['n_min_for_this_union_bound'] for s in samples],"binary_event_bound":binary['max_selected_event_difference_upper'],"source_checks":[out['source_replay']['verifier']['checks_passed'],out['source_replay']['verifier']['checks_total']],"plan_equivalent":out['plan_format_normalized_identical']},indent=2))


def itertools_product():
    return ((a,b) for a in range(12) for b in range(12))


if __name__=='__main__':main()
