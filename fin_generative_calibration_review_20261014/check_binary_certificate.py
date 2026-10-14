#!/usr/bin/env python3
"""Exact serialized binary-event no-go checker; no numerical dependencies."""
import argparse
import copy
import json
from fractions import Fraction as F
from pathlib import Path


def check(c):
    scale=c['integer_scale'];assert scale==10**15
    probes=c['probe_laws'];assert len(probes)==6
    assert [p['model'] for p in probes]==['FIN_A7']*3+['FLAT_P7_TRACE']*3
    assert [p['gain_point'] for p in probes]==['2.97','3','3.03']*2
    rows=[]
    for p in probes:
        rr=[]
        assert len(p['probability_boxes'])==12
        for lo,hi in p['probability_boxes']:
            assert 0<=F(lo)<=F(hi)<=1
            a,b=F(lo)*scale,F(hi)*scale
            assert a.denominator==b.denominator==1
            rr.append((int(a),int(b)))
        assert sum(a for a,b in rr)<=scale<=sum(b for a,b in rr)
        rows.append(rr)
    def event(row,mask):
        inside=[x for d,x in enumerate(row) if mask>>d&1]
        outside=[x for d,x in enumerate(row) if not mask>>d&1]
        return max(sum(a for a,b in inside),scale-sum(b for a,b in outside)),min(sum(b for a,b in inside),scale-sum(a for a,b in outside))
    choices=c['event_choice_rows'];assert len(choices)==4096
    assert {row[0] for row in choices}==set(range(4096))
    for mask,i,j,saved in choices:
        assert 0<=i<3<=j<6
        a,b=event(rows[i],mask),event(rows[j],mask)
        exact=max(abs(a[0]-b[1]),abs(a[1]-b[0]))
        assert exact<=saved and F(saved,scale)<2*F(c['class']['output_TV_ball_each'])
    assert c['proper_13_output_events_covered']==2*(4096-2)+2==8190
    assert F(c['max_selected_event_difference_upper'])==F(max(r[3] for r in choices),scale)
    return 8190


def main():
    ap=argparse.ArgumentParser();ap.add_argument('certificate',type=Path);ap.add_argument('--output',type=Path)
    args=ap.parse_args();c=json.loads(args.certificate.read_text());count=check(c)
    broken=copy.deepcopy(c);broken['event_choice_rows'].pop()
    try:check(broken)
    except AssertionError:rejected=True
    else:raise AssertionError('missing event coverage accepted')
    result={"events_checked":count,"exact_arithmetic_pass":True,"omitted_event_mutation_rejected":rejected,"limit":"input probability enclosures inherit audited E interval/backend assumptions; output law-TV class, not measured detector"}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
