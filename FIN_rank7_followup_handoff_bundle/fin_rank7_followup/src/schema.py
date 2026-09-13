from __future__ import annotations
from fractions import Fraction
import json, re

EXEC_STATUSES={'NOT_STARTED','RUNNING','DONE','RESOURCE_STOP','BLOCKED_DEPENDENCY','SUPERSEDED_BY_COUNTEREXAMPLE','OUT_OF_SCOPE'}
CLAIM_STATUSES={'EXACT_PROVED','INTERVAL_CERTIFIED','NUMERICAL_REPRODUCED','NUMERICAL_NEW','COUNTEREXAMPLE_CERTIFIED','COUNTEREXAMPLE_NUMERICAL','CONDITIONAL_LEMMA','UNRESOLVED','NOT_A_SCIENTIFIC_CLAIM',
'EXISTING_ACCEPTED_INPUT','ANALYTIC_IDENTITY_IMPORTED','EXISTING_ACCEPTED_SCOPED','EXACT_INTERVAL_CERTIFIED_UPSTREAM','EXISTING_ACCEPTED','NUMERICAL_REPRODUCED_STARTING','EXACT_STRICT_INTERVAL_ACCEPTED','NUMERICAL_IMPORTED','COUNTEREXAMPLE_CERTIFIED_STARTING','UNRESOLVED','EXACT_REPO_THEOREM','COUNTEREXAMPLE_CERTIFIED_REPO','EXACT_REPO_NOGO_SCOPED','EXISTING_SEPARATE_SCOPED'}

def _need(d,k,where):
    if k not in d or d[k] in (None,'',[]): raise ValueError(f'{where}: missing {k}')

def validate_task(t):
    for k in ['id','title','execution_status','claim_status','task_text','accept_output']:_need(t,k,'task')
    if not re.fullmatch(r'R7P-\d{3}',t['id']): raise ValueError('task: bad id')
    if t['execution_status'] not in EXEC_STATUSES: raise ValueError('task: bad execution_status')
    if t['claim_status'] not in CLAIM_STATUSES: raise ValueError('task: bad claim_status')
    return True

def validate_claim(c):
    for k in ['id','statement','status','dimension','domain','finite_or_limit','source','evidence_level']:_need(c,k,'claim')
    return True

def rational_endpoint(x):
    # Exact/certified endpoints must be integer ratio strings or {num,den}; decimal floats are forbidden.
    if isinstance(x,dict) and set(x)>= {'num','den'}:
        if not isinstance(x['num'],int) or not isinstance(x['den'],int) or x['den']==0: raise ValueError('bad rational endpoint')
        return Fraction(x['num'],x['den'])
    if isinstance(x,str) and re.fullmatch(r'-?\d+/[1-9]\d*',x): return Fraction(x)
    raise ValueError('certified endpoint must be rational, not float/decimal')

def validate_interval(iv, certified=True):
    _need(iv,'lo','interval'); _need(iv,'hi','interval')
    if certified:
        lo=rational_endpoint(iv['lo']); hi=rational_endpoint(iv['hi'])
    else:
        lo=float(iv['lo']); hi=float(iv['hi'])
    if lo>hi: raise ValueError('interval lo>hi')
    return True

def validate_certificate(c):
    for k in ['id','claim_id','domain','quantifiers','assumptions','proof_type','inputs','conclusion']:_need(c,k,'certificate')
    if c.get('global_pass') and c.get('unresolved_leaves'):
        raise ValueError('global PASS cannot have unresolved leaves')
    for iv in c.get('certified_intervals',[]): validate_interval(iv, certified=True)
    return True

def load_json(path):
    with open(path) as f:return json.load(f)
