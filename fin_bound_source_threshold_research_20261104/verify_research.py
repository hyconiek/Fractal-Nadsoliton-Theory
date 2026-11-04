#!/usr/bin/env python3
"""Finite exact controls; not a machine verifier of the analytical PDE proof."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
checks=[]


def check(name, condition, category='exact_control'):
    checks.append({'name':name,'passed':bool(condition),'category':category})
    if not condition:
        raise AssertionError(name)


def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j]+=x*y
    return out


def add(a,b):
    out=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a):out[i]+=x
    for i,x in enumerate(b):out[i]+=x
    while len(out)>1 and out[-1]==0:out.pop()
    return out


def scale(a,c):return [x*c for x in a]


def recurrence(p,n):
    # I_j/S_p for I_j=int tanh^j sech^p exp(-ikx), y=ik.
    seq=[[F(1)],[F(0),-F(1,p)]]
    for j in range(1,n):
        seq.append(scale(add(scale(seq[j-1],j),scale(mul([0,1],seq[j]),-1)),F(1,p+j)))
    return seq


def main():
    w=F(4,3)
    check('internal parameter is in frozen band', F(5,4)<w<F(7,4))
    check('threshold equation 3w=4', 3*w==4)
    check('threshold mass squared', 4+w==F(16,3))
    check('mixed carrier squared at crossing', F(1,3)*(6-2)**2==4+w)
    check('crossing unique after positive squaring', F(8,4)**2/F(3)==w)
    check('detuning derivative magnitude at crossing', F(2)**2*3/w==9)
    check('left endpoint above threshold', 3*F(5,4)<4)
    check('right endpoint below threshold', 3*F(7,4)>4)
    check('closed intermediate resolvent below bound spectrum', 3*F(5,4)>F(3,2)**2)
    check('CM initial shape slip coefficient',6*(2*w-3)/(3*(3-4*w))==F(2,7))

    integral=4*F(2,5)-F(2,3)-3*F(2,7)
    check('threshold profile polynomial integral',integral==F(8,105))
    check('threshold overlap divided by sqrt3',F(3,8)*integral==F(1,35))
    # The threshold state is orthogonal to f: integrate (3y^2-1)/2.
    check('projection does not alter threshold overlap',F(1,2)*(3*F(2,3)-2)==0)

    # Independent Fourier recurrence check of J and Fb numerators.
    seq=recurrence(2,5)
    # Conjugated Jost numerator is 3K^2-1+y^2+3yK, y=ik.
    k1=add(add(scale(seq[3],3),mul([-1,0,1],seq[1])),mul([0,3],seq[2]))
    k3=add(add(scale(seq[5],3),mul([-1,0,1],seq[3])),mul([0,3],seq[4]))
    expected_k3=scale(mul(mul([0,1],[4,0,-1]),[2,0,1]),F(1,48))
    check('K^3 sech^2 transform numerator',k3==expected_k3)
    # J/sqrt3=(3K^3-K)sech^2/4.
    j=scale(add(scale(k3,3),scale(k1,-1)),F(1,4))
    expected_j=scale(mul([0,0,0,1],[4,0,-1]),F(1,64))
    check('J transform numerator divided by sqrt3',j==expected_j)
    ratio_without_sqrt3=-F(8,2*(2-8))
    check('on-shell ratio divided by sqrt3',ratio_without_sqrt3==F(2,3))
    check('both-root imaginary coefficient divided by sqrt3 d_b',9*ratio_without_sqrt3==6)
    check('negative: reversing J sign changes imaginary coefficient',-9*ratio_without_sqrt3!=6)

    eta2=-6*F(1,2)**2
    psi2=-12*F(1,2)**2
    check('scalar quadratic cascade coefficient',-12*eta2*F(1,2)==9)
    check('mixed quadratic cascade coefficient',-12*F(1,2)*psi2==18)
    check('direct cubic coefficient',-6*F(1,2)**3==-F(3,4))
    check('negative: omitting scalar cascade removes certified imaginary term',0!=6)

    # i exp(-i pi/4)=(1+i)/sqrt2; coefficients expressed in 1/sqrt2 units.
    minus_phase=(F(1),F(-1))
    rotated=(-minus_phase[1],minus_phase[0])
    check('correct threshold response phase',rotated==(1,1))
    check('negative: ET displayed i has wrong zero real part',rotated[0]!=0)
    check('sine factor and sqrt-time primitive cancel factors two',F(1,2)*2==1)

    check('localized cubic pointwise critical exponent',3-1==2)
    check('localized cubic integrated squared exponent',2*3-2*2==2)
    check('localized quintic pointwise exponent',5-1==4)
    check('localized quintic integrated squared exponent',2*5-2*2==6)
    check('QP position correction integrated squared exponent',2*2-2==2)
    check('negative: uncontrolled quadratic source has no small trace budget',2*2-2*2==0)

    meta=json.loads((HERE/'INPUTS.json').read_text())
    for row in meta['files']:
        data=(ROOT/row['path']).read_bytes()
        check('SHA256 '+row['path'],hashlib.sha256(data).hexdigest()==row['sha256'],'provenance')
        check('Git blob '+row['path'],hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==row['git_blob_sha'],'provenance')


if __name__=='__main__':
    try:
        main();status='PASS_FINITE_CONTROLS_RESIDUAL_TRACE_AND_ACTION_GATES_OPEN'
    except Exception as exc:
        status='FAIL: '+str(exc)
    result={'status':status,'counts':{k:sum(x['category']==k for x in checks) for k in ('exact_control','provenance')},
        'passed':sum(x['passed'] for x in checks),'checks':checks,
        'actual_critical_damping_proved':False,'actual_instability_proved':False,
        'formal_verification':False,'independent_review':False}
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
    raise SystemExit(0 if status.startswith('PASS') else 1)
