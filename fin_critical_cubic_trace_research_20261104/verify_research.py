#!/usr/bin/env python3
"""Finite formula controls and previews; not a verifier of PDE estimates."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = []
previews = []


def check(name, value, category='exact_formula'):
    checks.append({'name':name, 'passed':bool(value), 'category':category})
    if not value:
        raise AssertionError(name)


def clean(d):
    return {k:v for k,v in d.items() if v}


def add(a, b):
    d = dict(a)
    for k,v in b.items():
        d[k] = d.get(k, F(0))+v
    return clean(d)


def dt(jet):
    # Key: (envelope derivative order, cos/sin, power of 1/T, power of m).
    out = {}
    for (n,trig,ti,mi),v in jet.items():
        key = (n+1,trig,ti+1,mi)
        out[key] = out.get(key,0)+v
        key = (n,'sin' if trig=='cos' else 'cos',ti,mi+1)
        out[key] = out.get(key,0)+(-v if trig=='cos' else v)
    return clean(out)


def simpson(fun, a, b, n=2048):
    step = (b-a)/n
    return step/3*(fun(a)+fun(b)+sum((4 if i%2 else 2)*fun(a+i*step) for i in range(1,n)))


def time_integral(M):
    # Exact elementary expression, evaluated in floating arithmetic.
    z = 2*M
    osc = math.sin(z)/z+4*math.cos(z)/z**2-12*math.sin(z)/z**3-24*math.cos(z)/z**4+24*math.sin(z)/z**5
    return 0.1+osc/2


def main():
    initial = {(0,'cos',0,0):F(1)}
    mass = {(0,'cos',0,2):F(1)}
    residual = add(dt(dt(initial)),mass)
    expected = {(2,'cos',2,0):F(1),(1,'sin',1,1):F(-2)}
    check('exact time residual and mass cancellation', residual==expected)
    wrong = {(2,'cos',2,0):F(1),(1,'sin',1,1):F(2)}
    check('negative: wrong oscillatory force sign rejected', residual!=wrong)
    check('negative: omitted envelope derivative rejected', residual!={(2,'cos',2,0):F(1)})
    # For h=(1+y^2)^-1, h'' numerator after denominator (1+y^2)^3.
    numerator = [-2,0,-2]
    numerator[2] += 8
    check('exact h second derivative numerator', numerator==[-2,0,6])
    check('zero initial position for a=s^2', F(0)**2==0)
    check('zero initial velocity for a=s^2', 2*F(0)==0)
    check('field L2 exponent', F(3,2)-F(1,2)==1)
    check('spatial derivative L2 exponent', F(3,2)+1-F(1,2)==2)
    check('forcing leading L2 exponent', F(3,2)+2-F(1,2)==3)
    check('forcing second-envelope term exponent', F(3,2)+4-F(1,2)==5)
    check('compact trace integrated exponent', 2*F(3,2)-2==1)
    check('CT5 critical exponent', 6-2*F(5,2)==1)
    check('negative: order epsilon is not order epsilon squared', 1-2<0)
    check('old pointwise integrated exponent', 6-3*2==0)
    check('low-frequency observation squared exponent', 1-F(1,2)==F(1,2))
    check('oscillatory integral constant', F(1,2)*F(1,5)==F(1,10))
    check('IBP boundary plus derivative error', F(1,4)+F(4,4)*F(1,4)==F(1,2))
    check('oscillatory lower bound at mT>=10', F(1,10)-F(1,20)==F(1,20))
    check('compact spatial coefficient lower bound', 2*F(1,2)*F(1,2)==F(1,2))
    check('trace lower constant', F(1,2)**2*F(1,20)==F(1,80))
    check('trace upper constant', 2**2*F(1,5)==F(4,5))
    check('epsilon<=1/4 implies mT>=32', F(2)/F(1,4)**2==32)
    check('leading action-error exponent stays insufficient', 2-1+F(1,2)==F(3,2))

    # Bound previews use the closed formula, not a simulated nonlinear PDE.
    for ep in (0.25,0.125,0.0625,0.03125):
        A = simpson(lambda x:1/(math.sqrt(1+x*x)*(1+(ep*x)**2)),-1,1)
        for m in (math.sqrt(5.25),math.sqrt(5.5),math.sqrt(5.75)):
            M=m/ep**2
            jt=time_integral(M)
            integral=ep*A*A*jt
            previews.append({'epsilon':ep,'mass':m,'spatial_coefficient':A,
                'time_integral':jt,'trace_integral':integral,
                'trace_over_epsilon':integral/ep,
                'trace_over_epsilon_squared':integral/ep**2})
            check('explicit trace preview '+str(ep)+' '+str(m),
                  ep/80<=integral<=4*ep/5,'floating_preview')
    M=math.sqrt(5.5)/(0.25**2)
    check('closed time integral versus separate quadrature',
          abs(time_integral(M)-simpson(lambda s:s**4*math.cos(M*s)**2,0,1,8192))<1e-9,
          'floating_preview')

    sources=json.loads((HERE/'INPUTS.json').read_text())
    for item in sources['files']:
        data=(ROOT/item['path']).read_bytes()
        check('SHA256 '+item['path'],hashlib.sha256(data).hexdigest()==item['sha256'],'provenance')
        blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        check('Git blob '+item['path'],blob==item['git_blob_sha'],'provenance')


if __name__=='__main__':
    try:
        main()
        status='PASS_FINITE_CONTROLS_ACTUAL_JOINT_GATE_OPEN'
    except Exception as exc:
        status='FAIL: '+str(exc)
    counts={k:sum(c['category']==k for c in checks) for k in ('exact_formula','floating_preview','provenance')}
    result={'status':status,'counts':counts,'passed':sum(c['passed'] for c in checks),
            'checks':checks,'previews':previews,'PDE_theorem_machine_checked':False,
            'actual_FIN_counterexample':False,'critical_damping_proved':False}
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('checks','previews')},indent=2))
    raise SystemExit(0 if status.startswith('PASS') else 1)
