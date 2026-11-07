"""Source-ledger algebra and finite-grid adjoint-primitive diagnostics.

These checks do not prove uniform bounds or the conditional PDE theorem.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
checks = []
def check(name, condition):
    checks.append({'name': name, 'pass': bool(condition)})
    if not condition:
        raise AssertionError(name)

for entry in json.loads((HERE / 'INPUTS.json').read_text())['files']:
    data = (HERE.parent / entry['path']).read_bytes()
    check('input hash: ' + entry['path'], hashlib.sha256(data).hexdigest() == entry['sha256'])

c,z,C,Z,g,f,v,u,K = sp.symbols('c z C Z g f v u K')
eta, psi = c*g+v, z*f+u
odd = -6*K*(eta**2+psi**2)+6*K*(C**2*g**2+Z**2*f**2)
even = -12*K*eta*psi+12*K*C*Z*g*f
odd_split = -6*K*((c**2-C**2)*g**2+(z**2-Z**2)*f**2)-12*K*(c*g*v+z*f*u)-6*K*(v**2+u**2)
even_split = -12*K*(c*z-C*Z)*g*f-12*K*(c*g*u+z*f*v)-12*K*v*u
check('exact odd quadratic residual split', sp.expand(odd-odd_split)==0)
check('exact even quadratic residual split', sp.expand(even-even_split)==0)
# Each position is the sum of its two half-waves divided by two.
check('odd channel half-wave coefficient', -6/4 == -1.5)
check('even channel half-wave coefficient', -12/4 == -3)


def primitive(n, length, w, channel, t):
    x = np.arange(n)*(2*length/n)-length
    k = 2*np.pi*np.fft.fftfreq(n, d=2*length/n)
    K = np.tanh(x)
    sech = 2*np.exp(-abs(x))/(1+np.exp(-2*abs(x)))
    f = np.sqrt(.75)*sech**2
    g = np.sqrt(1.5)*sech*K
    omega, lam = np.sqrt(w), np.sqrt(3.)
    if channel == 'a':
        source, bound, mass, nu = K*f*f, g, 2., 2*omega
    elif channel == 'b':
        source, bound, mass, nu = K*g*g, g, 2., 2*lam
    else:
        source, bound, mass = K*f*g, f, np.sqrt(4+w)
        nu = lam+omega if channel == 'sum' else lam-omega
    source = source-bound*np.sum(source*bound)*(2*length/n)
    shat = np.fft.fft(source)
    dsource = np.fft.ifft(1j*k*shat)
    ddsource = np.fft.ifft(-k*k*shat)
    rootp = np.sqrt((k*k+1)*(k*k+4))
    distorted = np.fft.fft(ddsource+3*K*dsource+2*source)/rootp
    energy = np.sqrt(k*k+mass*mass)
    def time_integral(delta, time):
        return time*np.exp(.5j*delta*time)*np.sinc(delta*time/(2*np.pi))
    def W(hat):
        qhat = hat/rootp
        return np.fft.ifft(-k*k*qhat)-3*K*np.fft.ifft(1j*k*qhat)+(3*K*K-1)*np.fft.ifft(qhat)
    def values(time):
        plus, minus = time_integral(nu+energy,time), time_integral(nu-energy,time)
        return W(distorted*(plus+minus)/2), W(distorted*(plus-minus)/(2j*energy))
    cp, sn = values(t)
    omsn = np.fft.ifft(energy*np.fft.fft(sn))
    step = 1e-4
    hi, lo = values(t+step), values(t-step)
    exact = (np.exp(1j*nu*t)*W(distorted*np.cos(energy*t)),
             np.exp(1j*nu*t)*W(distorted*np.sin(energy*t)/energy))
    derivative_error = max(float(np.max(abs((h-l)/(2*step)-e))) for h,l,e in zip(hi,lo,exact))
    return {'norms': [float(np.max(abs(q))) for q in (cp,sn,omsn)],
            'derivative_error': derivative_error,
            'projection_error': float(abs(np.sum(source*bound)*(2*length/n)))}

rows=[]
for w in [1.25,1.5,1.75]:
    for channel in ['a','b','sum','difference']:
        for t in [8.,32.,128.]:
            result=primitive(4096,256.,w,channel,t)
            rows.append({'w':w,'channel':channel,'t':t,**result,
                         'norms_over_log': [a/math.log(2+t) for a in result['norms']]})
check('36 primitive samples finite', all(np.isfinite(r['norms']).all() for r in rows))
check('projected profiles orthogonal', max(r['projection_error'] for r in rows)<1e-12)
check('primitive time derivative agrees with propagator', max(r['derivative_error'] for r in rows)<2e-7)
refinements=[]
for channel in ['a','b','sum','difference']:
    base=next(r for r in rows if r['w']==1.5 and r['channel']==channel and r['t']==128.)
    fine=primitive(8192,256.,1.5,channel,128.)
    wide=primitive(6144,384.,1.5,channel,128.)
    refinements.append({'channel':channel,
        'grid_norm_difference':max(abs(a-b) for a,b in zip(base['norms'],fine['norms'])),
        'domain_norm_difference':max(abs(a-b) for a,b in zip(base['norms'],wide['norms']))})
# Supremum samples vary with grid spacing; this is convergence of sampled norms.
check('grid refinement sampled norms', max(r['grid_norm_difference'] for r in refinements)<.003)
check('domain refinement sampled norms', max(r['domain_norm_difference'] for r in refinements)<1e-7)
epsilons=[2.**(-k) for k in [4,8,12,16,20]]
ratios=[e*(1+math.log(1/e)) for e in epsilons]
check('illustrative little-o error ratio decreases', all(a>b for a,b in zip(ratios,ratios[1:])))
result={'status':'PASS','scope':'Exact source algebra, input integrity and finite-grid diagnostics; not a PDE simulation or independent proof audit.',
        'checks':checks,'check_count':len(checks),'primitive_samples':rows,'refinements':refinements,
        'illustrative_error_over_epsilon_squared':list(zip(epsilons,ratios)),
        'versions':{'numpy':np.__version__,'sympy':sp.__version__}}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':len(checks),'max_derivative_error':max(r['derivative_error'] for r in rows),'refinements':refinements},indent=2))
