import numpy as np, math
from math import factorial
Q=12; D=7; j=np.arange(Q)
# basis/operator identical to C++
W=np.array([[0.0 if a==b else math.cos(0.18575*min(abs(a-b),Q-abs(a-b))+0.1625)/(1+min(abs(a-b),Q-abs(a-b))**1.8) for b in range(Q)] for a in range(Q)])
L=np.diag(W.sum(1))-W; lam=np.fft.fft(L[0]).real[:7]
E=[]; sectors=[]
for k in (3,4,5):
    E += [np.sqrt(2/Q)*np.cos(2*np.pi*k*j/Q),np.sqrt(2/Q)*np.sin(2*np.pi*k*j/Q)]; sectors += [k,k]
E += [((-1.)**j)/np.sqrt(Q)]; sectors += [6]
E=np.column_stack(E); Ad=np.array([lam[k] for k in sectors]);
PV=E@E.T; PH=np.eye(Q)-np.ones((Q,Q))/Q-PV
# degree4 alphas in same recursion order
def gen(rem,pos,cur,D):
    if pos==D-1:
        yield tuple(cur+[rem]); return
    for v in range(rem+1): yield from gen(rem-v,pos+1,cur+[v],D)
alphas=list(gen(4,0,[],D)); idx={a:i for i,a in enumerate(alphas)}
# read actual
rows=[]
with open('/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_06/quartic_tensor_loo_actual.txt') as f:
    head=f.readline()
    for line in f:
        left,*vals=line.split(); aa=tuple(map(int,left.split(','))); rows.append((aa,*map(float,vals)))
assert [r[0] for r in rows]==alphas
actual=np.array([[r[1],r[2],r[3]] for r in rows]).T
# quadratic polynomial coefficients in z for each node
# conventional coefficients: q(z)=sum_alpha q_alpha z^alpha, no multinomial factor
qalphas=list(gen(2,0,[],D)); qidx={a:i for i,a in enumerate(qalphas)}
Bphi=np.einsum('nr,ns->nrs',E,E)
BA=np.einsum('nr,ns,s->nrs',E,E,Ad); BA=(BA+BA.swapaxes(1,2))/2
H=np.einsum('nm,mrs->nrs',PH,Bphi)
K=np.einsum('nm,mrs->nrs',PH,BA)
def quad_coeff(Bnode):
    out=np.zeros(len(qalphas))
    for qi,a in enumerate(qalphas):
        inds=[]
        for r,p in enumerate(a): inds += [r]*p
        r,s=inds
        out[qi]=Bnode[r,s]*(2 if r!=s else 1)
    return out
def qmul(a,b):
    out=np.zeros(len(alphas))
    for ia,aa in enumerate(qalphas):
        ca=a[ia]
        if ca==0: continue
        for ib,bb in enumerate(qalphas):
            cb=b[ib]
            if cb==0: continue
            cc=tuple(aa[r]+bb[r] for r in range(D)); out[idx[cc]]+=ca*cb
    return out
pred=np.zeros((3,len(alphas)))
for n in range(Q):
    h=quad_coeff(H[n]); k=quad_coeff(K[n])
    pred[0]+=(-12/Q)*qmul(h,h)
    pred[1]+=(2/Q)*qmul(h,k)
res=actual-pred
print('lambda',lam[3:7])
for g in range(3):
    im=int(np.argmax(np.abs(res[g])))
    print('gdeg',g,'max_abs',format(np.max(np.abs(res[g])),'.17g'),'alpha',alphas[im],'actual',format(actual[g,im],'.17g'),'pred',format(pred[g,im],'.17g'))
print('g2_actual_max',format(np.max(np.abs(actual[2])),'.17g'))
# relative on nonzero
for g in (0,1):
    nz=np.abs(pred[g])>1e-10
    print('gdeg',g,'max_rel_nonzero',np.max(np.abs(res[g,nz])/np.abs(pred[g,nz])), 'nonzero',nz.sum())
np.savez('/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_06/quartic_tensor_loo_certificate.npz',actual=actual,pred=pred,res=res,alphas=np.array(alphas),lam=lam)
