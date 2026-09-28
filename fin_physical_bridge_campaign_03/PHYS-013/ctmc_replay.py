#!/usr/bin/env python3
import numpy as np
from pathlib import Path
import importlib.util

# Requires the campaign-I fin_core.py or the user's uploaded seed implementation.
FIN_CORE=Path(__file__).resolve().parents[2]/"fin_core.py"
if not FIN_CORE.exists():
    FIN_CORE=Path("/mnt/data/fin_core.py")
spec=importlib.util.spec_from_file_location("fc",FIN_CORE)
fc=importlib.util.module_from_spec(spec); spec.loader.exec_module(fc)

g=3.0
A,_=fc.A7_default()
E=-(g/2)*A
target=np.exp((g/2)*A[0]); target/=target.sum()

rows=[
(27.21616378414087,0.002124831534156879),(27.26326682273862,0.0061154645400471885),
(27.41337226127547,-0.008556839939843908),(26.906115479047102,0.0022846785213488374),
(27.557525159448726,0.006837583062976993),(26.805303960497984,-0.0017595959711329545),
(27.644025711899264,0.005340759246436733),(26.89875735630125,-0.00850552523438075),
(27.744031935187554,0.00492552647818556),(26.677421561963182,0.0029549696237285428),
(27.836526043158187,0.005667204477395771),(26.637860287428314,0.00861021470344836),
(28.056035735885136,0.0),(27.988627437412525,0.005096276013137002),
(26.33227853641043,-0.00036473768052680544),(26.447364855675033,-0.006952207273967614),
(26.451423068324804,0.009160266581268894),(26.20491248228989,-0.0028026893922706853),
(28.190333086266843,-0.008712116812113369),(28.26414117021333,-0.008080994075530867),
(26.24539960819259,0.008957948219677103),(26.15516857745329,0.004953831863102609),
(28.321715396935492,-0.007595471796703723),(28.587646568523265,0.00030121384198267265)]
T=np.array([x[0] for x in rows]); B=1+np.array([x[1] for x in rows])
T1,T2=T[:12],T[12:]; B1,B2=B[:12],B[12:]

Q=np.zeros((144,144))
idx=lambda i,j:12*i+j
r1=np.zeros((12,12)); r2=np.zeros((12,12))
for a in range(12):
    w=np.rint(T1[a]*E[a,:])
    r1[a,:]=B1[a]*np.exp(-w/T1[a])
for b in range(12):
    w=np.rint(T2[b]*E[:,b])
    r2[b,:]=B2[b]*np.exp(-w/T2[b])
for i in range(12):
    for j in range(12):
        s=idx(i,j); total=0.
        for a in range(12):
            if a!=i:
                rr=r1[a,j]; Q[s,idx(a,j)]+=rr; total+=rr
        for b in range(12):
            if b!=j:
                rr=r2[b,i]; Q[s,idx(i,b)]+=rr; total+=rr
        Q[s,s]-=total

M=Q.T.copy(); M[-1,:]=1
rhs=np.zeros(144); rhs[-1]=1
pi=np.linalg.solve(M,rhs)
diff=np.zeros(12)
for i in range(12):
    for j in range(12):
        diff[(j-i)%12]+=pi[idx(i,j)]
tv=.5*np.abs(diff-target).sum()
print("TV",tv)
print("difference histogram",diff.tolist())
