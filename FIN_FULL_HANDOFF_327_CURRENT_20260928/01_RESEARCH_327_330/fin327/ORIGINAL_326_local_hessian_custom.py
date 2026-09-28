import json, math
from pathlib import Path
import numpy as np, mpmath as mp
OUT=Path('/mnt/data/fin326'); mp.mp.dps=80
# read root/B by importing constants from existing verifier without execution side effects too hard; reconstruct numpy same as model
Q=12;G=mp.mpf('5.145228719489142')
W=np.array([[0.0 if i==j else math.cos(0.18575*min(abs(i-j),Q-abs(i-j))+0.1625)/(1+min(abs(i-j),Q-abs(i-j))**1.8) for j in range(Q)] for i in range(Q)],float)
L=np.diag(W.sum(1))-W;lam=np.fft.fft(L[0]).real[:7];jj=np.arange(Q);cols=[]
for k in (3,4,5):cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*jj/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*jj/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**jj];Xf=np.column_stack(cols)
Tc=np.zeros((7,6));Tc[0,0]=1;Tc[1,1]=1;Tc[2,2]=.5;Tc[3,2]=math.sqrt(3)/2;Tc[4,3]=1;Tc[5,4]=1;Tc[6,5]=1
Bf=Xf@Tc; B=[[mp.mpf(str(x)) for x in row] for row in Bf]
root=json.load(open(OUT/'LOCAL_PATH_INTERVAL_CERT_326.json'))['root_y']; y=[mp.mpf(x) for x in root];rad=mp.mpf('.02')
lo=[x-rad for x in y];hi=[x+rad for x in y]
classes=[[j for j in range(Q) if j%3==a] for a in range(3)];c=[mp.mpf(1 if j%3==0 else -1 if j%3==1 else 0) for j in range(Q)]
hL=[];hU=[]
for j in range(Q):
 a=b=mp.mpf('0')
 for k in range(6):
  if B[j][k]>=0:a+=B[j][k]*lo[k];b+=B[j][k]*hi[k]
  else:a+=B[j][k]*hi[k];b+=B[j][k]*lo[k]
 hL.append(a);hU.append(b)
eL=[mp.exp(x) for x in hL];eU=[mp.exp(x) for x in hU]
SL=[sum(eL[j] for j in I) for I in classes];SU=[sum(eU[j] for j in I) for I in classes]
qL=mp.sqrt(SL[0]*SL[1]);qU=mp.sqrt(SU[0]*SU[1]);mL=qL/(2*qL+SU[2]);mU=qU/(2*qU+SL[2]);m2L=SL[2]/(2*qU+SL[2]);m2U=SU[2]/(2*qL+SU[2])
pL=[mp.mpf(0)]*Q;pU=[mp.mpf(0)]*Q
for a,I in enumerate(classes):
 cmL=mL if a<2 else m2L;cmU=mU if a<2 else m2U
 for i in I:
  pL[i]=cmL*eL[i]/(eL[i]+sum(eU[j] for j in I if j!=i))
  pU[i]=cmU*eU[i]/(eU[i]+sum(eL[j] for j in I if j!=i))
def ext(coeff,maxi):
 x=pL.copy();rem=mp.mpf(1)-sum(x);order=sorted(range(Q),key=lambda i:coeff[i],reverse=maxi)
 for i in order:
  if rem<=0:break
  add=min(rem,pU[i]-x[i]);x[i]+=add;rem-=add
 return sum(coeff[i]*x[i] for i in range(Q))
def mul(I,J):
 vals=[I[0]*J[0],I[0]*J[1],I[1]*J[0],I[1]*J[1]];return min(vals),max(vals)
mu=[];uu=[]
for a in range(6):
 co=[B[j][a] for j in range(Q)];mu.append((ext(co,False),ext(co,True)))
 co=[B[j][a]*c[j] for j in range(Q)];uu.append((ext(co,False),ext(co,True)))
cc=(2*mL,2*mU);HL=np.zeros((6,6));HU=np.zeros((6,6))
for a in range(6):
 for b in range(6):
  co=[B[j][a]*B[j][b] for j in range(Q)];ebb=(ext(co,False),ext(co,True));mm=mul(mu[a],mu[b]);qq=mul(uu[a],uu[b]); vals=[qq[0]/cc[0],qq[0]/cc[1],qq[1]/cc[0],qq[1]/cc[1]];div=(min(vals),max(vals));klo=ebb[0]-mm[1]-div[1];khi=ebb[1]-mm[0]-div[0]
  HL[a,b]=float((1/G if a==b else 0)-khi);HU[a,b]=float((1/G if a==b else 0)-klo)
mid=(HL+HU)/2;rr=(HU-HL)/2;lb=float(np.linalg.eigvalsh((mid+mid.T)/2)[0]-np.linalg.norm(rr,'fro')-1e-10)
out={'radius':0.02,'class_mass_bounds':[str(mL),str(mU)],'hessian_eigen_lower_bound':lb,'method':'high-precision p-box + exact LP extrema + midpoint/radius spectral enclosure','positive':lb>0}
(OUT/'LOCAL_HESSIAN_INTERVAL_326.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
