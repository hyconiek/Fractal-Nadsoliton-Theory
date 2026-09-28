import mpmath as mp, numpy as np, json
mp.mp.dps=90; mp.iv.dps=70; iv=mp.iv
Q=12; g=mp.mpf('5.145228719489142')
# exact high-precision kernel
W=[[mp.mpf(0) if i==j else mp.cos(mp.mpf('0.18575')*min(abs(i-j),Q-abs(i-j))+mp.mpf('0.1625'))/(1+mp.mpf(min(abs(i-j),Q-abs(i-j)))**mp.mpf('1.8')) for j in range(Q)] for i in range(Q)]
r=sum(W[0]);l0=[(r if j==0 else mp.mpf(0))-W[0][j] for j in range(Q)]
lam=[sum(l0[j]*mp.cos(2*mp.pi*k*j/Q) for j in range(Q)) for k in range(7)]
X=[]
for j in range(Q):
    row=[]
    for k in (3,4,5): row += [mp.sqrt(lam[k]/6)*mp.cos(2*mp.pi*k*j/Q),mp.sqrt(lam[k]/6)*mp.sin(2*mp.pi*k*j/Q)]
    row += [mp.sqrt(lam[6]/12)*((-1)**j)];X.append(row)
A=[[sum(X[i][a]*X[j][a] for a in range(7)) for j in range(Q)] for i in range(Q)]
Tc=[[mp.mpf(0) for _ in range(6)] for __ in range(7)];Tc[0][0]=1;Tc[1][1]=1;Tc[2][2]=mp.mpf('.5');Tc[3][2]=mp.sqrt(3)/2;Tc[4][3]=1;Tc[5][4]=1;Tc[6][5]=1
B=[[sum(X[j][a]*Tc[a][k] for a in range(7)) for k in range(6)] for j in range(Q)]
c=[mp.mpf(1 if j%3==0 else -1 if j%3==1 else 0) for j in range(Q)];classes=[[j for j in range(Q) if j%3==a] for a in range(3)]
def soft(v):
    M=max(v);e=[mp.exp(x-M) for x in v];s=sum(e);return [x/s for x in e]
def py(y):
    h=[sum(B[j][k]*y[k] for k in range(6)) for j in range(Q)];S=[sum(mp.exp(h[j]) for j in classes[a]) for a in range(3)];z=mp.mpf('.5')*mp.log(S[1]/S[0]);return soft([h[j]+z*c[j] for j in range(Q)])
def Fy(*yy):
    p=py(yy);return tuple(yy[k]/g-sum(B[j][k]*p[j] for j in range(Q)) for k in range(6))
def jac(y):
    p=py(y);C=[[(p[i] if i==j else 0)-p[i]*p[j] for j in range(Q)] for i in range(Q)]
    def bil(a,b):return sum(a[i]*C[i][j]*b[j] for i in range(Q) for j in range(Q))
    u=[bil([B[i][a] for i in range(Q)],c) for a in range(6)];cc=bil(c,c);J=mp.matrix(6)
    for a in range(6):
      for b in range(6):J[a,b]=(1/g if a==b else 0)-(bil([B[i][a] for i in range(Q)],[B[i][b] for i in range(Q)])-u[a]*u[b]/cc)
    return J
seed=(mp.mpf('2.5703048563159232'),0,mp.mpf('1.2992516656905943'),mp.mpf('.6817926532921314'),mp.mpf('-1.1808995157291638'),mp.mpf('1.9759863787908111'))
y=list(mp.findroot(Fy,seed,tol=mp.mpf('1e-70'),maxsteps=100));J0=jac(y);R=J0**-1
# interval helper
def I(x):
    if isinstance(x,(tuple,list)):return iv.mpf([mp.nstr(x[0],100),mp.nstr(x[1],100)])
    return iv.mpf(mp.nstr(x,100))
def bd(x):
    s=mp.nstr(x,100)
    if s.startswith('['):return mp.mpf(s[1:s.index(',')]),mp.mpf(s[s.index(',')+1:s.rindex(']')].strip())
    v=mp.mpf(s);return v,v
Bi=[[I(x) for x in row] for row in B];ci=[I(x) for x in c]
def piv(Y):
    h=[sum((Bi[j][k]*Y[k] for k in range(6)),I(0)) for j in range(Q)];S=[sum((iv.exp(h[j]) for j in classes[a]),I(0)) for a in range(3)];z=I(mp.mpf('.5'))*iv.log(S[1]/S[0]);e=[iv.exp(h[j]+z*ci[j]) for j in range(Q)];den=sum(e,I(0));return [v/den for v in e]
def jiv(Y):
    p=piv(Y);C=[[(p[i] if i==j else I(0))-p[i]*p[j] for j in range(Q)] for i in range(Q)]
    def bil(a,b):return sum((a[i]*C[i][j]*b[j] for i in range(Q) for j in range(Q)),I(0))
    u=[bil([Bi[i][a] for i in range(Q)],ci) for a in range(6)];cc=bil(ci,ci);J=[[I(0) for _ in range(6)] for __ in range(6)]
    for a in range(6):
      for b in range(6):J[a][b]=(I(1)/I(g) if a==b else I(0))-(bil([Bi[i][a] for i in range(Q)],[Bi[i][b] for i in range(Q)])-u[a]*u[b]/cc)
    return J
# Krawczyk
rad=mp.mpf('1e-18');Y=[I((y[k]-rad,y[k]+rad)) for k in range(6)];JX=jiv(Y);F0=[I(Fy(*y)[k]) for k in range(6)];K0=[I(y[i])-sum((I(R[i,j])*F0[j] for j in range(6)),I(0)) for i in range(6)];E=[[(I(1) if i==j else I(0))-sum((I(R[i,k])*JX[k][j] for k in range(6)),I(0)) for j in range(6)] for i in range(6)];dx=[I((-rad,rad)) for _ in range(6)];KI=[K0[i]+sum((E[i][j]*dx[j] for j in range(6)),I(0)) for i in range(6)];strict=True
for i in range(6):xl,xu=bd(Y[i]);kl,ku=bd(KI[i]);strict &= (xl<kl and ku<xu)
# tighter local Hessian bound via p-box + LP as custom checker
hrad=mp.mpf('.02');lo=[v-hrad for v in y];hi=[v+hrad for v in y];hL=[];hU=[]
for j in range(Q):
 a=b=mp.mpf(0)
 for k in range(6):
  if B[j][k]>=0:a+=B[j][k]*lo[k];b+=B[j][k]*hi[k]
  else:a+=B[j][k]*hi[k];b+=B[j][k]*lo[k]
 hL.append(a);hU.append(b)
eL=[mp.exp(v) for v in hL];eU=[mp.exp(v) for v in hU];SL=[sum(eL[j] for j in C) for C in classes];SU=[sum(eU[j] for j in C) for C in classes];qL=mp.sqrt(SL[0]*SL[1]);qU=mp.sqrt(SU[0]*SU[1]);mL=qL/(2*qL+SU[2]);mU=qU/(2*qU+SL[2]);m2L=SL[2]/(2*qU+SL[2]);m2U=SU[2]/(2*qL+SU[2]);pL=[0]*Q;pU=[0]*Q
for a,C in enumerate(classes):
 cmL=mL if a<2 else m2L;cmU=mU if a<2 else m2U
 for i in C:
  pL[i]=cmL*eL[i]/(eL[i]+sum(eU[j] for j in C if j!=i));pU[i]=cmU*eU[i]/(eU[i]+sum(eL[j] for j in C if j!=i))
def ext(co,mx):
 xx=pL.copy();rem=mp.mpf(1)-sum(xx);order=sorted(range(Q),key=lambda i:co[i],reverse=mx)
 for i in order:
  if rem<=0:break
  add=min(rem,pU[i]-xx[i]);xx[i]+=add;rem-=add
 return sum(co[i]*xx[i] for i in range(Q))
def mul(U,V):
 vals=[U[0]*V[0],U[0]*V[1],U[1]*V[0],U[1]*V[1]];return min(vals),max(vals)
mu=[];uu=[]
for a in range(6):mu.append((ext([B[j][a] for j in range(Q)],False),ext([B[j][a] for j in range(Q)],True)));uu.append((ext([B[j][a]*c[j] for j in range(Q)],False),ext([B[j][a]*c[j] for j in range(Q)],True)))
cc=(2*mL,2*mU);HL=np.zeros((6,6));HU=np.zeros((6,6))
for a in range(6):
 for b in range(6):
  ebb=(ext([B[j][a]*B[j][b] for j in range(Q)],False),ext([B[j][a]*B[j][b] for j in range(Q)],True));mm=mul(mu[a],mu[b]);qq=mul(uu[a],uu[b]);vv=[qq[0]/cc[0],qq[0]/cc[1],qq[1]/cc[0],qq[1]/cc[1]];div=(min(vv),max(vv));klo=ebb[0]-mm[1]-div[1];khi=ebb[1]-mm[0]-div[0];HL[a,b]=float((1/g if a==b else 0)-khi);HU[a,b]=float((1/g if a==b else 0)-klo)
mid=(HL+HU)/2;rr=(HU-HL)/2;hesslb=float(np.linalg.eigvalsh((mid+mid.T)/2)[0]-np.linalg.norm(rr,'fro')-1e-10)
# solve full localized and saddle from simple seeds in theta
# d4 saddle reconstructed from y + normal z
def y_to_theta(y):
 h=[sum(B[j][k]*y[k] for k in range(6)) for j in range(Q)];S=[sum(mp.exp(h[j]) for j in classes[a]) for a in range(3)];z=mp.mpf('.5')*mp.log(S[1]/S[0]);
 # normal vector in k4 plane: n=(sqrt3/2,-1/2) orthogonal to Tc k4 tangent (.5,sqrt3/2)
 th=[mp.mpf(0)]*7
 for a in range(7):th[a]=sum(Tc[a][k]*y[k] for k in range(6))
 th[2]+=mp.sqrt(3)/2*z;th[3]+=-mp.mpf('.5')*z
 return th
sadth=y_to_theta(y)
def pth(th):return soft([sum(X[j][a]*th[a] for a in range(7)) for j in range(Q)])
def Ft(*tt):
 p=pth(tt);return tuple(tt[a]-g*sum(X[j][a]*p[j] for j in range(Q)) for a in range(7))
# localized start g X0
locseed=[g*X[0][a] for a in range(7)];locth=list(mp.findroot(Ft,tuple(locseed),tol=mp.mpf('1e-65'),maxsteps=100));ploc=pth(locth);psad=pth(sadth);ploc4=[ploc[(j-4)%12] for j in range(Q)]
def V(p):return sum(p[j]*mp.log(12*p[j]) for j in range(Q))-g*sum(p[i]*A[i][j]*p[j] for i in range(Q) for j in range(Q))/2
Vloc=V(ploc);Vsad=V(psad)
# path interval derivative
def pathcert(a,b,sign,nmid=2000,edge=mp.mpf('.002')):
 d=[b[i]-a[i] for i in range(Q)];dAd=sum(d[i]*A[i][j]*d[j] for i in range(Q) for j in range(Q))
 def d2(tI):
  p=[I(a[i])+tI*I(d[i]) for i in range(Q)];return sum((I(d[i]*d[i])/p[i] for i in range(Q)),I(0))-I(g*dAd)
 left=d2(I((0,edge)));right=d2(I((1-edge,1)));ll,lu=bd(left);rl,ru=bd(right);worst=mp.inf if sign>0 else -mp.inf;ok=True
 for rr0 in range(nmid):
  aa=edge+(1-2*edge)*rr0/nmid;bb=edge+(1-2*edge)*(rr0+1)/nmid;tI=I((aa,bb));p=[I(a[i])+tI*I(d[i]) for i in range(Q)];grad=[]
  for i in range(Q):grad.append(iv.log(I(12)*p[i])+I(1)-I(g)*sum((I(A[i][j])*p[j] for j in range(Q)),I(0)))
  dv=sum((I(d[i])*grad[i] for i in range(Q)),I(0));dl,du=bd(dv)
  if sign>0:worst=min(worst,dl);ok &= dl>0
  else:worst=max(worst,du);ok &= du<0
 edgeok=(ll>0 and ru<0) if sign>0 else (lu<0 and rl>0)
 return {'middle_strict':bool(ok),'endpoint_curvature_strict':bool(edgeok),'middle_worst_bound':str(worst),'left_d2':[str(ll),str(lu)],'right_d2':[str(rl),str(ru)]}
out={'exact_kernel':True,'root_y':[mp.nstr(v,50) for v in y],'krawczyk_strict':bool(strict),'local_hessian_eigen_lower_bound':hesslb,'localized_energy':mp.nstr(Vloc,60),'d4_energy':mp.nstr(Vsad,60),'B4':mp.nstr(Vsad-Vloc,60),'path_loc_to_d4':pathcert(ploc,psad,+1),'path_d4_to_shift4loc':pathcert(psad,ploc4,-1)}
print(json.dumps(out,indent=2));open('/mnt/data/fin326_replay/fin326/EXACT_KERNEL_LOCAL_PATH_327.json','w').write(json.dumps(out,indent=2)+'\n')
