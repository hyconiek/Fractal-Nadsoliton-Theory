import mpmath as mp, math, json
mp.mp.dps=90; mp.iv.dps=70; iv=mp.iv
Q=12;D=6;Gm=mp.mpf('5.145228719489142'); G=iv.mpf([mp.nstr(Gm,90),mp.nstr(Gm,90)])
# exact high-precision mathematical kernel / feature construction
W=[[mp.mpf(0) if i==j else mp.cos(mp.mpf('0.18575')*min(abs(i-j),Q-abs(i-j))+mp.mpf('0.1625'))/(1+mp.mpf(min(abs(i-j),Q-abs(i-j)))**mp.mpf('1.8')) for j in range(Q)] for i in range(Q)]
r=sum(W[0]); l0=[(r if j==0 else mp.mpf(0))-W[0][j] for j in range(Q)]; lam=[sum(l0[j]*mp.cos(2*mp.pi*k*j/Q) for j in range(Q)) for k in range(7)]
X=[]
for j in range(Q):
 row=[]
 for k in (3,4,5):row += [mp.sqrt(lam[k]/6)*mp.cos(2*mp.pi*k*j/Q),mp.sqrt(lam[k]/6)*mp.sin(2*mp.pi*k*j/Q)]
 row += [mp.sqrt(lam[6]/12)*((-1)**j)];X.append(row)
Tc=[[mp.mpf(0) for _ in range(6)] for __ in range(7)];Tc[0][0]=1;Tc[1][1]=1;Tc[2][2]=mp.mpf('.5');Tc[3][2]=mp.sqrt(3)/2;Tc[4][3]=1;Tc[5][4]=1;Tc[6][5]=1
Bm=[[sum(X[j][a]*Tc[a][k] for a in range(7)) for k in range(6)] for j in range(Q)]; B=[[iv.mpf([mp.nstr(x,90),mp.nstr(x,90)]) for x in row] for row in Bm]
cls=lambda j:j%3
c=[iv.mpf(1 if cls(j)==0 else -1 if cls(j)==1 else 0) for j in range(Q)]
def bd(x):
 s=mp.nstr(x,100)
 if s.startswith('['):return mp.mpf(s[1:s.index(',')]),mp.mpf(s[s.index(',')+1:s.rindex(']')].strip())
 v=mp.mpf(s);return v,v

def lin_ext(co,pL,pU,mx):
 cc=[bd(x)[1 if mx else 0] for x in co]; L=[bd(x)[0] for x in pL];U=[bd(x)[1] for x in pU]
 x=L[:];rem=mp.mpf(1)-sum(x);order=sorted(range(Q),key=lambda i:cc[i],reverse=mx)
 for i in order:
  if rem<=0:break
  add=min(rem,U[i]-x[i]);x[i]+=add;rem-=add
 return sum(cc[i]*x[i] for i in range(Q))
def pbounds(box):
 h=[];e=[];S=[iv.mpf(0),iv.mpf(0),iv.mpf(0)]
 for j in range(Q):
  z=sum((B[j][k]*iv.mpf([mp.nstr(box[k][0],90),mp.nstr(box[k][1],90)]) for k in range(D)),iv.mpf(0));h.append(z);e.append(iv.exp(z));S[cls(j)]+=e[-1]
 q=iv.sqrt(S[0]*S[1]);den=2*q+S[2];m=q/den;m2=S[2]/den;p=[None]*Q
 for a in range(3):
  mass=m if a<2 else m2
  for j in range(Q):
   if cls(j)==a:p[j]=mass*e[j]/S[a]
 return p,m

def phi_grad(cen):
 e=[];S=[iv.mpf(0),iv.mpf(0),iv.mpf(0)];h=[]
 for j in range(Q):
  z=sum((B[j][k]*iv.mpf(mp.nstr(cen[k],90)) for k in range(D)),iv.mpf(0));h.append(z);e.append(iv.exp(z));S[cls(j)]+=e[-1]
 z=iv.mpf('.5')*iv.log(S[1]/S[0]);ef=[];den=iv.mpf(0)
 for j in range(Q):ef.append(iv.exp(h[j]+z*c[j]));den+=ef[-1]
 p=[x/den for x in ef]; norm=sum((iv.mpf(mp.nstr(x,90))**2 for x in cen),iv.mpf(0));ph=norm/(2*G)-iv.log((2*iv.sqrt(S[0]*S[1])+S[2])/12)
 gr=[]
 for k in range(D):gr.append(iv.mpf(mp.nstr(cen[k],90))/G-sum((p[j]*B[j][k] for j in range(Q)),iv.mpf(0)))
 return ph,gr

def m_lower(pL,pU):
 normco=[]
 for j in range(Q):normco.append(sum((B[j][k]*B[j][k] for k in range(D)),iv.mpf(0)))
 enn=lin_ext(normco,pL,pU,True);minsq=mp.mpf(0)
 for k in range(D):
  co=[B[j][k] for j in range(Q)];lo=lin_ext(co,pL,pU,False);hi=lin_ext(co,pL,pU,True)
  if lo>0:minsq+=lo*lo
  elif hi<0:minsq+=hi*hi
 return bd(1/G-iv.mpf(enn-minsq))[0]
def qmin(gi,w,m):
 gl,gu=bd(gi); vals=[mp.mpf(0)]
 def ev(g,d):return g*d+mp.mpf('.5')*m*d*d
 vals += [ev(gu,-w),ev(gl,w)]
 if m>0:
  d=-gu/m
  if -w<=d<=0:vals.append(ev(gu,d))
  d=-gl/m
  if 0<=d<=w:vals.append(ev(gl,d))
 return min(vals)
def lower(box):
 p,mass=pbounds(box);cen=[(a+b)/2 for a,b in box];w=[(b-a)/2 for a,b in box];ph,gr=phi_grad(cen);ml=m_lower(p,p); # p same intervals passed both
 return bd(ph)[0]+sum(qmin(gr[k],w[k],ml) for k in range(D))
def fp_k1_gap(box):
 p,mass=pbounds(box);co=[B[j][1] for j in range(Q)];lo=lin_ext(co,p,p,False);hi=lin_ext(co,p,p,True);t=G*iv.mpf([mp.nstr(lo,90),mp.nstr(hi,90)]);tl,tu=bd(t)
 # reason 11 means x.u[1] < lower target
 return tl-box[1][1],(lo,hi),(tl,tu),bd(mass)
FB=[(mp.mpf('2.5281079863900867187'),mp.mpf('2.5740735861426337499')),(-mp.mpf('0.045965599752547031249'),-mp.mpf('0.022982799876273515625')),(mp.mpf('1.265584694959951875'),mp.mpf('1.2838383203680280079')),(mp.mpf('0.67176155908084078124'),mp.mpf('0.69664161682457562499')),(-mp.mpf('1.2440028871867421875'),-mp.mpf('1.2191228294430073438')),(mp.mpf('1.9889892617352887501'),mp.mpf('2.0067480944293538282'))]
EB=[(mp.mpf('2.5855649860807705079'),mp.mpf('2.5970563860189072656')),(-mp.mpf('0.011491399938136757812'),mp.mpf('0')),(mp.mpf('1.2975285394240851074'),mp.mpf('1.3020919457761041407')),(mp.mpf('0.67176155908084078124'),mp.mpf('0.68420158795270820309')),(-mp.mpf('1.1755827283914713671'),-mp.mpf('1.1693627139555376562')),(mp.mpf('1.9623510126941911329'),mp.mpf('1.9667907208677074024'))]
out={'critical_feasibility':{},'critical_energy':{}}
gap,mu,t,mass=fp_k1_gap(FB);out['critical_feasibility']={'gap_exact_kernel':str(gap),'mu1_range':[str(mu[0]),str(mu[1])],'target_y1_range':[str(t[0]),str(t[1])],'class_mass_range':[str(mass[0]),str(mass[1])]}
lb=lower(EB);target=mp.mpf('-0.1431003250148400');out['critical_energy']={'lower_bound_exact_kernel':str(lb),'target':str(target),'margin':str(lb-target)}
print(json.dumps(out,indent=2));open('/mnt/data/fin326_replay/fin326/CRITICAL_MP_INTERVAL_326.json','w').write(json.dumps(out,indent=2)+'\n')
