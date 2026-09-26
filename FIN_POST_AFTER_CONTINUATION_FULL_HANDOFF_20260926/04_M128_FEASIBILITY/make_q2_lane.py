import sys,math,json
import numpy as np
from scipy import fft
base=sys.argv[1]; alpha=float(sys.argv[2]); pref=sys.argv[3]
V=np.array([
[4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343],
[-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338],
[0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338],
[5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343],
[-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343],
[0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338],
[-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338],
[5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338]])
G=V@V.T
def dec(h):
 x=int(h,16);return tuple((x>>(5*i))&31 for i in range(8))
def norm(c):
 a=np.asarray(c,float);return float(a@G@a)
rows=[]
with open(base) as f:
 for line in f:
  p=line.split();c=dec(p[0]);vv=[]
  for k in range(4):
   z=p[1+3*k:1+3*k+3];vv.append(None if z[0]=='nan' else tuple(map(float,z)))
  rows.append((c,vv))
shape=(9,)*7; N=9**7
# q2 M16 actual dense + derivative arrays
B0=np.zeros(shape);B1=np.zeros(shape);B2=np.zeros(shape); logs=[]; tmp=[]
a16=4*2/16**2
for c,v in rows:
 if max(c)<=8 and v[1] is not None:
  lz,d1,d2=v[1];n=norm(c);lg=lz-alpha*a16*n;logs.append(lg);tmp.append((c,n,lg,d1,d2))
M=max(logs)
A0=np.zeros(shape);A1=np.zeros(shape);A2=np.zeros(shape)
for c,n,lg,d1,d2 in tmp:
 z=math.exp(lg-M);l1=d1-a16*n
 A0[c[:7]]=z;A1[c[:7]]=z*l1;A2[c[:7]]=z*(d2+l1*l1)
 # actual base q2
 za=math.exp(lg); B0[c[:7]]=za;B1[c[:7]]=za*l1;B2[c[:7]]=za*(d2+l1*l1)
for arr,suf in [(B0,'q2m16_g0.raw'),(B1,'q2m16_g1.raw'),(B2,'q2m16_g2.raw')]:arr.astype('<f8').tofile(pref+suf)
F0=fft.rfftn(A0,workers=5);F1=fft.rfftn(A1,workers=5);F2=fft.rfftn(A2,workers=5)
C0=fft.irfftn(F0*F0,s=shape,workers=5);C1=fft.irfftn(2*F0*F1,s=shape,workers=5);C2=fft.irfftn(2*F0*F2+2*F1*F1,s=shape,workers=5)
# q4 diag target c=2h, scaled e^-2M
D0=np.zeros(shape);D1=np.zeros(shape);D2=np.zeros(shape);a4=4*4/16**2
for c,v in rows:
 if max(c)<=4 and v[2] is not None:
  lz,d1,d2=v[2];n=norm(c);lg=lz-alpha*a4*n;z=math.exp(lg-2*M);l1=d1-a4*n;tc=tuple(2*x for x in c)
  D0[tc[:7]]=z;D1[tc[:7]]=z*l1;D2[tc[:7]]=z*(d2+l1*l1)
# valid M32 cap8 vectorized
flat=np.arange(N,dtype=np.int64);t=flat.copy();dig=[None]*7;s=np.zeros(N,dtype=np.int16)
for i in range(6,-1,-1):
 d=(t%9).astype(np.int8);t//=9;dig[i]=d;s+=d
last=32-s;mask=(last>=0)&(last<=8);idx=flat[mask]
nv=np.empty(idx.size);ch=400000
for st in range(0,idx.size,ch):
 en=min(st+ch,idx.size);cs=np.column_stack([dig[i][idx[st:en]] for i in range(7)]+[last[idx[st:en]]]).astype(float);nv[st:en]=np.einsum('bi,ij,bj->b',cs,G,cs,optimize=True)
a32=4*2/32**2;T0=C0.ravel()[idx]+D0.ravel()[idx];T1=C1.ravel()[idx]+D1.ravel()[idx];T2=C2.ravel()[idx]+D2.ravel()[idx];fac=.5*np.exp(alpha*a32*nv)
P0=np.zeros(N);P1=np.zeros(N);P2=np.zeros(N)
p0=fac*T0;p1=fac*(T1+a32*nv*T0);p2=fac*(T2+2*a32*nv*T1+(a32*nv)**2*T0)
# undo q2 scale e^-2M -> actual
sc=math.exp(2*M);P0[idx]=p0*sc;P1[idx]=p1*sc;P2[idx]=p2*sc
for arr,suf in [(P0,'q2m32_g0.raw'),(P1,'q2m32_g1.raw'),(P2,'q2m32_g2.raw')]:arr.astype('<f8').tofile(pref+suf)
# q4 M32 homogeneous (4^8) from q4 base + q8 diag
q4={};q8={};a4b=4*4/16**2;a8=4*8/16**2
for c,v in rows:
 if v[2] is not None:
  lz,d1,d2=v[2];n=norm(c);lg=lz-alpha*a4b*n;z=math.exp(lg);l1=d1-a4b*n;q4[c]=(z,z*l1,z*(d2+l1*l1))
 if v[3] is not None:
  lz,d1,d2=v[3];n=norm(c);lg=lz-alpha*a8*n;z=math.exp(lg);l1=d1-a8*n;q8[c]=(z,z*l1,z*(d2+l1*l1))
def parent_hom(child,diag,rootval,halfeq,q,parent_m):
 o0=o1=o2=0.;cnt=0
 for c,v in child.items():
  if max(c)>rootval:continue
  b=tuple(rootval-x for x in c)
  if b not in child:continue
  w=child[b];o0+=v[0]*w[0];o1+=v[1]*w[0]+v[0]*w[1];o2+=v[2]*w[0]+2*v[1]*w[1]+v[0]*w[2];cnt+=1
 e=diag[halfeq];T0=o0+e[0];T1=o1+e[1];T2=o2+e[2];root=(rootval,)*8;n=norm(root);a=4*q/parent_m**2;f=.5*math.exp(alpha*a*n)
 return (f*T0,f*(T1+a*n*T0),f*(T2+2*a*n*T1+(a*n)**2*T0),cnt)
Q4_32=parent_hom(q4,q8,4,(2,)*8,4,32)
# q2 M64 homogeneous root (8^8) using dense q2m32 actual arrays + q4M32 equal
# sum over all 9^7 grid and reverse (complements 8-digits); invalid zeros harmless
r0=P0.reshape(shape);r1=P1.reshape(shape);r2=P2.reshape(shape);rev=(slice(None,None,-1),)*7
o0=float(np.sum(r0*r0[rev]));o1=float(np.sum(r1*r0[rev]+r0*r1[rev]));o2=float(np.sum(r2*r0[rev]+2*r1*r1[rev]+r0*r2[rev]))
e=Q4_32[:3];T0=o0+e[0];T1=o1+e[1];T2=o2+e[2];n=norm((8,)*8);a=4*2/64**2;f=.5*math.exp(alpha*a*n);Q2_64=(f*T0,f*(T1+a*n*T0),f*(T2+2*a*n*T1+(a*n)**2*T0))
meta={'alpha':alpha,'M_q2_16':M,'q4m32_hom':Q4_32[:3],'q2m64_hom':Q2_64,'q2m32_nonzero':int(idx.size)}
open(pref+'q2meta.json','w').write(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2))
