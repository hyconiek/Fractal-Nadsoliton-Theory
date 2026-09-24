#!/usr/bin/env python3
import numpy as np, json
T=(2*np.pi)**2
rows=[]
for a in [1,2,3,4]:
 for m in [16,32,64,128,256]:
  nx=a*m;ny=m;hx=1/nx;hy=1/ny;N=nx*ny
  dx=2*(1-np.cos(2*np.pi/nx));dy=2*(1-np.cos(2*np.pi/ny))
  rawx=N*dx;rawy=N*dy;mcell=hx*hy;cx=hy/hx;cy=hx/hy
  wx=cx*dx/mcell;wy=cy*dy/mcell
  rows.append((a,m,rawy/rawx,wy/wx,max(abs(wx-T),abs(wy-T))))
slopes={}
for a in [1,2,3,4]:
 r=[x for x in rows if x[0]==a]; slopes[str(a)]=float(np.polyfit(np.log([x[1] for x in r]),np.log([x[4] for x in r]),1)[0])
print(json.dumps({'last_rows':[{'a':a,'m':m,'raw_ratio':float(rr),'weighted_ratio':float(wr),'error':float(e)} for a,m,rr,wr,e in rows if m==256],'slopes':slopes},indent=2,sort_keys=True))
