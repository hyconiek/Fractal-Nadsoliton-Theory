"""Training outputs only. Reserved times and delta_2 are NEVER evaluated here."""
import sys,time
from model import *
TIMES=(0.15,0.5,1.,1.5)
PREPS=prep_names()

def diagnostic_twotime(m,T,B,p,t1,t2):
    U=m['U'];sqrtpiY=np.sqrt(m['piY'])
    z0=z_from_p(m,p)
    w0=np.r_[np.asarray(U.T@z0),B.T@z0]
    prop1=expm(T*t1);prop2=expm(T*(t2-t1));w1=prop1@w0
    # Observed y=(0,0) vs complementary unclassified event, projector at microscopic t1
    VE=np.column_stack((U[:12,:].toarray(),B[:12,:]));wE=VE.T@(VE@w1)
    a= sqrtpiY*(prop2@wE)[:Y];b=sqrtpiY*(prop2@(w1-wE))[:Y]
    return np.r_[a,b]

def train_one(name,A,g):
    m=build(A,g); checks=model_diagnostics(m)
    refs={p:{str(t):observed(m,uniform_tensor(m,make_prep(m,p),t)) for t in TIMES} for p in PREPS}
    errs={}
    for order in (1,2,4,6):
        B,labels=orthogonal_basis(m,order,True);T,w_init,w_to_y,value=reduced(m,B)
        rows=[];nonneg=[];maxsum=[]
        for prep in PREPS:
          for t in TIMES:
            pred=value(prep,t);ref=refs[prep][str(t)]
            rows.append({'preparation':prep,'t':t,'TV_vs_parent':tv(ref,pred),'TV_markov':tv(ref,markov(m,t)),'min_prediction':float(pred.min()),'sum_prediction':float(pred.sum())})
        errs[str(order)]={'rank':int(B.shape[1]),'basis_source_vectors':labels,'max_train_TV':float(max(r['TV_vs_parent'] for r in rows)),'min_probability':float(min(r['min_prediction'] for r in rows)),'max_mass_error':float(max(abs(r['sum_prediction']-1) for r in rows)),'rows':rows,'eig_max':float(np.linalg.eigvalsh(T)[-1]),'K0':memory_operator_stats(m,B)['k0_00']}
    # train two-time selected readout only, 0.15->0.5, Y00 versus complement
    p0=make_prep(m,'conditional_00');pb=uniform_tensor(m,p0,0.15)
    mask=np.zeros(N);mask[:12]=1
    refJ=np.r_[observed(m,uniform_tensor(m,pb*mask,.35)), observed(m,uniform_tensor(m,pb*(1-mask),.35))]
    tw={}
    for order in (1,2,4,6):
        B,_=orthogonal_basis(m,order,True);T,_,_,_=reduced(m,B)
        jp=diagnostic_twotime(m,T,B,p0,.15,.5)
        tw[str(order)]={'TV':tv(refJ,jp),'mass_error':abs(float(jp.sum())-1),'min_joint':float(jp.min())}
    return m,{'name':name,'g':g,'checks':checks,'training_by_order':errs,'train_twotime':tw,'training_reference_cond00':{str(t):ref[0] for t,ref in []},'reference_normalization_max':float(max(abs(refs[p][str(t)].sum()-1) for p in PREPS for t in TIMES)),'piY_event00':float(m['piY'][0]),'P0_'+name:float(m['P'].data.min())},refs

def main():
    t=time.time();mods=models();out={};refs={}
    for name,A in mods.items():
        m,res,refl=train_one(name,A,G);out[name]=res
        if name=='FIN_A7': refs={p:{tt:v.tolist() for tt,v in entries.items()} for p,entries in refl.items()}
        print('trained',name,'time',time.time()-t,'order6 max',res['training_by_order']['6']['max_train_TV'],flush=True)
    m0,res,_=train_one('FIN_A7_g0',mods['FIN_A7'],0.);out['FIN_A7_g0']=res
    write('RG-D03/TRAINING_PREDICTIONS.json',out)
    write('RG-D02/TRAINING_REFERENCE.json',{'model_id':'FIN_A7','training_times':list(TIMES),'training_preparations':list(PREPS),'all_144_outcomes':refs,'arithmetic_scope':'binary64 uniformization; all-operation allowance separately in ERROR_CERTIFICATE.json'})
    print('training done seconds',time.time()-t,flush=True)
if __name__=='__main__':main()
