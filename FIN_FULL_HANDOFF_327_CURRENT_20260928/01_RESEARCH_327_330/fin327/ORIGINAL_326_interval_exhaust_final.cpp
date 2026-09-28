#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cstdint>
using ld=long double;
constexpr int Q=12,D=6;
const ld G=5.145228719489142L;
// B matrix generated from exact same spectral construction as report 325, printed at long-double input precision.
ld B[Q][D] = {
  {0.57175269449539189992L,0L,0.30273536637722164855L,0.6189515694676895885L,0L,0.44179388493073518118L},
  {3.500975536088280713e-17L,0.57175269449539189992L,0.30273536637722109344L,-0.53602778287126795487L,0.30947578473384473874L,-0.44179388493073518118L},
  {-0.57175269449539189992L,7.0019510721765614259e-17L,-0.60547073275444263096L,0.30947578473384484976L,-0.53602778287126784384L,0.44179388493073518118L},
  {-1.0502926608264842139e-16L,-0.57175269449539189992L,0.30273536637722153753L,7.392386914346839704e-16L,0.6189515694676895885L,-0.44179388493073518118L},
  {0.57175269449539189992L,-1.4003902144353122852e-16L,0.3027353663772214265L,-0.30947578473384468323L,-0.53602778287126795487L,0.44179388493073518118L},
  {6.8286718146061965487e-16L,0.57175269449539189992L,-0.60547073275444263096L,0.5360277828712676218L,0.30947578473384512732L,-0.44179388493073518118L},
  {-0.57175269449539189992L,2.1005853216529684278e-16L,0.30273536637722137099L,-0.6189515694676895885L,1.4784773828693679408e-15L,0.44179388493073518118L},
  {-2.4506828752617964991e-16L,-0.57175269449539189992L,0.30273536637722164855L,0.53602778287126784384L,-0.30947578473384484976L,-0.44179388493073518118L},
  {0.57175269449539189992L,-2.8007804288706245704e-16L,-0.60547073275444263096L,-0.30947578473384501629L,0.53602778287126784384L,0.44179388493073518118L},
  {3.1508779824794531347e-16L,0.57175269449539189992L,0.30273536637722125997L,5.3098105989957209857e-16L,-0.6189515694676895885L,-0.44179388493073518118L},
  {-0.57175269449539189992L,1.3657343629212393097e-15L,0.30273536637722081588L,0.3094757847338440726L,0.53602778287126839896L,0.44179388493073518118L},
  {6.3052950034270040933e-16L,-0.57175269449539189992L,-0.60547073275444263096L,-0.53602778287126728873L,-0.30947578473384573794L,-0.44179388493073518118L}
};
struct Box{std::array<ld,D> l,u; uint16_t depth;};
ld initL[D],initU[D];
ld roots4[4][D];
ld targetU=-0.1431003250148460L; // deliberately above all four numerical root values
const ld LOCAL_R=0.02L;
const ld SAFE=1e-6L;
int cls(int j){return j%3;}
ld cval(int j){return cls(j)==0?1.0L:(cls(j)==1?-1.0L:0.0L);}

void p_bounds(const Box& x, ld pL[Q], ld pU[Q], ld &mL, ld &mU){
    ld hL[Q],hU[Q],eL[Q],eU[Q];
    for(int j=0;j<Q;j++){
        ld a=0,b=0;
        for(int k=0;k<D;k++){
            ld coef=B[j][k];
            if(coef>=0){a+=coef*x.l[k];b+=coef*x.u[k];}
            else {a+=coef*x.u[k];b+=coef*x.l[k];}
        }
        hL[j]=a;hU[j]=b;eL[j]=expl(a);eU[j]=expl(b);
    }
    ld SL[3]={0,0,0},SU[3]={0,0,0};
    for(int j=0;j<Q;j++){SL[cls(j)]+=eL[j];SU[cls(j)]+=eU[j];}
    ld qL=sqrtl(SL[0]*SL[1]), qU=sqrtl(SU[0]*SU[1]);
    mL=qL/(2*qL+SU[2]); mU=qU/(2*qU+SL[2]);
    ld m2L=SL[2]/(2*qU+SL[2]),m2U=SU[2]/(2*qL+SU[2]);
    for(int a=0;a<3;a++){
        ld sumL=SL[a],sumU=SU[a];
        ld cmL=(a<2?mL:m2L),cmU=(a<2?mU:m2U);
        for(int j=0;j<Q;j++) if(cls(j)==a){
            ld wL=eL[j]/(eL[j]+sumU-eU[j]);
            ld wU=eU[j]/(eU[j]+sumL-eL[j]);
            pL[j]=cmL*wL; pU[j]=cmU*wU;
            pL[j]=std::max((ld)0,pL[j]-1e-18L); pU[j]=std::min((ld)1,pU[j]+1e-18L);
        }
    }
}
ld lin_ext(const ld coeff[Q],const ld pL[Q],const ld pU[Q],bool mx){
    ld x[Q],rem=1; int idx[Q];
    for(int i=0;i<Q;i++){x[i]=pL[i];rem-=x[i];idx[i]=i;}
    std::sort(idx,idx+Q,[&](int a,int b){return mx?coeff[a]>coeff[b]:coeff[a]<coeff[b];});
    for(int t=0;t<Q && rem>0;t++){
        int i=idx[t]; ld add=std::min(rem,pU[i]-x[i]); if(add>0){x[i]+=add;rem-=add;}
    }
    ld s=0;for(int i=0;i<Q;i++)s+=coeff[i]*x[i];return s;
}
void phi_grad(const std::array<ld,D>& y,ld &phi,std::array<ld,D>& gr){
    ld h[Q], e[Q], S[3]={0,0,0};
    for(int j=0;j<Q;j++){h[j]=0;for(int k=0;k<D;k++)h[j]+=B[j][k]*y[k]; e[j]=expl(h[j]);S[cls(j)]+=e[j];}
    ld z=.5L*logl(S[1]/S[0]); ld ef[Q],den=0;
    for(int j=0;j<Q;j++){ef[j]=expl(h[j]+z*cval(j));den+=ef[j];}
    ld p[Q];for(int j=0;j<Q;j++)p[j]=ef[j]/den;
    ld norm=0;for(int k=0;k<D;k++)norm+=y[k]*y[k];
    ld logz=logl((2*sqrtl(S[0]*S[1])+S[2])/12.0L);phi=norm/(2*G)-logz;
    for(int k=0;k<D;k++){ld mu=0;for(int j=0;j<Q;j++)mu+=p[j]*B[j][k];gr[k]=y[k]/G-mu;}
}
ld local_m_bound(const Box& x){
    ld pL[Q],pU[Q],mL,mU;p_bounds(x,pL,pU,mL,mU);
    ld normcoef[Q];for(int j=0;j<Q;j++){normcoef[j]=0;for(int k=0;k<D;k++)normcoef[j]+=B[j][k]*B[j][k];}
    ld enn=lin_ext(normcoef,pL,pU,true);
    ld minsq=0;
    for(int k=0;k<D;k++){
        ld co[Q];for(int j=0;j<Q;j++)co[j]=B[j][k];
        ld a=lin_ext(co,pL,pU,false),b=lin_ext(co,pL,pU,true);
        if(a>0)minsq+=a*a;else if(b<0)minsq+=b*b;
    }
    ld tr=std::max((ld)0,enn-minsq);
    return 1.0L/G-tr;
}

bool root_possible(const Box& x){
    ld pL[Q],pU[Q],mL,mU;p_bounds(x,pL,pU,mL,mU);
    // Representative half-wall P0=P1>=P2 is exactly m=P0=P1 >= 1/3.
    if(mU < 1.0L/3.0L - 1e-14L) return false;
    // Quotient exact internal C4 translation symmetry by the k=3 wedge y0 >= |y1|.
    // Every C4 orbit intersects this closed wedge; translation by 3 preserves the half-wall.
    if(x.u[0]-x.l[1] < -1e-14L) return false; // cannot satisfy y0>=y1
    if(x.u[0]+x.u[1] < -1e-14L) return false; // cannot satisfy y0>=-y1
    for(int k=0;k<D;k++){
        ld co[Q];for(int j=0;j<Q;j++)co[j]=B[j][k];
        ld tlo=G*lin_ext(co,pL,pU,false);
        ld thi=G*lin_ext(co,pL,pU,true);
        // generous outward padding
        tlo-=1e-12L*(1+fabsl(tlo)); thi+=1e-12L*(1+fabsl(thi));
        if(x.u[k]<tlo || x.l[k]>thi) return false;
    }
    return true;
}

ld lower_bound(const Box& x){
    std::array<ld,D> cen,w,gr;for(int k=0;k<D;k++){cen[k]=(x.l[k]+x.u[k])/2;w[k]=(x.u[k]-x.l[k])/2;}
    ld val;phi_grad(cen,val,gr);ld m=local_m_bound(x),s=val;
    for(int k=0;k<D;k++){
        ld d;
        if(m>0){d=std::max(-w[k],std::min(w[k],-gr[k]/m));}
        else d=(gr[k]>=0?-w[k]:w[k]);
        s+=gr[k]*d+.5L*m*d*d;
    }
    return s-1e-8L; // very conservative outward safety inflation
}
bool inside_local(const Box& x){
    for(int r=0;r<4;r++){
        bool ok=true;for(int k=0;k<D;k++) if(x.l[k]<roots4[r][k]-LOCAL_R || x.u[k]>roots4[r][k]+LOCAL_R){ok=false;break;}
        if(ok)return true;
    }
    return false;
}
int main(){
    // generated constants
    initL[0]=-2.9417983841629919972L; initU[0]=2.9417983841629919972L;
    initL[1]=-2.9417983841629919972L; initU[1]=2.9417983841629919972L;
    initL[2]=-3.115285402978293483L; initU[2]=1.5576427014891485179L;
    initL[3]=-3.1846473911980353044L; initU[3]=3.1846473911980353044L;
    initL[4]=-3.1846473911980353044L; initU[4]=3.1846473911980353044L;
    initL[5]=-2.2731305848403002834L; initU[5]=2.2731305848403002834L;
    roots4[0][0]=2.5703048563159236473L;
    roots4[0][1]=-2.7061686225238190673e-16L;
    roots4[0][2]=1.2992516656905954697L;
    roots4[0][3]=0.68179265329213267766L;
    roots4[0][4]=-1.1808995157291635181L;
    roots4[0][5]=1.9759863787908105159L;
    roots4[1][0]=-5.898059818321144121e-16L;
    roots4[1][1]=-2.5703048563159227591L;
    roots4[1][2]=1.2992516656905956918L;
    roots4[1][3]=1.1808995157291650724L;
    roots4[1][4]=0.68179265329213201152L;
    roots4[1][5]=-1.975986378790810738L;
    roots4[2][0]=-2.5703048563159232032L;
    roots4[2][1]=3.4937330806172894881e-15L;
    roots4[2][2]=1.2992516656905934713L;
    roots4[2][3]=-0.6817926532921338989L;
    roots4[2][4]=1.1808995157291677369L;
    roots4[2][5]=1.9759863787908105159L;
    roots4[3][0]=9.1072982488782372457e-16L;
    roots4[3][1]=2.5703048563159232032L;
    roots4[3][2]=1.2992516656905939154L;
    roots4[3][3]=-1.1808995157291624079L;
    roots4[3][4]=-0.68179265329213234459L;
    roots4[3][5]=-1.975986378790810738L;
    Box init;for(int k=0;k<D;k++){init.l[k]=initL[k];init.u[k]=initU[k];}init.depth=0;
    std::vector<Box> st;st.reserve(1000000);st.push_back(init);
    uint64_t nodes=0,discard=0,local=0,discardFeas=0,discardEnergy=0; ld minMargin=1e100L; int maxDepth=0; ld minWidth=1e100L;
    const uint64_t MAXN=200000000ULL;
    while(!st.empty() && nodes<MAXN){
        Box x=st.back();st.pop_back();nodes++;maxDepth=std::max(maxDepth,(int)x.depth);
        if(inside_local(x)){local++;continue;}
        if(!root_possible(x)){discard++;discardFeas++;continue;}
        ld lb=lower_bound(x);ld margin=lb-targetU;
        if(margin>SAFE){discard++;discardEnergy++;minMargin=std::min(minMargin,margin);continue;}
        int k=0;ld best=-1;
        for(int j=0;j<D;j++){ld score=(x.u[j]-x.l[j])/(initU[j]-initL[j]);if(score>best){best=score;k=j;}}
        ld width=x.u[k]-x.l[k];minWidth=std::min(minWidth,width);
        if(x.depth>120 || width<1e-7L){
            std::cerr<<"UNRESOLVED depth="<<x.depth<<" lb="<<(double)lb<<" width="<<(double)width<<" center=";
            for(int j=0;j<D;j++)std::cerr<<(double)((x.l[j]+x.u[j])/2)<<",";
            std::cerr<<"\n";return 3;
        }
        ld mid=(x.l[k]+x.u[k])/2;Box a=x,b=x;a.u[k]=mid;b.l[k]=mid;a.depth=b.depth=x.depth+1;st.push_back(b);st.push_back(a);
        if(nodes%250000==0)std::cerr<<"nodes "<<nodes<<" stack "<<st.size()<<" discard "<<discard<<" local "<<local<<" depth "<<maxDepth<<"\n";
    }
    std::cout<<std::setprecision(18);
    std::cout<<"nodes "<<nodes<<"\ndiscard "<<discard<<"\ndiscardFeas "<<discardFeas<<"\ndiscardEnergy "<<discardEnergy<<"\nlocal "<<local<<"\nremaining "<<st.size()<<"\nmaxDepth "<<maxDepth<<"\nminDiscardMargin "<<(double)minMargin<<"\n";
    return st.empty()?0:2;
}
