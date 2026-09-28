#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/transc.hpp>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstdint>
#include <limits>
namespace bn=boost::numeric; namespace il=boost::numeric::interval_lib;
using ld=long double;
using I=bn::interval<ld, il::policies<il::save_state<il::rounded_transc_std<ld>>, il::checking_base<ld>>>;
constexpr int Q=12,D=6;
const I G=I(5.145228719489141L,5.145228719489143L); // encloses declared decimal gain
const ld GLO=5.145228719489141L, GHI=5.145228719489143L;
// Each printed B coefficient is widened by 5e-18 to cover decimal serialization/long-double conversion.
const ld BR=2e-15L;
ld B0[Q][D] = {
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
I B[Q][D];
struct Box{std::array<ld,D> l,u; uint16_t depth;};
ld initL[D],initU[D], roots4[4][D];
const ld TARGET=-0.1431003250148400L; // strictly above certified local minimum value
const ld LOCAL_R=0.02L;
int cls(int j){return j%3;}
I cval(int j){return I(cls(j)==0?1.0L:(cls(j)==1?-1.0L:0.0L));}
I dot_row_box(int j,const Box& x){I s(0);for(int k=0;k<D;k++)s += B[j][k]*I(x.l[k],x.u[k]);return s;}

struct PB { I p[Q]; I m; };
PB p_bounds(const Box& x){
    I e[Q],S[3]={I(0),I(0),I(0)};
    for(int j=0;j<Q;j++){e[j]=exp(dot_row_box(j,x));S[cls(j)]+=e[j];}
    I q=sqrt(S[0]*S[1]); I den=I(2)*q+S[2];
    I m=q/den, m2=S[2]/den;
    PB out; out.m=m;
    // Tighter conditional-within-class bounds using interval softmax shares.
    for(int a=0;a<3;a++){
      I mass=(a<2?m:m2);
      for(int j=0;j<Q;j++) if(cls(j)==a){
        I rest=S[a]-e[j]; // dependency-safe inclusion; may be wider but rigorous
        I share=e[j]/(e[j]+rest);
        out.p[j]=mass*share;
        out.p[j]=intersect(out.p[j],I(0,1));
      }
    }
    return out;
}

// Rigorous extrema for linear form under p_i in [L_i,U_i], sum p=1.
// Coeff interval uncertainty is absorbed by using lower endpoints for min, upper for max.
ld lin_ext(const std::array<I,Q>& co,const PB& pb,bool mx){
    struct Item{ld c,L,U;}; Item a[Q]; ld rem=1;
    for(int i=0;i<Q;i++){a[i].c=mx?co[i].upper():co[i].lower();a[i].L=pb.p[i].lower();a[i].U=pb.p[i].upper();rem-=a[i].L;}
    std::sort(a,a+Q,[&](const Item&x,const Item&y){return mx?x.c>y.c:x.c<y.c;});
    I sum(0);
    for(int i=0;i<Q;i++) sum += I(a[i].c)*I(a[i].L);
    for(int i=0;i<Q && rem>0;i++){
      ld cap=std::max((ld)0,a[i].U-a[i].L), add=std::min(rem,cap);
      if(add>0){sum += I(a[i].c)*I(add);rem-=add;}
    }
    // If interval bounds are so loose that lower masses exceed 1 or uppers do not reach 1,
    // return a maximally conservative endpoint.
    if(rem>1e-18L) return mx?std::numeric_limits<ld>::infinity():-std::numeric_limits<ld>::infinity();
    return mx?sum.upper():sum.lower();
}

I phi_center(const std::array<ld,D>& y,std::array<I,D>& gr){
    I e[Q],S[3]={I(0),I(0),I(0)};
    for(int j=0;j<Q;j++){I h(0);for(int k=0;k<D;k++)h += B[j][k]*I(y[k]);e[j]=exp(h);S[cls(j)]+=e[j];}
    I z=I(0.5L)*log(S[1]/S[0]); I ef[Q],den(0);
    for(int j=0;j<Q;j++){I h(0);for(int k=0;k<D;k++)h += B[j][k]*I(y[k]);ef[j]=exp(h+z*cval(j));den+=ef[j];}
    I p[Q];for(int j=0;j<Q;j++)p[j]=ef[j]/den;
    I norm(0);for(int k=0;k<D;k++)norm+=I(y[k])*I(y[k]);
    I logz=log((I(2)*sqrt(S[0]*S[1])+S[2])/I(12));
    I phi=norm/(I(2)*G)-logz;
    for(int k=0;k<D;k++){I mu(0);for(int j=0;j<Q;j++)mu+=p[j]*B[j][k];gr[k]=I(y[k])/G-mu;}
    return phi;
}

ld m_lower(const PB& pb){
    std::array<I,Q> normco;
    for(int j=0;j<Q;j++){I s(0);for(int k=0;k<D;k++)s += B[j][k]*B[j][k];normco[j]=s;}
    ld enn=lin_ext(normco,pb,true);
    if(!std::isfinite((double)enn)) return -1e100L;
    I minsq(0);
    for(int k=0;k<D;k++){
      std::array<I,Q> co;for(int j=0;j<Q;j++)co[j]=B[j][k];
      ld lo=lin_ext(co,pb,false),hi=lin_ext(co,pb,true);
      if(lo>0)minsq += I(lo)*I(lo); else if(hi<0)minsq += I(hi)*I(hi);
    }
    I tr=I(enn)-minsq; if(tr.lower()<0) tr=I(0,tr.upper());
    I mm=I(1)/G-tr;
    return mm.lower();
}

ld quad_min_1d(const I& gi,ld w,ld m){
    ld gl=gi.lower(), gu=gi.upper();
    auto eval=[&](ld g,ld d){I v=I(g)*I(d)+I(0.5L)*I(m)*I(d)*I(d);return v.lower();};
    ld best=0; // d=0
    auto upd=[&](ld v){if(v<best)best=v;};
    // negative half: minimizing gradient endpoint gu
    upd(eval(gu,-w));
    if(m>0){ld d=-gu/m;if(d>=-w && d<=0)upd(eval(gu,d));}
    // positive half: minimizing gradient endpoint gl
    upd(eval(gl,w));
    if(m>0){ld d=-gl/m;if(d>=0 && d<=w)upd(eval(gl,d));}
    return best;
}

ld lower_bound(const Box& x,const PB& pb){
    std::array<ld,D> c,w;for(int k=0;k<D;k++){c[k]=(x.l[k]+x.u[k])/2;w[k]=(x.u[k]-x.l[k])/2;}
    std::array<I,D> gr; I ph=phi_center(c,gr); ld m=m_lower(pb); I s(ph.lower());
    for(int k=0;k<D;k++)s += I(quad_min_1d(gr[k],w[k],m));
    return s.lower();
}

bool root_possible(const Box& x,const PB& pb, int &reason, ld &gap){
    const ld one3=1.0L/3.0L;
    if(pb.m.upper() < one3){reason=1;gap=one3-pb.m.upper();return false;}
    if(x.u[0] < x.l[1]){reason=2;gap=x.l[1]-x.u[0];return false;}
    if(x.u[0] < -x.u[1]){reason=3;gap=-x.u[1]-x.u[0];return false;}
    // Do not use interior stationary equations on boxes that can touch triple junction.
    if(pb.m.lower() <= one3) return true;
    for(int k=0;k<D;k++){
      std::array<I,Q> co;for(int j=0;j<Q;j++)co[j]=B[j][k];
      ld muL=lin_ext(co,pb,false), muU=lin_ext(co,pb,true);
      if(!std::isfinite((double)muL)||!std::isfinite((double)muU))return true;
      I t=G*I(muL,muU); ld lo=t.lower(),hi=t.upper();
      if(x.u[k]<lo){reason=10+k;gap=lo-x.u[k];return false;}
      if(x.l[k]>hi){reason=20+k;gap=x.l[k]-hi;return false;}
    }
    return true;
}

bool inside_local(const Box& x){for(int r=0;r<4;r++){bool ok=true;for(int k=0;k<D;k++)if(x.l[k]<roots4[r][k]-LOCAL_R||x.u[k]>roots4[r][k]+LOCAL_R){ok=false;break;}if(ok)return true;}return false;}

int main(){
  for(int j=0;j<Q;j++)for(int k=0;k<D;k++)B[j][k]=I(B0[j][k]-BR,B0[j][k]+BR);
  initL[0]=-2.94179838416301L; initU[0]=2.94179838416301L;
  initL[1]=-2.94179838416301L; initU[1]=2.94179838416301L;
  initL[2]=-3.11528540297832L; initU[2]=1.55764270148917L;
  initL[3]=-3.18464739119806L; initU[3]=3.18464739119806L;
  initL[4]=-3.18464739119806L; initU[4]=3.18464739119806L;
  initL[5]=-2.27313058484033L; initU[5]=2.27313058484033L;
  ld rr[4][D]={
   {2.5703048563159236473L,-2.7061686225238190673e-16L,1.2992516656905954697L,0.68179265329213267766L,-1.1808995157291635181L,1.9759863787908105159L},
   {-5.898059818321144121e-16L,-2.5703048563159227591L,1.2992516656905956918L,1.1808995157291650724L,0.68179265329213201152L,-1.975986378790810738L},
   {-2.5703048563159232032L,3.4937330806172894881e-15L,1.2992516656905934713L,-0.6817926532921338989L,1.1808995157291677369L,1.9759863787908105159L},
   {9.1072982488782372457e-16L,2.5703048563159232032L,1.2992516656905939154L,-1.1808995157291624079L,-0.68179265329213234459L,-1.975986378790810738L}};
  for(int r=0;r<4;r++)for(int k=0;k<D;k++)roots4[r][k]=rr[r][k];
  Box init;for(int k=0;k<D;k++){init.l[k]=initL[k];init.u[k]=initU[k];}init.depth=0;
  std::vector<Box> st;st.reserve(2000000);st.push_back(init);
  uint64_t nodes=0,df=0,de=0,local=0; int md=0; ld minFG=1e100L,minEM=1e100L; int minFR=0;
  while(!st.empty() && nodes<100000000ULL){
    Box x=st.back();st.pop_back();nodes++;md=std::max(md,(int)x.depth);
    if(inside_local(x)){local++;continue;}
    PB pb=p_bounds(x);int reason=0;ld gap=0;
    if(!root_possible(x,pb,reason,gap)){df++;if(gap<minFG){minFG=gap;minFR=reason;}continue;}
    ld lb=lower_bound(x,pb);ld margin=lb-TARGET;
    if(margin>0){de++;minEM=std::min(minEM,margin);continue;}
    int k=0;ld best=-1;for(int j=0;j<D;j++){ld sc=(x.u[j]-x.l[j])/(initU[j]-initL[j]);if(sc>best){best=sc;k=j;}}
    ld width=x.u[k]-x.l[k];if(x.depth>130||width<2e-8L){std::cerr<<"UNRESOLVED depth="<<x.depth<<" lb="<<(double)lb<<" width="<<(double)width<<"\n";return 3;}
    ld mid=(x.l[k]+x.u[k])/2;Box a=x,b=x;a.u[k]=mid;b.l[k]=mid;a.depth=b.depth=x.depth+1;st.push_back(b);st.push_back(a);
    if(nodes%250000==0)std::cerr<<"nodes "<<nodes<<" stack "<<st.size()<<" df "<<df<<" de "<<de<<" local "<<local<<" depth "<<md<<"\n";
  }
  std::cout<<std::setprecision(20)
    <<"nodes "<<nodes<<"\ndiscardFeas "<<df<<"\ndiscardEnergy "<<de<<"\nlocal "<<local<<"\nremaining "<<st.size()<<"\nmaxDepth "<<md
    <<"\nminFeasGap "<<minFG<<"\nminFeasReason "<<minFR<<"\nminEnergyMargin "<<minEM<<"\n";
  return st.empty()?0:2;
}
