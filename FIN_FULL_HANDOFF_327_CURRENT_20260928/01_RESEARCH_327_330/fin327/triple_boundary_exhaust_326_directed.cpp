#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/transc.hpp>
#include <array>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cstdint>
namespace bn=boost::numeric; namespace il=boost::numeric::interval_lib; using ld=long double;
using I=bn::interval<ld, il::policies<il::save_state<il::rounded_transc_std<ld>>, il::checking_base<ld>>>;
constexpr int Q=12,D=5; const I G=I(5.145228719489141L,5.145228719489143L); const ld TARGET=-0.1431003250148460L,BR=5e-18L;
ld B0[Q][D]={
  {0.57175269449539189992L,0L,0.6189515694676895885L,0L,0.44179388493073518118L},
  {3.500975536088280713e-17L,0.57175269449539189992L,-0.53602778287126795487L,0.30947578473384473874L,-0.44179388493073518118L},
  {-0.57175269449539189992L,7.0019510721765614259e-17L,0.30947578473384484976L,-0.53602778287126784384L,0.44179388493073518118L},
  {-1.0502926608264842139e-16L,-0.57175269449539189992L,7.392386914346839704e-16L,0.6189515694676895885L,-0.44179388493073518118L},
  {0.57175269449539189992L,-1.4003902144353122852e-16L,-0.30947578473384468323L,-0.53602778287126795487L,0.44179388493073518118L},
  {6.8286718146061965487e-16L,0.57175269449539189992L,0.5360277828712676218L,0.30947578473384512732L,-0.44179388493073518118L},
  {-0.57175269449539189992L,2.1005853216529684278e-16L,-0.6189515694676895885L,1.4784773828693679408e-15L,0.44179388493073518118L},
  {-2.4506828752617964991e-16L,-0.57175269449539189992L,0.53602778287126784384L,-0.30947578473384484976L,-0.44179388493073518118L},
  {0.57175269449539189992L,-2.8007804288706245704e-16L,-0.30947578473384501629L,0.53602778287126784384L,0.44179388493073518118L},
  {3.1508779824794531347e-16L,0.57175269449539189992L,5.3098105989957209857e-16L,-0.6189515694676895885L,-0.44179388493073518118L},
  {-0.57175269449539189992L,1.3657343629212393097e-15L,0.3094757847338440726L,0.53602778287126839896L,0.44179388493073518118L},
  {6.3052950034270040933e-16L,-0.57175269449539189992L,-0.53602778287126728873L,-0.30947578473384573794L,-0.44179388493073518118L}
}; I B[Q][D]; int cls(int j){return j%3;} struct Box{std::array<ld,D>l,u;uint16_t dep;};
struct PB{I p[Q];};
PB pbounds(const Box&x){I e[Q],S[3]={I(0),I(0),I(0)};for(int j=0;j<Q;j++){I h(0);for(int k=0;k<D;k++)h+=B[j][k]*I(x.l[k],x.u[k]);e[j]=exp(h);S[cls(j)]+=e[j];}PB pb;for(int a=0;a<3;a++)for(int j=0;j<Q;j++)if(cls(j)==a)pb.p[j]=I(1.0L/3.0L)*e[j]/S[a];return pb;}
ld extclass(const std::array<I,Q>&co,const PB&pb,int a,bool mx){struct T{ld c,L,U;};T t[4];int n=0;ld rem=1.0L/3.0L;I sum(0);for(int j=0;j<Q;j++)if(cls(j)==a){t[n++]={mx?co[j].upper():co[j].lower(),pb.p[j].lower(),pb.p[j].upper()};rem-=pb.p[j].lower();}std::sort(t,t+n,[&](auto&A,auto&B){return mx?A.c>B.c:A.c<B.c;});for(int i=0;i<n;i++)sum+=I(t[i].c)*I(t[i].L);for(int i=0;i<n&&rem>0;i++){ld add=std::min(rem,std::max((ld)0,t[i].U-t[i].L));sum+=I(t[i].c)*I(add);rem-=add;}if(rem>1e-18L)return mx?1e100L:-1e100L;return mx?sum.upper():sum.lower();}
ld ext(const std::array<I,Q>&co,const PB&pb,bool mx){I s(0);for(int a=0;a<3;a++)s+=I(extclass(co,pb,a,mx));return mx?s.upper():s.lower();}
I phic(const std::array<ld,D>&y,std::array<I,D>&gr){I e[Q],S[3]={I(0),I(0),I(0)};for(int j=0;j<Q;j++){I h(0);for(int k=0;k<D;k++)h+=B[j][k]*I(y[k]);e[j]=exp(h);S[cls(j)]+=e[j];}I p[Q];for(int j=0;j<Q;j++)p[j]=I(1.0L/3.0L)*e[j]/S[cls(j)];I norm(0);for(int k=0;k<D;k++)norm+=I(y[k])*I(y[k]);I phi=norm/(I(2)*G)-(log(S[0])+log(S[1])+log(S[2]))/I(3)+log(I(4));for(int k=0;k<D;k++){I mu(0);for(int j=0;j<Q;j++)mu+=p[j]*B[j][k];gr[k]=I(y[k])/G-mu;}return phi;}
ld mlower(const PB&pb){I tr(0);for(int a=0;a<3;a++){std::array<I,Q> normco;for(int j=0;j<Q;j++){I s(0);for(int k=0;k<D;k++)s+=B[j][k]*B[j][k];normco[j]=s;}tr+=I(extclass(normco,pb,a,true));}return (I(1)/G-tr).lower();}
ld qmin(const I&g,ld w,ld m){ld gl=g.lower(),gu=g.upper(),best=0;auto ev=[&](ld gg,ld d){return (I(gg)*I(d)+I(.5L)*I(m)*I(d)*I(d)).lower();};best=std::min(best,ev(gu,-w));best=std::min(best,ev(gl,w));if(m>0){ld d=-gu/m;if(d>=-w&&d<=0)best=std::min(best,ev(gu,d));d=-gl/m;if(d>=0&&d<=w)best=std::min(best,ev(gl,d));}return best;}
ld lower(const Box&x,const PB&pb){std::array<ld,D>c,w;for(int k=0;k<D;k++){c[k]=(x.l[k]+x.u[k])/2;w[k]=(x.u[k]-x.l[k])/2;}std::array<I,D>gr;I ph=phic(c,gr);ld m=mlower(pb);I s(ph.lower());for(int k=0;k<D;k++)s+=I(qmin(gr[k],w[k],m));return s.lower();}
bool possible(const Box&x,const PB&pb,int&reason,ld&gap){if(x.u[0]<x.l[1]){reason=1;gap=x.l[1]-x.u[0];return false;}if(x.u[0]<-x.u[1]){reason=2;gap=-x.u[1]-x.u[0];return false;}for(int k=0;k<D;k++){std::array<I,Q>co;for(int j=0;j<Q;j++)co[j]=B[j][k];ld lo=ext(co,pb,false),hi=ext(co,pb,true);I t=G*I(lo,hi);if(x.u[k]<t.lower()){reason=10+k;gap=t.lower()-x.u[k];return false;}if(x.l[k]>t.upper()){reason=20+k;gap=x.l[k]-t.upper();return false;}}return true;}
int main(){for(int j=0;j<Q;j++)for(int k=0;k<D;k++)B[j][k]=I(B0[j][k]-BR,B0[j][k]+BR);ld L[D]={-2.94179838416301L,-2.94179838416301L,-3.18464739119806L,-3.18464739119806L,-2.27313058484033L};ld U[D]={2.94179838416301L,2.94179838416301L,3.18464739119806L,3.18464739119806L,2.27313058484033L};Box z;for(int k=0;k<D;k++){z.l[k]=L[k];z.u[k]=U[k];}z.dep=0;std::vector<Box>st;st.push_back(z);uint64_t n=0,df=0,de=0;ld minfg=1e100L,minem=1e100L;int minr=0,md=0;while(!st.empty()&&n<100000000){Box x=st.back();st.pop_back();n++;md=std::max(md,(int)x.dep);PB pb=pbounds(x);int r=0;ld gap=0;if(!possible(x,pb,r,gap)){df++;if(gap<minfg){minfg=gap;minr=r;}continue;}ld lb=lower(x,pb),mar=lb-TARGET;if(mar>0){de++;minem=std::min(minem,mar);continue;}int k=0;ld best=-1;for(int j=0;j<D;j++){ld sc=(x.u[j]-x.l[j])/(U[j]-L[j]);if(sc>best){best=sc;k=j;}}ld width=x.u[k]-x.l[k];if(x.dep>130||width<2e-8L){std::cerr<<"UNRESOLVED "<<(double)lb<<"\n";return 3;}ld m=(x.l[k]+x.u[k])/2;Box a=x,b=x;a.u[k]=m;b.l[k]=m;a.dep=b.dep=x.dep+1;st.push_back(b);st.push_back(a);}std::cout<<std::setprecision(20)<<"nodes "<<n<<"\ndiscardFeas "<<df<<"\ndiscardEnergy "<<de<<"\nremaining "<<st.size()<<"\nmaxDepth "<<md<<"\nminFeasGap "<<minfg<<"\nminFeasReason "<<minr<<"\nminEnergyMargin "<<minem<<"\n";return st.empty()?0:2;}
