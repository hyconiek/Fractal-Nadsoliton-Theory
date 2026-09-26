#include <bits/stdc++.h>
#include <cfenv>
#include <cmath>
#pragma STDC FENV_ACCESS ON
using namespace std;
struct Term{array<int,4> e; double c;};
struct Box{array<double,3> l,u; double ub;};
static const Term T[]={
 {{{2,2,0,0}},0.04822142787987652},
 {{{2,0,2,0}},0.05333608476307321},
 {{{1,2,1,0}},0.11385707390148081},
 {{{1,1,1,1}},0.3428679997023158},
 {{{1,0,3,0}},0.06204319383393708},
 {{{0,2,2,0}},0.0656356460216043},
 {{{0,2,0,2}},0.13577211811522824},
 {{{0,1,2,1}},0.2940921349879693},
 {{{0,0,4,0}},0.017687575726200243},
 {{{0,0,2,2}},0.1460014318816216}
};
static inline double mul_up(double a,double b){volatile double x=a*b;return x;}
static inline double add_up(double a,double b){volatile double x=a+b;return x;}
static inline double sqrt_up(double a){volatile double x=sqrt(a);return x;}
static inline double ph(double y,int e){if(!e)return 1.;if(y<=0)return 0.;if(e==1)return sqrt_up(y);if(e==2)return y;if(e==3)return mul_up(y,sqrt_up(y));if(e==4)return mul_up(y,y);abort();}
static inline bool derive(const Box&b,array<double,4>&l,array<double,4>&u){for(int i=0;i<3;i++){l[i]=b.l[i];u[i]=b.u[i];}double sl=b.l[0]+b.l[1]+b.l[2],su=b.u[0]+b.u[1]+b.u[2];if(sl>1)return false;l[3]=max(0.,1-su);u[3]=min(1.,1-sl);return l[3]<=u[3];}
static double mub(const Term&t,const array<double,4>&l,const array<double,4>&u){double box=1.;double inactive=0;for(int i=0;i<4;i++){if(t.e[i])box=mul_up(box,ph(u[i],t.e[i]));else inactive+=l[i];}if(box==0)return 0;double S=1-inactive,am=1.;if(S<0)return 0;for(int i=0;i<4;i++)if(t.e[i])am=mul_up(am,ph(S*(.25*t.e[i]),t.e[i]));return min(box,am);}
static double UB(const Box&b){array<double,4>l,u;if(!derive(b,l,u))return -INFINITY;double s=0;for(auto&t:T)s=add_up(s,mul_up(t.c,mub(t,l,u)));return s;}
int main(int argc,char**argv){fesetround(FE_UPWARD);double best=0.0932487158813161;double width=argc>1?strtod(argv[1],0):0.002;vector<Box> stack;Box r{{0,0,0},{1,1,1},0};r.ub=UB(r);stack.push_back(r);array<double,4>GL={1,1,1,1},GU={0,0,0,0};long long nodes=0,leaves=0;vector<array<double,4>> ctrs;while(!stack.empty()){Box b=stack.back();stack.pop_back();nodes++;if(b.ub<best)continue;array<double,4>l,u;if(!derive(b,l,u))continue;double w=0;int k=0;for(int i=0;i<3;i++)if(b.u[i]-b.l[i]>w){w=b.u[i]-b.l[i];k=i;}if(w<=width){leaves++;for(int i=0;i<4;i++){GL[i]=min(GL[i],l[i]);GU[i]=max(GU[i],u[i]);}array<double,4>c;for(int i=0;i<4;i++)c[i]=(l[i]+u[i])*.5;ctrs.push_back(c);continue;}double m=(b.l[k]+b.u[k])*.5;for(int s=0;s<2;s++){Box q=b;if(!s)q.u[k]=m;else q.l[k]=m;q.ub=UB(q);if(q.ub>=best)stack.push_back(q);}}
printf("nodes %lld leaves %lld width %.17g\n",nodes,leaves,width);for(int i=0;i<4;i++)printf("y%d [%.17g, %.17g]\n",i+3,GL[i],GU[i]);
// crude max center spread relative to known optimum
array<double,4> y0={0.08598024,0.23653999,0.42523964,0.25224014};double md=0;for(auto&c:ctrs){double d=0;for(int i=0;i<4;i++)d=max(d,abs(c[i]-y0[i]));md=max(md,d);}printf("max_center_Linf_from_reference %.17g\n",md);
}
