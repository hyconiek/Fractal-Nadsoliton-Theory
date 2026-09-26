#include <bits/stdc++.h>
#include <cfenv>
#include <cmath>
#pragma STDC FENV_ACCESS ON
using namespace std;
struct Term{array<int,4> e; double c;};
struct Box{array<double,3> l,u; double ub;};
struct Cmp{bool operator()(Box const&a,Box const&b)const{return a.ub<b.ub;}};
static const Term T[]={
 {{{2,2,0,0}},0.04822142787987652},{{{2,0,2,0}},0.05333608476307321},{{{1,2,1,0}},0.11385707390148081},
 {{{1,1,1,1}},0.3428679997023158},{{{1,0,3,0}},0.06204319383393708},{{{0,2,2,0}},0.0656356460216043},
 {{{0,2,0,2}},0.13577211811522824},{{{0,1,2,1}},0.2940921349879693},{{{0,0,4,0}},0.017687575726200243},
 {{{0,0,2,2}},0.1460014318816216}};
static inline double mul_up(double a,double b){volatile double x=a*b;return x;}
static inline double add_up(double a,double b){volatile double x=a+b;return x;}
static inline double sqrt_up(double a){volatile double x=sqrt(a);return x;}
static inline double ph(double y,int e){if(e==0)return 1.;if(y<=0)return 0.;if(e==1)return sqrt_up(y);if(e==2)return y;if(e==3)return mul_up(y,sqrt_up(y));if(e==4)return mul_up(y,y);abort();}
static inline bool derive(const Box&b,array<double,4>&l,array<double,4>&u){for(int i=0;i<3;i++){l[i]=b.l[i];u[i]=b.u[i];}double sl=b.l[0]+b.l[1]+b.l[2],su=b.u[0]+b.u[1]+b.u[2];if(sl>1.)return false;l[3]=max(0.,1.-su);u[3]=min(1.,1.-sl);return l[3]<=u[3];}
static double mono(const Term&t,const array<double,4>&l,const array<double,4>&u){double box=1.;double inactive=0.;for(int i=0;i<4;i++){if(t.e[i])box=mul_up(box,ph(u[i],t.e[i]));else inactive+=l[i];}if(box==0)return 0.;double S=1.-inactive;if(S<0)return 0.;double am=1.;for(int i=0;i<4;i++)if(t.e[i])am=mul_up(am,ph(S*(0.25*t.e[i]),t.e[i]));return min(box,am);}
static double UB(const Box&b){array<double,4>l,u;if(!derive(b,l,u))return -INFINITY;double s=0.;for(auto&t:T)s=add_up(s,mul_up(t.c,mono(t,l,u)));return s;}
int main(int argc,char**argv){fesetround(FE_UPWARD);double best=0.0932487158813161;double tol=argc>1?atof(argv[1]):1e-5;long long maxnodes=argc>2?atoll(argv[2]):50000000LL;
 priority_queue<Box,vector<Box>,Cmp> pq;Box root{{0,0,0},{1,1,1},0};root.ub=UB(root);pq.push(root);
 array<double,4> glo{1,1,1,1},ghi{0,0,0,0}; long long frontier=0,n=0; double maxfront=best;
 auto keep=[&](const Box&q){if(q.ub<best)return; if(q.ub<=best+tol){array<double,4>l,u;if(!derive(q,l,u))return;for(int i=0;i<4;i++){glo[i]=min(glo[i],l[i]);ghi[i]=max(ghi[i],u[i]);}frontier++;maxfront=max(maxfront,q.ub);}else pq.push(q);};
 auto st=chrono::steady_clock::now();
 while(!pq.empty()&&n<maxnodes){Box b=pq.top();pq.pop();int k=0;double w=-1;for(int i=0;i<3;i++){double z=b.u[i]-b.l[i];if(z>w){w=z;k=i;}}double mid=(b.l[k]+b.u[k])*.5;for(int side=0;side<2;side++){Box q=b;if(!side)q.u[k]=mid;else q.l[k]=mid;array<double,4>l,u;if(!derive(q,l,u))continue;q.ub=UB(q);keep(q);}n++;if(n%1000000==0)fprintf(stderr,"%lld top %.17g heap %zu frontier %lld sec %.2f\n",n,pq.empty()?best:pq.top().ub,pq.size(),frontier,chrono::duration<double>(chrono::steady_clock::now()-st).count());}
 printf("nodes %lld heap %zu frontier %lld best %.17g tol %.17g maxfront %.17g\n",n,pq.size(),frontier,best,tol,maxfront);
 for(int i=0;i<4;i++)printf("y%d [%.17g, %.17g]  c%d [%.17g, %.17g]\n",i,glo[i],ghi[i],i,sqrt(glo[i]),sqrt(ghi[i]));
 return pq.empty()?0:2;}
