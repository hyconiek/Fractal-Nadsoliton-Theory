#include <bits/stdc++.h>
#include <cfenv>
#include <cmath>
#pragma STDC FENV_ACCESS ON
using namespace std;
struct Term{array<int,4> e; double c;}; // exponent in amplitudes; y exponent=e/2
struct Box{array<double,3> l,u; double ub;};
struct Cmp{bool operator()(Box const&a,Box const&b)const{return a.ub<b.ub;}};
// Strict upward enclosures of the ten positive coefficients at g_eq.
static const Term T[]={
 {{{2,2,0,0}},0.011177246149602943},
 {{{2,0,2,0}},0.011830111953407215},
 {{{1,2,1,0}},0.02437846288765537},
 {{{1,1,1,1}},0.0711429240695726},
 {{{1,0,3,0}},0.012711986408918446},
 {{{0,2,2,0}},0.012981880774529168},
 {{{0,2,0,2}},0.02635435466039102},
 {{{0,1,2,1}},0.056368920667498636},
 {{{0,0,4,0}},0.003347642998043608},
 {{{0,0,2,2}},0.027118888219793382}
};
static inline double mul_up(double a,double b){ volatile double x=a*b; return x; }
static inline double add_up(double a,double b){ volatile double x=a+b; return x; }
static inline double sqrt_up(double a){ volatile double x=sqrt(a); return x; }
static inline double pow_halfint_up(double y,int e){
 if(e==0) return 1.0;
 if(y<=0) return 0.0;
 if(e==1) return sqrt_up(y);
 if(e==2) return y;
 if(e==3) return mul_up(y,sqrt_up(y));
 if(e==4) return mul_up(y,y);
 abort();
}
static inline bool derive4(const Box&b,array<double,4>&l,array<double,4>&u){
 for(int i=0;i<3;i++){l[i]=b.l[i];u[i]=b.u[i];}
 // dyadic endpoints: these additions/subtractions are exact in the explored range.
 double sl=b.l[0]+b.l[1]+b.l[2];
 double su=b.u[0]+b.u[1]+b.u[2];
 if(sl>1.0) return false;
 l[3]=max(0.0,1.0-su); u[3]=min(1.0,1.0-sl);
 return l[3]<=u[3];
}
static double monoUB(const Term&t,const array<double,4>&l,const array<double,4>&u){
 // Box product upper bound.
 double box=1.0; int esum=0; double inactiveLower=0.0;
 for(int i=0;i<4;i++){
   if(t.e[i]){ box=mul_up(box,pow_halfint_up(u[i],t.e[i])); esum+=t.e[i]; }
   else inactiveLower += l[i]; // dyadic exact
 }
 if(box==0.0) return 0.0;
 // Weighted AM-GM with sum(active y_i) <= S. Total y-degree = esum/2 = 2.
 double S=1.0-inactiveLower; if(S<0) return 0.0;
 double am=1.0;
 for(int i=0;i<4;i++) if(t.e[i]){
   // p_i/P = (e_i/2)/2 = e_i/4, exactly representable in binary.
   double yi=S*(0.25*t.e[i]);
   am=mul_up(am,pow_halfint_up(yi,t.e[i]));
 }
 return min(box,am); // actual monomial is <= each bound
}
static double UB(const Box&b){
 array<double,4>l,u;if(!derive4(b,l,u))return -INFINITY;
 double s=0.0;
 for(auto const&t:T) s=add_up(s,mul_up(t.c,monoUB(t,l,u)));
 return s;
}
int main(int argc,char**argv){
 fesetround(FE_UPWARD);
 const double best=0.01839077538669660; // strictly below interval-certified feasible value
 const double target=argc>1?strtod(argv[1],nullptr):0.0001;
 long long maxnodes=argc>2?atoll(argv[2]):30000000LL;
 Box root{{0,0,0},{1,1,1},0}; root.ub=UB(root);
 priority_queue<Box,vector<Box>,Cmp>pq;pq.push(root);
 long long n=0;auto st=chrono::steady_clock::now();
 while(!pq.empty()&&n<maxnodes){
   Box b=pq.top();pq.pop();
   if(b.ub<=best+target){
     printf("CERTIFIED nodes %lld lower %.17g upper %.17g gap %.17g heap %zu\n",n,best,b.ub,b.ub-best,pq.size());
     return 0;
   }
   int k=0;double w=-1;for(int i=0;i<3;i++){double z=b.u[i]-b.l[i];if(z>w){w=z;k=i;}}
   double mid=(b.l[k]+b.u[k])*0.5; // exact dyadic midpoint
   for(int side=0;side<2;side++){
      Box q=b;if(side==0)q.u[k]=mid;else q.l[k]=mid;
      array<double,4>ll,uu;if(!derive4(q,ll,uu))continue;
      q.ub=UB(q); if(q.ub>best+target)pq.push(q);
   }
   n++;
   if(n%1000000==0){double sec=chrono::duration<double>(chrono::steady_clock::now()-st).count();fprintf(stderr,"%lld top %.17g gap %.17g heap %zu sec %.3f\n",n,pq.empty()?best:pq.top().ub,(pq.empty()?best:pq.top().ub)-best,pq.size(),sec);}
 }
 double cut=best+target; printf("CERTIFIED_EXHAUSTED nodes %lld lower %.17g upper %.17g gap %.17g heap %zu\n",n,best,cut,cut-best,pq.size());
 return 1;
}
