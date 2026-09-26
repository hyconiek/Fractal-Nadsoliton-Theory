#include <bits/stdc++.h>
#include <cfenv>
#include <cmath>
#pragma STDC FENV_ACCESS ON
using namespace std;
struct Term{array<int,4> e; double cU,hU,hL,qU;};
struct Box{array<double,3> l,u;double ub;}; struct Cmp{bool operator()(Box const&a,Box const&b)const{return a.ub<b.ub;}};
// C coefficients are strict upward enclosures copied from the certified quartic BnB.
// H coefficients bracket exact positive coefficients of ||P_H(phi^2)||_u^2.
static const Term T[]={
 {{{2,2,0,0}},0.04822142787987652,0.013888888888888891,0.013888888888888886,0.1674220396877376},
 {{{2,0,2,0}},0.05333608476307321,0.013888888888888891,0.013888888888888886,0.2048211315254686},
 {{{1,2,1,0}},0.11385707390148081,0.027777777777777783,0.027777777777777773,0.46668359798666154},
 {{{1,1,1,1}},0.3428679997023158,0.07856742013183864,0.07856742013183860,1.4962749829713182},
 {{{1,0,3,0}},0.06204319383393708,0.013888888888888891,0.013888888888888886,0.27715376888031523},
 {{{0,2,2,0}},0.0656356460216043,0.013888888888888891,0.013888888888888886,0.31017873806448054},
 {{{0,2,0,2}},0.13577211811522824,0.027777777777777783,0.027777777777777773,0.6636264500698378},
 {{{0,1,2,1}},0.2940921349879693,0.05892556509887898,0.05892556509887894,1.4677870923536958},
 {{{0,0,4,0}},0.017687575726200243,0.003472222222222223,0.0034722222222222215,0.0901008965001796},
 {{{0,0,2,2}},0.1460014318816216,0.027777777777777783,0.027777777777777773,0.7673910520134166}
};
static inline bool derive4(const Box&b,array<double,4>&l,array<double,4>&u){for(int i=0;i<3;i++){l[i]=b.l[i];u[i]=b.u[i];}double sl=b.l[0]+b.l[1]+b.l[2],su=b.u[0]+b.u[1]+b.u[2];if(sl>1)return false;l[3]=max(0.0,1-su);u[3]=min(1.0,1-sl);return l[3]<=u[3];}
static inline double pwh(double y,int e){if(e==0)return 1;if(y<=0)return 0;if(e==1)return sqrt(y);if(e==2)return y;if(e==3)return y*sqrt(y);if(e==4)return y*y;abort();}
static double monoUB(const Term&t,const array<double,4>&l,const array<double,4>&u){double box=1;int es=0;double inact=0;for(int i=0;i<4;i++){if(t.e[i]){box*=pwh(u[i],t.e[i]);es+=t.e[i];}else inact+=l[i];}double S=max(0.0,1-inact),am=1;for(int i=0;i<4;i++)if(t.e[i])am*=pwh(S*(0.25*t.e[i]),t.e[i]);return min(box,am);}
static double monoLB(const Term&t,const array<double,4>&l){double z=1;for(int i=0;i<4;i++)if(t.e[i])z*=pwh(l[i],t.e[i]);return z;}
static double UB(const Box&b){array<double,4>l,u;if(!derive4(b,l,u))return -INFINITY;long double CU=0,HU=0,HL=0,QU=0;for(auto&t:T){double mu=monoUB(t,l,u);CU+=(long double)t.cU*mu;HU+=(long double)t.hU*mu;QU+=(long double)t.qU*mu;HL+=(long double)t.hL*monoLB(t,l);}long double b0=(4.0L/3.7183L)*sqrt(QU);long double b1=5.829L*sqrt((long double)HU);long double b2=INFINITY;if(HL>0)b2=4.0L*CU/(3.7183L*sqrt((long double)HL));return (double)min(b0,min(b1,b2));}
int main(int argc,char**argv){const double best=argc>1?atof(argv[1]):0.7175330;const double target=argc>2?atof(argv[2]):0.0001;long long maxn=argc>3?atoll(argv[3]):30000000LL;Box r{{0,0,0},{1,1,1},0};r.ub=UB(r);priority_queue<Box,vector<Box>,Cmp>pq;pq.push(r);long long n=0;auto st=chrono::steady_clock::now();while(!pq.empty()&&n<maxn){Box b=pq.top();pq.pop();if(b.ub<=best+target){cout<<setprecision(17)<<"CERTIFIED nodes "<<n<<" lower "<<best<<" upper "<<b.ub<<" gap "<<b.ub-best<<" heap "<<pq.size()<<"\n";return 0;}int k=0;double w=-1;for(int i=0;i<3;i++){double q=b.u[i]-b.l[i];if(q>w){w=q;k=i;}}double m=(b.l[k]+b.u[k])*0.5;for(int side=0;side<2;side++){Box q=b;if(!side)q.u[k]=m;else q.l[k]=m;array<double,4>ll,uu;if(!derive4(q,ll,uu))continue;q.ub=UB(q);if(q.ub>best+target)pq.push(q);}n++;if(n%1000000==0){double sec=chrono::duration<double>(chrono::steady_clock::now()-st).count();cerr<<n<<" top "<<(pq.empty()?best:pq.top().ub)<<" heap "<<pq.size()<<" sec "<<sec<<"\n";}}
cout<<setprecision(17)<<"INCOMPLETE nodes "<<n<<" top "<<(pq.empty()?best:pq.top().ub)<<" heap "<<pq.size()<<"\n";return 1;}
