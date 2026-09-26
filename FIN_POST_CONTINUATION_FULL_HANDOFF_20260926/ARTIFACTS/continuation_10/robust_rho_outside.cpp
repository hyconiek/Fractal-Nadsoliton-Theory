#include <bits/stdc++.h>
using namespace std;
struct Term{array<int,4> e; double c,h,q;};
struct Box{array<double,3> l,u; double ub;}; struct Cmp{bool operator()(Box const&a,Box const&b)const{return a.ub<b.ub;}};
static const Term T[]={
 {{{2,2,0,0}},0.04822142787987652,0.013888888888888891,0.1674220396877376},
 {{{2,0,2,0}},0.05333608476307321,0.013888888888888891,0.2048211315254686},
 {{{1,2,1,0}},0.11385707390148081,0.027777777777777783,0.46668359798666154},
 {{{1,1,1,1}},0.3428679997023158,0.07856742013183864,1.4962749829713182},
 {{{1,0,3,0}},0.06204319383393708,0.013888888888888891,0.27715376888031523},
 {{{0,2,2,0}},0.0656356460216043,0.013888888888888891,0.31017873806448054},
 {{{0,2,0,2}},0.13577211811522824,0.027777777777777783,0.6636264500698378},
 {{{0,1,2,1}},0.2940921349879693,0.05892556509887898,1.4677870923536958},
 {{{0,0,4,0}},0.017687575726200243,0.003472222222222223,0.0901008965001796},
 {{{0,0,2,2}},0.1460014318816216,0.027777777777777783,0.7673910520134166}};
static inline bool derive(const Box&b,array<double,4>&l,array<double,4>&u){for(int i=0;i<3;i++){l[i]=b.l[i];u[i]=b.u[i];}double sl=b.l[0]+b.l[1]+b.l[2],su=b.u[0]+b.u[1]+b.u[2];if(sl>1)return false;l[3]=max(0.0,1-su);u[3]=min(1.0,1-sl);return l[3]<=u[3];}
static inline double pw(double y,int e){if(e==0)return 1;if(y<=0)return 0;if(e==1)return sqrt(y);if(e==2)return y;if(e==3)return y*sqrt(y);return y*y;}
static double monoUB(const Term&t,const array<double,4>&l,const array<double,4>&u){double box=1;double inactive=0;for(int i=0;i<4;i++){if(t.e[i])box*=pw(u[i],t.e[i]);else inactive+=l[i];}double S=max(0.0,1-inactive),am=1;for(int i=0;i<4;i++)if(t.e[i])am*=pw(S*(.25*t.e[i]),t.e[i]);return min(box,am);} 
static double monoLB(const Term&t,const array<double,4>&l){double z=1;for(int i=0;i<4;i++)if(t.e[i])z*=pw(l[i],t.e[i]);return z;}
static double UB(const Box&b){array<double,4>l,u;if(!derive(b,l,u))return -INFINITY;long double CU=0,HU=0,HL=0,QU=0;for(auto&t:T){double mu=monoUB(t,l,u);CU+=(long double)t.c*mu;HU+=(long double)t.h*mu;QU+=(long double)t.q*mu;HL+=(long double)t.h*monoLB(t,l);}long double b0=(4.0L/3.7183448981203875L)*sqrt(QU);long double b1=5.654184L*sqrt((long double)HU);long double b2=INFINITY;if(HL>0)b2=4.0L*CU/(3.7183448981203875L*sqrt((long double)HL));return (double)min(b0,min(b1,b2));}
static array<double,3> RL={0.058,0.215,0.427};
static array<double,3> RU={0.070,0.237,0.450};
static bool insideR(const Box&b){for(int i=0;i<3;i++) if(b.l[i]<RL[i] || b.u[i]>RU[i]) return false; return true;}
static bool disjointR(const Box&b){for(int i=0;i<3;i++) if(b.u[i]<=RL[i] || b.l[i]>=RU[i]) return true; return false;}
int main(int argc,char**argv){double best=0.7175332680535;long long maxn=argc>1?atoll(argv[1]):30000000LL;Box root{{0,0,0},{1,1,1},0};root.ub=UB(root);priority_queue<Box,vector<Box>,Cmp>pq;pq.push(root);long long n=0,skipped=0;double maxskip=0;auto st=chrono::steady_clock::now();while(!pq.empty() && n<maxn){Box b=pq.top();pq.pop(); if(b.ub<=best){cout<<setprecision(17)<<"CERTIFIED_OUTSIDE nodes "<<n<<" outside_upper "<<b.ub<<" skipped "<<skipped<<" maxskip "<<maxskip<<"\n"; return 0;} if(insideR(b)){skipped++;maxskip=max(maxskip,b.ub);continue;} int k=0; double w=-1; // prefer split along a boundary-crossing coordinate
 for(int i=0;i<3;i++){double z=b.u[i]-b.l[i]; bool cross=!(b.u[i]<=RL[i]||b.l[i]>=RU[i]) && !(b.l[i]>=RL[i]&&b.u[i]<=RU[i]); double score=z*(cross?10.0:1.0); if(score>w){w=score;k=i;}}
 double mid=.5*(b.l[k]+b.u[k]); // if crossing region boundary, split exactly there when interior
 if(b.l[k]<RL[k] && RL[k]<b.u[k]) mid=RL[k]; else if(b.l[k]<RU[k] && RU[k]<b.u[k]) mid=RU[k];
 for(int s=0;s<2;s++){Box q=b;if(!s)q.u[k]=mid;else q.l[k]=mid;array<double,4>l,u;if(!derive(q,l,u))continue;q.ub=UB(q);if(q.ub>best)pq.push(q);} n++; if(n%1000000==0)cerr<<n<<" top "<<(pq.empty()?best:pq.top().ub)<<" heap "<<pq.size()<<" skip "<<skipped<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"\n";
 }
 cout<<setprecision(17)<<"INCOMPLETE nodes "<<n<<" top "<<(pq.empty()?best:pq.top().ub)<<" heap "<<pq.size()<<" skipped "<<skipped<<" maxskip "<<maxskip<<"\n";return 1;}
