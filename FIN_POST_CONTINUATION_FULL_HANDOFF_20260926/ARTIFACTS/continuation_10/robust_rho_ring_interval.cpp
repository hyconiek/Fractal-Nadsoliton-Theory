#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/utility.hpp>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <queue>
using namespace std; using boost::numeric::interval; using boost::numeric::lower; using boost::numeric::upper; using I=interval<double>;
struct Term{array<int,4>e;double cL,cU,hL,hU,qU;};
static const Term T[]={
{{2,2,0,0},0.04822142787987651,0.04822142787987652,0.013888888888888886,0.013888888888888891,0.1674220396877376},
{{2,0,2,0},0.05333608476307319,0.0533360847630732,0.013888888888888886,0.013888888888888891,0.2048211315254686},
{{1,2,1,0},0.11385707390148077,0.1138570739014808,0.027777777777777773,0.027777777777777783,0.46668359798666154},
{{1,1,1,1},0.34286799970231563,0.34286799970231574,0.07856742013183860,0.07856742013183864,1.4962749829713182},
{{1,0,3,0},0.06204319383393707,0.06204319383393708,0.013888888888888886,0.013888888888888891,0.27715376888031523},
{{0,2,2,0},0.06563564602160425,0.06563564602160428,0.013888888888888886,0.013888888888888891,0.31017873806448054},
{{0,2,0,2},0.13577211811522816,0.1357721181152282,0.027777777777777773,0.027777777777777783,0.6636264500698378},
{{0,1,2,1},0.2940921349879692,0.2940921349879693,0.05892556509887894,0.05892556509887898,1.4677870923536958},
{{0,0,4,0},0.017687575726200232,0.01768757572620024,0.0034722222222222215,0.003472222222222223,0.0901008965001796},
{{0,0,2,2},0.14600143188162154,0.1460014318816216,0.027777777777777773,0.027777777777777783,0.7673910520134166}};
struct Box{array<double,3>l,u;double ub;};struct Cmp{bool operator()(Box const&a,Box const&b)const{return a.ub<b.ub;}};
static array<double,3> RL={0.063915514399,0.225833540948,0.438532296078},RU={0.063935514399,0.225853540948,0.438552296078};
I powi(I x,int n){I r(1.0);for(int k=0;k<n;k++)r*=x;return r;}
bool variables(const Box&b,array<I,4>&Y){for(int i=0;i<3;i++)Y[i]=I(b.l[i],b.u[i]);I y3=I(1.0)-Y[0]-Y[1]-Y[2];if(upper(y3)<0)return false;Y[3]=intersect(y3,I(0.0,1.0));return lower(Y[3])<=upper(Y[3]);}
double monUB(const Term&t,const array<I,4>&Y){I prod(1.0);double inactiveL=0;for(int i=0;i<4;i++){if(t.e[i])prod*=powi(sqrt(Y[i]),t.e[i]);else inactiveL+=lower(Y[i]);}double box=upper(prod);double S=max(0.0,1.0-inactiveL);I am(1.0);for(int i=0;i<4;i++)if(t.e[i]){I z(S*(double)t.e[i]/4.0);am*=powi(sqrt(z),t.e[i]);}return min(box,upper(am));}
double monLB(const Term&t,const array<I,4>&Y){I prod(1.0);for(int i=0;i<4;i++)if(t.e[i])prod*=powi(sqrt(Y[i]),t.e[i]);return lower(prod);}
double UB(const Box&b){array<I,4>Y;if(!variables(b,Y))return -INFINITY;I C(0),H(0),Q(0);double HL=0;for(auto&t:T){double mu=monUB(t,Y),ml=monLB(t,Y);C+=I(t.cL,t.cU)*I(0,mu);H+=I(t.hL,t.hU)*I(0,mu);Q+=I(0,t.qU)*I(0,mu);HL += t.hL*ml;}I b0=I(4.0)/I(3.7183448981203875)*sqrt(Q);double ans=upper(b0);if(HL>0){I b2=I(4.0)*C/(I(3.7183448981203875)*sqrt(I(HL)));ans=min(ans,upper(b2));}return ans;}
bool insideR(const Box&b){for(int i=0;i<3;i++)if(b.l[i]<RL[i]||b.u[i]>RU[i])return false;return true;}
int main(){const double best=0.71753326806113;Box r{{0.02,0.15,0.35},{0.12,0.30,0.55},0};r.ub=UB(r);priority_queue<Box,vector<Box>,Cmp>pq;pq.push(r);long long n=0,skip=0;double skipub=0;while(!pq.empty()){Box b=pq.top();pq.pop();if(b.ub<=best){cout<<setprecision(17)<<"CERTIFIED_OUTSIDE nodes "<<n<<" upper "<<b.ub<<" skipped "<<skip<<" skipped_ub "<<skipub<<"\n";return 0;}if(insideR(b)){skip++;skipub=max(skipub,b.ub);continue;}int k=0;double score=-1,mid=0;for(int i=0;i<3;i++){double w=b.u[i]-b.l[i];bool cross=(b.l[i]<RL[i]&&RL[i]<b.u[i])||(b.l[i]<RU[i]&&RU[i]<b.u[i]);double s=w*(cross?10:1);if(s>score){score=s;k=i;}}mid=.5*(b.l[k]+b.u[k]);if(b.l[k]<RL[k]&&RL[k]<b.u[k])mid=RL[k];else if(b.l[k]<RU[k]&&RU[k]<b.u[k])mid=RU[k];for(int side=0;side<2;side++){Box q=b;if(!side)q.u[k]=mid;else q.l[k]=mid;q.ub=UB(q);if(q.ub>best)pq.push(q);}n++;if(n>5000000){cout<<"INCOMPLETE "<<n<<" top "<<(pq.empty()?best:pq.top().ub)<<" heap "<<pq.size()<<"\n";return 2;}}
cout<<"CERTIFIED_OUTSIDE_EMPTY\n";}
