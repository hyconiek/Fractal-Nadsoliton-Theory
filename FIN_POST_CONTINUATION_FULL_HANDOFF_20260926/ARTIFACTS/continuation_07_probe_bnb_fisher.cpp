#include <bits/stdc++.h>
using namespace std;
struct Term { array<double,4> p; double c; };
struct Box { array<double,4> l,u; double ub; long long id; };
struct Cmp { bool operator()(Box const&a, Box const&b) const { return a.ub < b.ub; } };
static vector<Term> terms={
 {{{0,0,1,1}},0.027118888219793386},
 {{{0,0,2,0}},0.003347642998043605},
 {{{0,.5,1,.5}},0.056368920667498615},
 {{{0,1,0,1}},0.02635435466039104},
 {{{0,1,1,0}},0.012981880774529166},
 {{{.5,0,1.5,0}},0.01271198640891842},
 {{{.5,.5,.5,.5}},0.07114292406957254},
 {{{.5,1,.5,0}},0.024378462887655362},
 {{{1,0,1,0}},0.011830111953407206},
 {{{1,1,0,0}},0.011177246149602929}
};
static bool tighten(array<double,4>&l,array<double,4>&u){
 for(int it=0;it<8;it++){
  bool ch=false; double sl=accumulate(l.begin(),l.end(),0.0), su=accumulate(u.begin(),u.end(),0.0);
  if(sl>1+1e-14||su<1-1e-14) return false;
  for(int i=0;i<4;i++){
   double nl=max(l[i],1-(su-u[i]));
   double nu=min(u[i],1-(sl-l[i]));
   nl=max(0.0,nl);nu=min(1.0,nu);
   if(nl>nu+1e-14) return false;
   if(abs(nl-l[i])>1e-15||abs(nu-u[i])>1e-15)ch=true;
   l[i]=nl;u[i]=nu;
  }
  if(!ch) break;
 }
 return true;
}
static double monoMax(Term const&t,array<double,4> const&l,array<double,4> const&u){
 vector<int>a,ina; double linact=0, uact=0;
 for(int i=0;i<4;i++) if(t.p[i]>0){a.push_back(i);uact+=u[i];} else {ina.push_back(i);linact+=l[i];}
 if(a.empty()) return 1.0;
 double S=min(uact,1-linact);
 double minS=0;for(int i:a)minS+=l[i];
 if(S<minS-1e-14) return 0;
 // solve sum clip(p/lambda,l,u)=S. Handle S close upper/lower.
 double sumU=0,sumL=0;for(int i:a){sumU+=u[i];sumL+=l[i];}
 array<double,4> y{};
 if(S>=sumU-1e-14){for(int i:a)y[i]=u[i];}
 else if(S<=sumL+1e-14){for(int i:a)y[i]=l[i];}
 else {
   double lo=1e-14,hi=1e14;
   for(int it=0;it<70;it++){
    double lam=sqrt(lo*hi), ss=0;
    for(int i:a) ss+=min(u[i],max(l[i],t.p[i]/lam));
    if(ss>S) lo=lam; else hi=lam;
   }
   for(int i:a)y[i]=min(u[i],max(l[i],t.p[i]/hi));
 }
 double prod=1;
 for(int i:a){ if(y[i]<=0) return 0; prod*=pow(y[i],t.p[i]); }
 return prod;
}
static double UB(array<double,4> const&l,array<double,4> const&u){double s=0;for(auto&t:terms)s+=t.c*monoMax(t,l,u);return s;}
static double P(array<double,4> const&y){double s=0;for(auto&t:terms){double m=t.c;for(int i=0;i<4;i++) if(t.p[i])m*=pow(max(y[i],0.0),t.p[i]);s+=m;}return s;}
int main(int argc,char**argv){
 double best=0.018390775386696624; double tol=argc>1?stod(argv[1]):1e-7; long long maxnodes=argc>2?stoll(argv[2]):5000000;
 array<double,4> l={0,0,0,0},u={1,1,1,1};tighten(l,u);
 priority_queue<Box,vector<Box>,Cmp> pq; long long id=0;pq.push({l,u,UB(l,u),id++});
 auto st=chrono::steady_clock::now(); long long n=0;
 while(!pq.empty()&&n<maxnodes){
  Box b=pq.top();pq.pop();
  if(b.ub<=best+tol){cerr<<"CERTIFIED n="<<n<<" ub="<<setprecision(17)<<b.ub<<" best="<<best<<" gap="<<b.ub-best<<"\n";return 0;}
  // choose widest relative interval, but prefer variables contributing candidate mass
  int k=0;double bw=-1;for(int i=0;i<4;i++){double w=b.u[i]-b.l[i];if(w>bw){bw=w;k=i;}}
  double mid=(b.l[k]+b.u[k])/2;
  for(int side=0;side<2;side++){
   auto lc=b.l,uc=b.u; if(side==0)uc[k]=mid;else lc[k]=mid;
   if(!tighten(lc,uc))continue;double ub=UB(lc,uc);if(ub>best+tol)pq.push({lc,uc,ub,id++});
  }
  n++;
  if(n%100000==0){auto sec=chrono::duration<double>(chrono::steady_clock::now()-st).count();cerr<<n<<" top="<<setprecision(12)<<(pq.empty()?best:pq.top().ub)<<" gap="<<(pq.empty()?0:pq.top().ub-best)<<" heap="<<pq.size()<<" sec="<<sec<<"\n";}
 }
 cerr<<"STOP n="<<n<<" top="<<setprecision(17)<<(pq.empty()?best:pq.top().ub)<<" gap="<<(pq.empty()?0:pq.top().ub-best)<<" heap="<<pq.size()<<"\n";
}
