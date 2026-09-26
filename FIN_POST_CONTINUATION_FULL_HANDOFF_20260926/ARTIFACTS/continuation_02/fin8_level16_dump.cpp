#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using Arr=array<uint8_t,8>;
static const double V[8][7]={
{4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343},
{-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338},
{0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338},
{5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343},
{-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343},
{0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338},
{-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338},
{5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338}
};
static const int P[8][8]={{0,1,2,3,4,5,6,7},{2,3,4,5,6,7,0,1},{4,5,6,7,0,1,2,3},{6,7,0,1,2,3,4,5},{7,6,5,4,3,2,1,0},{1,0,7,6,5,4,3,2},{3,2,1,0,7,6,5,4},{5,4,3,2,1,0,7,6}};
inline uint32_t enc(const Arr&a){uint32_t x=0;for(int i=0;i<8;i++)x|=(uint32_t)a[i]<<(4*i);return x;}
inline Arr dec(uint32_t x){Arr a{};for(int i=0;i<8;i++)a[i]=(x>>(4*i))&15;return a;}
inline int total(const Arr&a){int s=0;for(auto x:a)s+=x;return s;}
inline bool lexle(const Arr&a,const Arr&b){for(int i=0;i<8;i++){if(a[i]<b[i])return true;if(a[i]>b[i])return false;}return true;}
inline double lse2(double a,double b){if(!isfinite(a))return b;if(!isfinite(b))return a;double m=max(a,b);return m+log(exp(a-m)+exp(b-m));}
uint32_t canon_slow(const Arr&a){uint32_t best=UINT32_MAX;Arr t{};for(int g=0;g<8;g++){t.fill(0);for(int i=0;i<8;i++)t[P[g][i]]=a[i];best=min(best,enc(t));}return best;}

double alpha0;
unordered_map<uint32_t,uint32_t> CAN;
unordered_map<uint32_t,double> NORM;
unordered_map<uint64_t,double> VAL; // key q<<32 | canon code; total encoded in code
inline uint64_t key(uint32_t c,int q){return ((uint64_t)q<<32)|c;}
inline uint32_t canon_code(uint32_t x){auto it=CAN.find(x); if(it!=CAN.end()) return it->second; Arr a=dec(x); return canon_slow(a);} // fallback should not occur
inline double norm_code(uint32_t x){auto it=NORM.find(x); if(it!=NORM.end())return it->second; Arr a=dec(x);double t[7]={0};for(int i=0;i<8;i++)for(int d=0;d<7;d++)t[d]+=a[i]*V[i][d];double s=0;for(double z:t)s+=z*z;return s;}
inline double getv(uint32_t raw,int q){int m=total(dec(raw)); if(m==1)return 0.; uint32_t c=canon_code(raw); auto it=VAL.find(key(c,q)); if(it==VAL.end()){cerr<<"MISSING m="<<m<<" q="<<q<<" code="<<hex<<c<<dec<<"\n"; abort();} return it->second;}

void gen_raw_rec(int idx,int rem,int cap,Arr &a,vector<uint32_t>&out){
 if(idx==8){if(rem==0)out.push_back(enc(a));return;}
 int tail=(7-idx)*cap; int lo=max(0,rem-tail), hi=min(cap,rem);
 for(int x=lo;x<=hi;x++){a[idx]=x;gen_raw_rec(idx+1,rem-x,cap,a,out);}a[idx]=0;
}
vector<uint32_t> gen_raw(int m,int cap){vector<uint32_t> out;Arr a{};gen_raw_rec(0,m,cap,a,out);return out;}
vector<uint32_t> unique_can(const vector<uint32_t>&raw){vector<uint32_t> v;v.reserve(raw.size());for(auto x:raw)v.push_back(CAN.at(x));sort(v.begin(),v.end());v.erase(unique(v.begin(),v.end()),v.end());return v;}

struct EvalCtx{Arr c,a;int m,h,q;double beta,nc,acc;uint64_t splits;};
void split_rec(EvalCtx &ctx,int idx,int rem){
 if(idx==8){if(rem)return;Arr b{};for(int i=0;i<8;i++)b[i]=ctx.c[i]-ctx.a[i];if(!lexle(ctx.a,b))return;ctx.splits++;
   uint32_t ra=enc(ctx.a),rb=enc(b);double de=(2*norm_code(ra)+2*norm_code(rb)-ctx.nc)/ctx.m;
   double la=getv(ra,ctx.q),term;
   if(ctx.a!=b){double lb=getv(rb,ctx.q);term=-ctx.beta*de+la+lb;}
   else{double l2=getv(ra,ctx.q*2);term=-ctx.beta*de+lse2(2*la,l2)-log(2.0);}
   ctx.acc=lse2(ctx.acc,term);return;
 }
 int tail=0;for(int i=idx+1;i<8;i++)tail+=ctx.c[i];int lo=max(0,rem-tail),hi=min((int)ctx.c[idx],rem);
 for(int x=lo;x<=hi;x++){ctx.a[idx]=x;split_rec(ctx,idx+1,rem-x);}ctx.a[idx]=0;
}

double eval_state(uint32_t code,int q,uint64_t &splits){Arr c=dec(code);int m=total(c);EvalCtx ctx{c,{},m,m/2,q,8.0*alpha0*q/m,norm_code(code),-INFINITY,0};split_rec(ctx,0,m/2);splits=ctx.splits;return ctx.acc;}

int main(int argc,char**argv){alpha0=argc>1?stod(argv[1]):0.745; int threads=argc>2?stoi(argv[2]):5; omp_set_num_threads(threads);
 CAN.reserve(4000000);NORM.reserve(4000000);VAL.reserve(1000000); CAN.max_load_factor(.7);NORM.max_load_factor(.7);VAL.max_load_factor(.7);
 // Build canonical/norm cache for every raw composition at totals used, cap 8.
 for(int m: {1,2,4,8,16}){
   auto raw=gen_raw(m,8); cerr<<"raw m "<<m<<" count "<<raw.size()<<"\n";
   for(auto x:raw){Arr a=dec(x);uint32_t c=canon_slow(a);CAN.emplace(x,c); if(!NORM.count(x)){double t[7]={0};for(int i=0;i<8;i++)for(int d=0;d<7;d++)t[d]+=a[i]*V[i][d];double s=0;for(double z:t)s+=z*z;NORM.emplace(x,s);} }
 }
 cerr<<"cache built CAN "<<CAN.size()<<"\n";
 // Compute levels bottom-up. For q, cap=8/q. Only feasible if m<=8*cap.
 for(int m: {2,4,8,16}){
   for(int q: {8,4,2,1}){int cap=8/q;if(m>8*cap)continue;
     auto raw=gen_raw(m,cap); auto states=unique_can(raw); cerr<<"LEVEL m "<<m<<" q "<<q<<" cap "<<cap<<" states "<<states.size()<<" raw "<<raw.size()<<"\n";
     vector<double> vals(states.size()); vector<uint64_t> sp(states.size());
     auto t0=chrono::steady_clock::now();
     #pragma omp parallel for schedule(dynamic,32)
     for(long long i=0;i<(long long)states.size();i++) vals[i]=eval_state(states[i],q,sp[i]);
     for(size_t i=0;i<states.size();i++)VAL.emplace(key(states[i],q),vals[i]);
     uint64_t nsp=0;for(auto x:sp)nsp+=x;double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();cerr<<"DONE m "<<m<<" q "<<q<<" splits "<<nsp<<" sec "<<sec<<" VAL "<<VAL.size()<<"\n";
   }
 }
 auto raw16=gen_raw(16,8); cout<<setprecision(17);
for(auto x:raw16){ Arr a=dec(x); uint32_t c=CAN.at(x); double z1=VAL.at(key(c,1)); cout<<hex<<x<<dec<<" "<<z1; bool ok2=true,ok4=true; for(auto v:a){if(v>4)ok2=false;if(v>2)ok4=false;} if(ok2)cout<<" "<<VAL.at(key(c,2)); else cout<<" nan"; if(ok4)cout<<" "<<VAL.at(key(c,4)); else cout<<" nan"; cout<<"\n";} cerr<<"TOTAL VAL "<<VAL.size()<<"\n";
}
