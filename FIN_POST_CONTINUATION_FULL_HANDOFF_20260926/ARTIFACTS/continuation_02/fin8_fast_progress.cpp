#include <bits/stdc++.h>
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

unordered_map<uint32_t,uint32_t> ccache;
unordered_map<uint32_t,double> ncache;
unordered_map<uint64_t,double> memo;
uint64_t splitcount=0, calls=0;
double alpha0;

uint32_t canon_code(uint32_t x){auto it=ccache.find(x); if(it!=ccache.end()) return it->second; Arr a=dec(x),t{}; uint32_t best=UINT32_MAX; for(int g=0;g<8;g++){t.fill(0); for(int i=0;i<8;i++)t[P[g][i]]=a[i]; best=min(best,enc(t));} ccache.emplace(x,best); return best;}
inline uint32_t canon_arr(const Arr&a){return canon_code(enc(a));}
double norm_code(uint32_t x){auto it=ncache.find(x); if(it!=ncache.end())return it->second;Arr a=dec(x);double t[7]={0};for(int i=0;i<8;i++) if(a[i]) for(int d=0;d<7;d++) t[d]+=a[i]*V[i][d];double s=0;for(double z:t)s+=z*z;ncache.emplace(x,s);return s;}
inline uint64_t mkey(uint32_t s,uint16_t q){return ((uint64_t)q<<32)|s;}
double logZ_code(uint32_t cs,uint16_t qmult);
inline double logZ_arr(const Arr&a,uint16_t qmult){return logZ_code(canon_arr(a),qmult);}

struct AccCtx { const Arr* c; int m,h; uint16_t q; double beta,nc,acc; Arr a;};
void recsplit(AccCtx &ctx,int idx,int rem){
 if(idx==8){if(rem) return; Arr b{};for(int i=0;i<8;i++)b[i]=(*ctx.c)[i]-ctx.a[i]; if(!lexle(ctx.a,b)) return; splitcount++; if((splitcount%10000000ull)==0){cerr<<"PROGRESS splits "<<splitcount<<" calls "<<calls<<" memo "<<memo.size()<<" canon "<<ccache.size()<<"\n"<<flush;}
   uint32_t ra=enc(ctx.a), rb=enc(b); double de=(2*norm_code(ra)+2*norm_code(rb)-ctx.nc)/ctx.m;
   double la=logZ_code(canon_code(ra),ctx.q), term;
   if(ctx.a!=b){ double lb=logZ_code(canon_code(rb),ctx.q); term=-ctx.beta*de+la+lb; }
   else { double l2=logZ_code(canon_code(ra),ctx.q*2); term=-ctx.beta*de+lse2(2*la,l2)-log(2.0); }
   ctx.acc=lse2(ctx.acc,term); return; }
 int tail=0;for(int i=idx+1;i<8;i++)tail+=(*ctx.c)[i];int lo=max(0,rem-tail),hi=min((int)(*ctx.c)[idx],rem);
 for(int x=lo;x<=hi;x++){ctx.a[idx]=x;recsplit(ctx,idx+1,rem-x);}ctx.a[idx]=0;
}

double logZ_code(uint32_t cs,uint16_t qmult){
 uint64_t key=mkey(cs,qmult);auto it=memo.find(key);if(it!=memo.end())return it->second;calls++;Arr c=dec(cs);int m=total(c);if(m==1){memo.emplace(key,0.0);return 0.0;}double beta=8.0*alpha0*qmult/m;AccCtx ctx{&c,m,m/2,qmult,beta,norm_code(cs),-INFINITY,{}};recsplit(ctx,0,m/2);memo.emplace(key,ctx.acc);return ctx.acc;
}
int main(int argc,char**argv){int n=argc>1?stoi(argv[1]):4;alpha0=argc>2?stod(argv[2]):0.745;ccache.reserve(3000000);ncache.reserve(3000000);memo.reserve(1000000);ccache.max_load_factor(.7);ncache.max_load_factor(.7);memo.max_load_factor(.7);Arr root{};root.fill(n);auto t0=chrono::steady_clock::now();double z=logZ_arr(root,1);double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();cerr<<"states "<<memo.size()<<" canon "<<ccache.size()<<" norm "<<ncache.size()<<" splits "<<splitcount<<" sec "<<sec<<"\n";cout<<setprecision(17)<<"n "<<n<<" alpha "<<alpha0<<" logZ "<<z<<"\n";}
