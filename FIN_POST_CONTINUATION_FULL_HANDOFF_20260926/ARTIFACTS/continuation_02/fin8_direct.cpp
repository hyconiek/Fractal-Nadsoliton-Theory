#include <bits/stdc++.h>
using namespace std;
struct Key { uint32_t s; uint16_t q; bool operator==(Key const&o)const{return s==o.s&&q==o.q;} };
struct KH { size_t operator()(Key const&k) const { return ((uint64_t)k.s<<16)^k.q; } };
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
static const int P[8][8]={
{0,1,2,3,4,5,6,7},{2,3,4,5,6,7,0,1},{4,5,6,7,0,1,2,3},{6,7,0,1,2,3,4,5},
{7,6,5,4,3,2,1,0},{1,0,7,6,5,4,3,2},{3,2,1,0,7,6,5,4},{5,4,3,2,1,0,7,6}
};
using Arr=array<uint8_t,8>;
uint32_t enc(const Arr&a){uint32_t x=0;for(int i=0;i<8;i++)x|=(uint32_t)a[i]<<(4*i);return x;}
Arr dec(uint32_t x){Arr a{};for(int i=0;i<8;i++)a[i]=(x>>(4*i))&15;return a;}
uint32_t canon(const Arr&a){uint32_t best=UINT32_MAX;Arr t{};for(int g=0;g<8;g++){t.fill(0);for(int i=0;i<8;i++)t[P[g][i]]=a[i];best=min(best,enc(t));}return best;}
bool lexle(const Arr&a,const Arr&b){for(int i=0;i<8;i++){if(a[i]<b[i])return true;if(a[i]>b[i])return false;}return true;}
int total(const Arr&a){int s=0;for(auto x:a)s+=x;return s;}
double alpha0;
unordered_map<Key,double,KH> memo;
long long calls=0, splitcount=0;

double lse2(double a,double b){ if(!isfinite(a))return b;if(!isfinite(b))return a; double m=max(a,b);return m+log(exp(a-m)+exp(b-m));}
double deltaE(const Arr&a,const Arr&b){int na=total(a),nb=total(b);double ma[7]={0},mb[7]={0};for(int i=0;i<8;i++)for(int d=0;d<7;d++){ma[d]+=a[i]*V[i][d];mb[d]+=b[i]*V[i][d];}for(int d=0;d<7;d++){ma[d]/=na;mb[d]/=nb;}double ss=0;for(int d=0;d<7;d++){double z=ma[d]-mb[d];ss+=z*z;}return (double)na*nb/(na+nb)*ss;}

double logZ_key(uint32_t cs,uint16_t qmult);
double logZ_arr(const Arr&c,uint16_t qmult){return logZ_key(canon(c),qmult);} 

void enumSplitsRec(const Arr&c,int idx,int rem,Arr&a, vector<pair<Arr,Arr>>& out){
 if(idx==8){if(rem==0){Arr b{};for(int i=0;i<8;i++)b[i]=c[i]-a[i]; if(lexle(a,b)) out.push_back({a,b});}return;}
 int tail=0;for(int i=idx+1;i<8;i++)tail+=c[i];int lo=max(0,rem-tail), hi=min((int)c[idx],rem);
 for(int x=lo;x<=hi;x++){a[idx]=x;enumSplitsRec(c,idx+1,rem-x,a,out);}a[idx]=0;
}

double logZ_key(uint32_t cs,uint16_t qmult){
 Key K{cs,qmult};auto it=memo.find(K);if(it!=memo.end())return it->second;calls++;
 Arr c=dec(cs);int m=total(c); if(m==1)return memo[K]=0.0;
 int h=m/2; vector<pair<Arr,Arr>> sp; sp.reserve(1024); Arr a{};enumSplitsRec(c,0,h,a,sp); splitcount+=sp.size();
 double beta=8.0*alpha0*qmult/m; double acc=-INFINITY;
 for(auto &pr:sp){auto &aa=pr.first;auto &bb=pr.second;double d=deltaE(aa,bb);double la=logZ_arr(aa,qmult);double term;
   if(aa!=bb){double lb=logZ_arr(bb,qmult);term=-beta*d+la+lb;}
   else {double l2=logZ_arr(aa,qmult*2);term=-beta*d+lse2(2*la,l2)-log(2.0);}
   acc=lse2(acc,term);
 }
 memo[K]=acc;return acc;
}
int main(int argc,char**argv){int n=argc>1?stoi(argv[1]):4;alpha0=argc>2?stod(argv[2]):0.745;Arr root{};root.fill(n);auto t0=chrono::steady_clock::now();double z=logZ_arr(root,1);auto t1=chrono::steady_clock::now();cerr<<"states "<<memo.size()<<" calls "<<calls<<" split_terms "<<splitcount<<" sec "<<chrono::duration<double>(t1-t0).count()<<"\n";cout<<setprecision(17)<<"n "<<n<<" alpha "<<alpha0<<" logZ "<<z<<"\n";
}
