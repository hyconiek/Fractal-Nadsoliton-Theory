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
{5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338}};
static const int P[8][8]={{0,1,2,3,4,5,6,7},{2,3,4,5,6,7,0,1},{4,5,6,7,0,1,2,3},{6,7,0,1,2,3,4,5},{7,6,5,4,3,2,1,0},{1,0,7,6,5,4,3,2},{3,2,1,0,7,6,5,4},{5,4,3,2,1,0,7,6}};
using Arr=array<uint8_t,8>;
uint32_t enc(const Arr&a){uint32_t x=0;for(int i=0;i<8;i++)x|=(uint32_t)a[i]<<(4*i);return x;}
Arr dec(uint32_t x){Arr a{};for(int i=0;i<8;i++)a[i]=(x>>(4*i))&15;return a;}
uint32_t canon(const Arr&a){uint32_t best=UINT32_MAX;Arr t{};for(int g=0;g<8;g++){t.fill(0);for(int i=0;i<8;i++)t[P[g][i]]=a[i];best=min(best,enc(t));}return best;}
bool lexle(const Arr&a,const Arr&b){for(int i=0;i<8;i++){if(a[i]<b[i])return true;if(a[i]>b[i])return false;}return true;}
int total(const Arr&a){int s=0;for(auto x:a)s+=x;return s;}
double alpha0;
struct Stat { double lz=0,d1=0,d2=0; vector<double> lev; };
unordered_map<Key,Stat,KH> memo;
long long splitsN=0;
double deltaE(const Arr&a,const Arr&b){int na=total(a),nb=total(b);double ma[7]={0},mb[7]={0};for(int i=0;i<8;i++)for(int d=0;d<7;d++){ma[d]+=a[i]*V[i][d];mb[d]+=b[i]*V[i][d];}for(int d=0;d<7;d++){ma[d]/=na;mb[d]/=nb;}double ss=0;for(int d=0;d<7;d++){double z=ma[d]-mb[d];ss+=z*z;}return (double)na*nb/(na+nb)*ss;}
void enumRec(const Arr&c,int idx,int rem,Arr&a,vector<pair<Arr,Arr>>&out){if(idx==8){if(rem==0){Arr b{};for(int i=0;i<8;i++)b[i]=c[i]-a[i];if(lexle(a,b))out.push_back({a,b});}return;}int tail=0;for(int i=idx+1;i<8;i++)tail+=c[i];for(int x=max(0,rem-tail);x<=min((int)c[idx],rem);x++){a[idx]=x;enumRec(c,idx+1,rem-x,a,out);}a[idx]=0;}
Stat solve(uint32_t cs,uint16_t q){Key K{cs,q};auto it=memo.find(K);if(it!=memo.end())return it->second;Arr c=dec(cs);int m=total(c);if(m==1){Stat z;memo[K]=z;return z;}vector<pair<Arr,Arr>> sp;Arr aa{};enumRec(c,0,m/2,aa,sp);splitsN+=sp.size();double fac=8.0*q/m;
 struct Br{double lw,d1,d2;vector<double> lev;}; vector<Br> br;br.reserve(sp.size()+32);
 for(auto &pr:sp){auto a=pr.first,b=pr.second;double de=deltaE(a,b), rootd=-fac*de;auto A=solve(canon(a),q);
  if(a!=b){auto B=solve(canon(b),q);Br r;r.lw=-alpha0*fac*de+A.lz+B.lz;r.d1=rootd+A.d1+B.d1;r.d2=A.d2+B.d2;size_t L=max(A.lev.size(),B.lev.size());r.lev.assign(L,0);for(size_t k=0;k<A.lev.size();k++)r.lev[k]+=A.lev[k];for(size_t k=0;k<B.lev.size();k++)r.lev[k]+=B.lev[k];br.push_back(move(r));}
  else { // exact unordered correction as two separate mixture branches
    {Br r;r.lw=-alpha0*fac*de+2*A.lz-log(2.0);r.d1=rootd+2*A.d1;r.d2=2*A.d2;r.lev=A.lev;for(auto &x:r.lev)x*=2;br.push_back(move(r));}
    auto A2=solve(canon(a),q*2);{Br r;r.lw=-alpha0*fac*de+A2.lz-log(2.0);r.d1=rootd+A2.d1;r.d2=A2.d2;r.lev=A2.lev;br.push_back(move(r));}
  }
 }
 double mx=-INFINITY;for(auto&r:br)mx=max(mx,r.lw);double sw=0;for(auto&r:br)sw+=exp(r.lw-mx);double lz=mx+log(sw),mean=0;for(auto&r:br){double p=exp(r.lw-lz);mean+=p*r.d1;}
 double between=0,within=0;size_t maxL=0;for(auto&r:br){double p=exp(r.lw-lz);between+=p*(r.d1-mean)*(r.d1-mean);within+=p*r.d2;maxL=max(maxL,r.lev.size());}
 Stat out;out.lz=lz;out.d1=mean;out.d2=between+within;out.lev.assign(maxL+1,0);out.lev[0]=between;for(auto&r:br){double p=exp(r.lw-lz);for(size_t k=0;k<r.lev.size();k++)out.lev[k+1]+=p*r.lev[k];}
 double check=accumulate(out.lev.begin(),out.lev.end(),0.0); if(fabs(check-out.d2)>1e-8*max(1.0,fabs(out.d2))){cerr<<"CHECKFAIL "<<m<<" "<<q<<" "<<check<<" "<<out.d2<<"\n";exit(3);}memo[K]=out;return out;
}
int main(int argc,char**argv){int n=argc>1?stoi(argv[1]):4;alpha0=argc>2?stod(argv[2]):.745;Arr r{};r.fill(n);auto t=chrono::steady_clock::now();auto S=solve(canon(r),1);double sec=chrono::duration<double>(chrono::steady_clock::now()-t).count();cout<<setprecision(17)<<"n "<<n<<" alpha "<<alpha0<<" logZ "<<S.lz<<" d2 "<<S.d2<<" CH "<<alpha0*alpha0*S.d2<<" states "<<memo.size()<<" splits "<<splitsN<<" sec "<<sec<<"\n";double sum=0;for(size_t k=0;k<S.lev.size();k++){double ch=alpha0*alpha0*S.lev[k];sum+=ch;cout<<"depth "<<k<<" CH "<<ch<<" fraction "<<ch/(alpha0*alpha0*S.d2)<<"\n";}cout<<"sumCH "<<sum<<" residual "<<sum-alpha0*alpha0*S.d2<<"\n";}
