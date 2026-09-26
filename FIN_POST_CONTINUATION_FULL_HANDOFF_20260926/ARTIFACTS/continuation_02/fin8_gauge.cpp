#include <bits/stdc++.h>
using namespace std;
struct Key { uint32_t s; uint64_t betaBits; bool operator==(Key const&o)const{return s==o.s&&betaBits==o.betaBits;} };
struct KH { size_t operator()(Key const&k)const{return ((uint64_t)k.s*11400714819323198485ull)^k.betaBits;} };
using Arr=array<uint8_t,8>;
static const double V[8][7]={
{4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343},
{-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338},
{0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338},
{5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343},
{-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343},
{0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338},
{-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338},
{5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338,0}
};
// NOTE: last row above deliberately not used; replace below with runtime constants for row 7.
static const double V7[7]={5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338};
static const int P[8][8]={{0,1,2,3,4,5,6,7},{2,3,4,5,6,7,0,1},{4,5,6,7,0,1,2,3},{6,7,0,1,2,3,4,5},{7,6,5,4,3,2,1,0},{1,0,7,6,5,4,3,2},{3,2,1,0,7,6,5,4},{5,4,3,2,1,0,7,6}};
uint32_t enc(const Arr&a){uint32_t x=0;for(int i=0;i<8;i++)x|=(uint32_t)a[i]<<(4*i);return x;}
Arr dec(uint32_t x){Arr a{};for(int i=0;i<8;i++)a[i]=(x>>(4*i))&15;return a;}
uint32_t canon(const Arr&a){uint32_t best=UINT32_MAX;Arr t{};for(int g=0;g<8;g++){t.fill(0);for(int i=0;i<8;i++)t[P[g][i]]=a[i];best=min(best,enc(t));}return best;}
int total(const Arr&a){int s=0;for(auto x:a)s+=x;return s;}
double vv(int i,int d){return i==7?V7[d]:V[i][d];}
double normT2(const Arr&c){double t[7]={0};for(int i=0;i<8;i++)for(int d=0;d<7;d++)t[d]+=c[i]*vv(i,d);double s=0;for(double x:t)s+=x*x;return s;}
double lse2(double a,double b){if(!isfinite(a))return b;if(!isfinite(b))return a;double m=max(a,b);return m+log(exp(a-m)+exp(b-m));}
unordered_map<Key,double,KH> memo;
uint64_t bits(double x){uint64_t b;memcpy(&b,&x,8);return b;}
double logY(const Arr&c,double beta);
void enumOrd(const Arr&c,int idx,int rem,Arr&a,double beta,double &acc){if(idx==8){if(rem==0){Arr b{};for(int i=0;i<8;i++)b[i]=c[i]-a[i];double term=logY(a,2*beta)+logY(b,2*beta);acc=lse2(acc,term);}return;}int tail=0;for(int i=idx+1;i<8;i++)tail+=c[i];for(int x=max(0,rem-tail);x<=min((int)c[idx],rem);x++){a[idx]=x;enumOrd(c,idx+1,rem-x,a,beta,acc);}a[idx]=0;}
double logY(const Arr&cin,double beta){uint32_t cs=canon(cin);Arr c=dec(cs);Key k{cs,bits(beta)};auto it=memo.find(k);if(it!=memo.end())return it->second;int m=total(c);double n2=normT2(c);if(m==1)return memo[k]=-beta*n2/(2*m);double acc=-INFINITY;Arr a{};enumOrd(c,0,m/2,a,beta,acc);bool even=true;Arr h{};for(int i=0;i<8;i++){if(c[i]&1)even=false;h[i]=c[i]/2;}if(even)acc=lse2(acc,logY(h,4*beta));double ans=beta*n2/(2*m)-log(2.0)+acc;memo[k]=ans;return ans;}
int main(int argc,char**argv){int n=argc>1?stoi(argv[1]):4;double alpha=argc>2?stod(argv[2]):0.745;Arr c{};c.fill(n);double beta=alpha/n;double z=logY(c,beta)+beta*normT2(c)/(2*total(c));cout<<setprecision(17)<<z<<" states "<<memo.size()<<"\n";}
