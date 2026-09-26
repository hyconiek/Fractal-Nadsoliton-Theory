#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std; using boost::multiprecision::cpp_int; using Arr=array<uint8_t,8>;
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
inline uint32_t enc(const Arr&a){uint32_t x=0;for(int i=0;i<8;i++)x|=(uint32_t)a[i]<<(4*i);return x;} inline Arr dec(uint32_t x){Arr a{};for(int i=0;i<8;i++)a[i]=(x>>(4*i))&15;return a;} inline int total(const Arr&a){int s=0;for(auto x:a)s+=x;return s;} inline bool lexle(const Arr&a,const Arr&b){for(int i=0;i<8;i++){if(a[i]<b[i])return true;if(a[i]>b[i])return false;}return true;}
unordered_map<uint32_t,uint32_t> CC; unordered_map<uint32_t,double> NC;
uint32_t canon(uint32_t x){auto it=CC.find(x);if(it!=CC.end())return it->second;Arr a=dec(x),t{};uint32_t best=UINT32_MAX;for(int g=0;g<8;g++){t.fill(0);for(int i=0;i<8;i++)t[P[g][i]]=a[i];best=min(best,enc(t));}CC[x]=best;return best;}
double normc(uint32_t x){auto it=NC.find(x);if(it!=NC.end())return it->second;Arr a=dec(x);double t[7]={0};for(int i=0;i<8;i++)for(int d=0;d<7;d++)t[d]+=a[i]*V[i][d];double s=0;for(double z:t)s+=z*z;NC[x]=s;return s;}
struct R{double e;cpp_int g;}; unordered_map<uint32_t,R> MM;
R solve(uint32_t cs);
void rec(const Arr&c,int idx,int rem,Arr&a,int m,double nc,double &best,cpp_int&deg){if(idx==8){if(rem)return;Arr b{};for(int i=0;i<8;i++)b[i]=c[i]-a[i];if(!lexle(a,b))return;uint32_t ra=enc(a),rb=enc(b);double de=(2*normc(ra)+2*normc(rb)-nc)/m;R A=solve(canon(ra)),B=solve(canon(rb));double E=de+2*(A.e+B.e);cpp_int g; if(a==b) g=A.g*(A.g+1)/2; else g=A.g*B.g; if(E<best-1e-10){best=E;deg=g;}else if(fabs(E-best)<=1e-10)deg+=g;return;}int tail=0;for(int i=idx+1;i<8;i++)tail+=c[i];for(int x=max(0,rem-tail);x<=min((int)c[idx],rem);x++){a[idx]=x;rec(c,idx+1,rem-x,a,m,nc,best,deg);}a[idx]=0;}
R solve(uint32_t cs){auto it=MM.find(cs);if(it!=MM.end())return it->second;Arr c=dec(cs);int m=total(c);if(m==1)return MM[cs]={0,1};double best=1e300;cpp_int deg=0;Arr a{};rec(c,0,m/2,a,m,normc(cs),best,deg);return MM[cs]={best,deg};}
int main(){for(int n:{8}){MM.clear();CC.clear();NC.clear();Arr r{};r.fill(n);auto x=solve(canon(enc(r)));cout<<n<<" "<<setprecision(17)<<x.e<<" "<<x.g<<" states "<<MM.size()<<"\n";}}
