#include <bits/stdc++.h>
using namespace std;
constexpr int Q=12, D=7;
struct State { array<int8_t,Q> a{}; bool operator==(State const&o) const{return a==o.a;} };
struct Hash { size_t operator()(State const&s) const noexcept { uint64_t h=1469598103934665603ULL; for(auto v:s.a){h^=(uint8_t)(v+8);h*=1099511628211ULL;} return h;} };
struct Poly { double c[3]{0,0,0}; };
static inline Poly addp(Poly a,const Poly&b,double sign=1){for(int k=0;k<3;k++)a.c[k]+=sign*b.c[k];return a;}
static inline Poly mulp(const Poly&a,const Poly&b){Poly z; for(int i=0;i<3;i++)for(int j=0;j+i<3;j++)z.c[i+j]+=a.c[i]*b.c[j];return z;}

double E[Q][D], PV[Q][Q], PH[Q][Q], Aop[Q][Q], lam[7];
double THETA=1.0;
const double u=1.0/Q;

void init(){
  double W[Q][Q]{}; double L[Q][Q]{};
  for(int i=0;i<Q;i++) for(int j=0;j<Q;j++) if(i!=j){int dd=abs(i-j);dd=min(dd,Q-dd);W[i][j]=cos(0.18575*dd+0.1625)/(1+pow((double)dd,1.8));}
  for(int i=0;i<Q;i++){double s=0;for(int j=0;j<Q;j++)s+=W[i][j];for(int j=0;j<Q;j++)L[i][j]=(i==j?s:0)-W[i][j];}
  // circulant eigenvalues from first row, real DFT
  for(int k=0;k<=6;k++){double s=0;for(int j=0;j<Q;j++)s+=L[0][j]*cos(2*M_PI*k*j/Q);lam[k]=s;}
  int col=0;
  for(int k: {3,4,5}){
    for(int j=0;j<Q;j++) E[j][col]=sqrt(2.0/Q)*cos(2*M_PI*k*j/Q); col++;
    for(int j=0;j<Q;j++) E[j][col]=sqrt(2.0/Q)*sin(2*M_PI*k*j/Q); col++;
  }
  for(int j=0;j<Q;j++) E[j][col]=((j%2)?-1.0:1.0)/sqrt((double)Q);
  for(int i=0;i<Q;i++)for(int j=0;j<Q;j++){double s=0;for(int r=0;r<D;r++)s+=E[i][r]*E[j][r];PV[i][j]=s;PH[i][j]=(i==j?1.0:0.0)-1.0/Q-s;}
  double ad[D]={lam[3],lam[3],lam[4],lam[4],lam[5],lam[5],lam[6]};
  for(int i=0;i<Q;i++)for(int j=0;j<Q;j++){double s=0;for(int r=0;r<D;r++)s+=E[i][r]*ad[r]*E[j][r];Aop[i][j]=s;}
}

struct Coeffs { double v[Q],w[Q],aa[Q],maa; };
Coeffs coeffs(const State&s){
  Coeffs o{}; double d[Q];for(int i=0;i<Q;i++)d[i]=s.a[i];
  for(int i=0;i<Q;i++){for(int j=0;j<Q;j++){o.v[i]+=PV[i][j]*d[j];o.aa[i]+=Aop[i][j]*d[j];}}
  double vv[Q];for(int i=0;i<Q;i++)vv[i]=o.v[i]*o.v[i];
  for(int i=0;i<Q;i++){double h=0;for(int j=0;j<Q;j++)h+=PH[i][j]*vv[j];o.w[i]=6*h;o.maa+=o.aa[i]*o.aa[i]/Q;}
  return o;
}
Poly rate(const State&s,bool closed,int hr,int i,int j,const Coeffs&c){
  Poly z; double p1=closed?c.v[i]:(double)s.a[i]; double p2=closed?c.w[i]:0.0;
  double t[Q], mt2=0.0;
  for(int l=0;l<Q;l++){ t[l]=c.aa[l]-THETA*Aop[l][i]; mt2 += t[l]*t[l]/Q; }
  if(hr==0) z.c[0]=u*u;
  else if(hr==1){z.c[0]=p1*u; z.c[1]=u*t[j]/Q;}
  else {z.c[0]=p2*u; z.c[1]=p1*t[j]/Q; z.c[2]=u*(t[j]*t[j]-mt2)/(2*Q);}
  return z;
}
using Map=unordered_map<State,Poly,Hash>;
Map apply(const Map&m,bool closed,int hr){
  Map out; out.reserve(m.size()*8+100);
  for(auto const&kv:m){const State&s=kv.first; const Poly&wp=kv.second; Coeffs c=coeffs(s);
    for(int i=0;i<Q;i++)for(int j=0;j<Q;j++)if(i!=j){Poly rp=rate(s,closed,hr,i,j,c), cp=mulp(wp,rp); State t=s;t.a[i]--;t.a[j]++;
      auto &zt=out[t]; auto &zs=out[s]; for(int g=0;g<3;g++){zt.c[g]+=cp.c[g];zs.c[g]-=cp.c[g];}
    }
  }
  return out;
}
Map final_measure(bool closed){
  vector<array<int,3>> seq={{2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}};
  Map total; State z{};
  for(auto ss:seq){Map m; m[z].c[0]=1; for(int t=0;t<3;t++)m=apply(m,closed,ss[t]); cerr<<(closed?"C":"F")<<" seq "<<ss[0]<<ss[1]<<ss[2]<<" stepdone states="<<m.size()<<"\n"; for(auto &kv:m){auto&o=total[kv.first];for(int g=0;g<3;g++)o.c[g]+=kv.second.c[g];}}
  return total;
}

void enum_alpha_rec(int rem,int pos,array<int,D>&a,vector<array<int,D>>&out){ if(pos==D-1){a[pos]=rem;out.push_back(a);return;} for(int v=0;v<=rem;v++){a[pos]=v;enum_alpha_rec(rem-v,pos+1,a,out);} }
double facti(int n){double r=1;for(int i=2;i<=n;i++)r*=i;return r;}

int main(int argc,char**argv){ios::sync_with_stdio(false); if(argc>1)THETA=stod(argv[1]); init();cerr<<setprecision(17)<<"theta="<<THETA<<" lam5="<<lam[5]<<"\n";
  Map F=final_measure(false), C=final_measure(true), Wm=F; for(auto &kv:C){auto&o=Wm[kv.first];for(int g=0;g<3;g++)o.c[g]-=kv.second.c[g];}
  cerr<<"counts F="<<F.size()<<" C="<<C.size()<<" diff="<<Wm.size()<<"\n";
  vector<array<int,D>> alphas; array<int,D>a{};enum_alpha_rec(4,0,a,alphas); int M=alphas.size(); vector<array<double,3>> actual(M); // store per monomial [g]
  for(auto const&kv:Wm){double x[D]{};for(int r=0;r<D;r++)for(int i=0;i<Q;i++)x[r]+=E[i][r]*kv.first.a[i];
    for(int m=0;m<M;m++){double multi=facti(4),val=1;for(int r=0;r<D;r++){multi/=facti(alphas[m][r]); if(alphas[m][r])val*=pow(x[r],alphas[m][r]);} val*=multi;for(int g=0;g<3;g++)actual[m][g]+=kv.second.c[g]*val;}
  }
  // candidate via direct polynomial evaluation on enough? Instead output actual coefficients for python candidate comparison.
  cout<<setprecision(17); cout<<"# sectors lambda "<<lam[3]<<" "<<lam[4]<<" "<<lam[5]<<" "<<lam[6]<<"\n";
  for(int m=0;m<M;m++){for(int r=0;r<D;r++){if(r)cout<<',';cout<<alphas[m][r];}cout<<" "<<actual[m][0]<<" "<<actual[m][1]<<" "<<actual[m][2]<<"\n";}
}
