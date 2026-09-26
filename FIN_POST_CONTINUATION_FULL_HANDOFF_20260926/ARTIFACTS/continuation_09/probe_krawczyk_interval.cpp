#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/utility.hpp>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace std; using boost::numeric::interval; using boost::numeric::lower; using boost::numeric::upper;
using I=interval<double>;
struct Term{array<int,4> e; double lo,hi;};
static double ipow(double x,int p){double r=1;for(int k=0;k<p;k++)r*=x;return r;}
static I ipow(I x,int p){I r(1.0);for(int k=0;k<p;k++)r*=x;return r;}
static bool inv4(double A[4][4],double R[4][4]){double M[4][8]{};for(int i=0;i<4;i++){for(int j=0;j<4;j++)M[i][j]=A[i][j];M[i][4+i]=1;}for(int c=0;c<4;c++){int p=c;for(int r=c+1;r<4;r++)if(fabs(M[r][c])>fabs(M[p][c]))p=r;if(fabs(M[p][c])<1e-18)return false;for(int j=0;j<8;j++)swap(M[c][j],M[p][j]);double d=M[c][c];for(int j=0;j<8;j++)M[c][j]/=d;for(int r=0;r<4;r++)if(r!=c){double q=M[r][c];for(int j=0;j<8;j++)M[r][j]-=q*M[c][j];}}for(int i=0;i<4;i++)for(int j=0;j<4;j++)R[i][j]=M[i][4+j];return true;}
struct Model{vector<Term>T; array<double,4>x0,lo,hi; string name;};
static I coeff(const Term&t){return I(t.lo,t.hi);} 
static I mono(const array<I,4>&x,const array<int,4>&e){I r(1.0);for(int i=0;i<4;i++)r*=ipow(x[i],e[i]);return r;}
static I grad_i(const Model&m,const array<I,4>&x,int i){I s(0.0);for(auto&t:m.T)if(t.e[i]){auto e=t.e;e[i]--;s+=coeff(t)*double(t.e[i])*mono(x,e);}return s;}
static I hess_ij(const Model&m,const array<I,4>&x,int i,int j){I s(0.0);for(auto&t:m.T){auto e=t.e;int fac=0;if(i==j){if(e[i]>=2){fac=e[i]*(e[i]-1);e[i]-=2;}}else{if(e[i]>=1&&e[j]>=1){fac=e[i]*e[j];e[i]--;e[j]--;}}if(fac)s+=coeff(t)*double(fac)*mono(x,e);}return s;}
static array<I,4> fvec(const Model&m,const array<I,4>&x){array<I,4>f;f[0]=x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]-1.0;auto g0=grad_i(m,x,0);for(int r=1;r<4;r++){auto gj=grad_i(m,x,r);f[r]=x[r]*g0-x[0]*gj;}return f;}
static array<array<I,4>,4> jac(const Model&m,const array<I,4>&x){array<array<I,4>,4>J;for(int k=0;k<4;k++)J[0][k]=2.0*x[k];auto g0=grad_i(m,x,0);array<I,4>gg;for(int j=1;j<4;j++)gg[j]=grad_i(m,x,j);for(int r=1;r<4;r++){int j=r;for(int k=0;k<4;k++){I v(0.0);if(k==j)v+=g0;v+=x[j]*hess_ij(m,x,0,k);if(k==0)v-=gg[j];v-=x[0]*hess_ij(m,x,j,k);J[r][k]=v;}}return J;}
static void run(const Model&m){array<I,4>X,X0;for(int i=0;i<4;i++){X[i]=I(m.lo[i],m.hi[i]);X0[i]=I(m.x0[i]);}auto J0i=jac(m,X0);double J0[4][4],C[4][4];for(int i=0;i<4;i++)for(int j=0;j<4;j++)J0[i][j]=(lower(J0i[i][j])+upper(J0i[i][j]))/2;if(!inv4(J0,C)){cerr<<"singular\n";return;}auto f0=fvec(m,X0);array<I,4>Y;for(int i=0;i<4;i++){I s(m.x0[i]);for(int k=0;k<4;k++)s-=C[i][k]*f0[k];Y[i]=s;}auto JX=jac(m,X);array<array<I,4>,4>R;for(int i=0;i<4;i++)for(int j=0;j<4;j++){I s(i==j?1.0:0.0);for(int k=0;k<4;k++)s-=C[i][k]*JX[k][j];R[i][j]=s;}array<I,4>K;for(int i=0;i<4;i++){I s=Y[i];for(int j=0;j<4;j++)s+=R[i][j]*(X[j]-m.x0[j]);K[i]=s;}
 cout<<m.name<<"\n"<<setprecision(17);bool ok=true;for(int i=0;i<4;i++){cout<<"K"<<i<<" ["<<lower(K[i])<<", "<<upper(K[i])<<"] box ["<<m.lo[i]<<", "<<m.hi[i]<<"]\n";if(!(lower(K[i])>m.lo[i]&&upper(K[i])<m.hi[i]))ok=false;}cout<<"STRICT_INCLUSION "<<(ok?"YES":"NO")<<"\n";}
int main(){
 vector<array<int,4>> E={{{2,2,0,0}},{{2,0,2,0}},{{1,2,1,0}},{{1,1,1,1}},{{1,0,3,0}},{{0,2,2,0}},{{0,2,0,2}},{{0,1,2,1}},{{0,0,4,0}},{{0,0,2,2}}};
 vector<pair<double,double>> A={{0.04822142787987651,0.04822142787987652},{0.05333608476307319,0.0533360847630732},{0.11385707390148077,0.1138570739014808},{0.34286799970231563,0.34286799970231574},{0.06204319383393707,0.06204319383393708},{0.06563564602160425,0.06563564602160428},{0.13577211811522816,0.1357721181152282},{0.2940921349879692,0.2940921349879693},{0.017687575726200232,0.01768757572620024},{0.14600143188162154,0.1460014318816216}};
 Model c; c.name="coefficient_sphere";for(int i=0;i<10;i++)c.T.push_back({E[i],A[i].first,A[i].second});c.x0={0.2932238705,0.4863537663,0.6521040063,0.5022351449};double n=0;for(double x:c.x0)n+=x*x;n=sqrt(n);for(double&x:c.x0)x/=n;c.lo={0.29248046875,0.485107421875,0.650634765625,0.50070260273916145};c.hi={0.2939453125,0.487548828125,0.653564453125,0.50375470388786237};run(c);
 vector<pair<double,double>> F={{0.011177246149602938,0.011177246149602941},{0.011830111953407212,0.011830111953407215},{0.024378462887655362,0.02437846288765537},{0.07114292406957255,0.07114292406957258},{0.01271198640891844,0.012711986408918444},{0.012981880774529162,0.012981880774529166},{0.026354354660391013,0.02635435466039102},{0.056368920667498615,0.05636892066749863},{0.003347642998043607,0.003347642998043608},{0.027118888219793375,0.027118888219793382}};
 Model f;f.name="fisher_sphere";for(int i=0;i<10;i++)f.T.push_back({E[i],F[i].first,F[i].second});f.x0={0.3150280854,0.4923890678,0.6461252850,0.4907468058};n=0;for(double x:f.x0)n+=x*x;n=sqrt(n);for(double&x:f.x0)x/=n;f.lo={0.314453125,0.490966796875,0.644775390625,0.48899325529568444};f.hi={0.315673828125,0.493896484375,0.6474609375,0.49221644240765933};run(f);
}
