#include <bits/stdc++.h>
using namespace std;
struct Term{array<int,4> e; double a;};
static const Term T[]={
 {{{2,2,0,0}},0.011177246149602943},
 {{{2,0,2,0}},0.011830111953407213},
 {{{1,2,1,0}},0.024378462887655366},
 {{{1,1,1,1}},0.0711429240695726},
 {{{1,0,3,0}},0.012711986408918444},
 {{{0,2,2,0}},0.01298188077452917},
 {{{0,2,0,2}},0.026354354660391023},
 {{{0,1,2,1}},0.05636892066749864},
 {{{0,0,4,0}},0.003347642998043608},
 {{{0,0,2,2}},0.027118888219793382}
};
struct B{array<double,3> l,u;};
bool derive(const B&b,array<double,4>&l,array<double,4>&u){double sl=0,su=0;for(int i=0;i<3;i++){l[i]=b.l[i];u[i]=b.u[i];sl+=l[i]*l[i];su+=u[i]*u[i];}if(sl>1)return false;l[3]=sqrt(max(0.0,1-su));u[3]=sqrt(max(0.0,1-sl));return l[3]<=u[3];}
double mono(const array<double,4>&x,const array<int,4>&e){double z=1;for(int i=0;i<4;i++)for(int k=0;k<e[i];k++)z*=x[i];return z;}
double UB(const B&b){array<double,4>l,u;if(!derive(b,l,u))return -1e300;double s=0;for(auto&t:T){double box=mono(u,t.e);double inactive=0;for(int i=0;i<4;i++)if(!t.e[i])inactive+=l[i]*l[i];double S=max(0.0,1-inactive),am=1;for(int i=0;i<4;i++)if(t.e[i]){double yi=S*(t.e[i]/4.0); am*=pow(yi,t.e[i]/2.0);}s+=t.a*min(box,am);}return s;}
pair<double,double> Hrange(const array<double,4>&l,const array<double,4>&u,int i,int j){double lo=0,hi=0;for(auto&t:T){ // c_j*d_i - c_i*d_j
  if(t.e[i]){auto e=t.e;e[i]--;e[j]++;double a=t.a*t.e[i];double mn=mono(l,e),mx=mono(u,e);lo+=a*mn;hi+=a*mx;}
  if(t.e[j]){auto e=t.e;e[j]--;e[i]++;double a=-t.a*t.e[j];double mn=mono(l,e),mx=mono(u,e);lo+=a*mx;hi+=a*mn;}
 }return {lo,hi};}
int main(int argc,char**argv){double width=argc>1?atof(argv[1]):0.002;double best=0.0183907753866966;vector<B> st(1,B{{0,0,0},{1,1,1}});long long nodes=0,leaves=0,prub=0,prk=0;array<double,4> GL={1,1,1,1},GU={0,0,0,0}; vector<B> keep;while(!st.empty()){B b=st.back();st.pop_back();nodes++;if(UB(b)<best){prub++;continue;}array<double,4>l,u;if(!derive(b,l,u))continue;bool fail=false;for(int j=1;j<4;j++){auto r=Hrange(l,u,0,j);if(r.first>0||r.second<0){fail=true;break;}}if(fail){prk++;continue;}double w=0;int k=0;for(int i=0;i<3;i++)if(b.u[i]-b.l[i]>w){w=b.u[i]-b.l[i];k=i;}if(w<=width){leaves++; keep.push_back(b); for(int i=0;i<4;i++){GL[i]=min(GL[i],l[i]);GU[i]=max(GU[i],u[i]);}continue;}double m=(b.l[k]+b.u[k])*.5;B a=b,c=b;a.u[k]=m;c.l[k]=m;st.push_back(a);st.push_back(c);}cout<<setprecision(17);cout<<"nodes "<<nodes<<" leaves "<<leaves<<" prUB "<<prub<<" prKKT "<<prk<<"\n";for(int i=0;i<4;i++)cout<<"c"<<i+3<<" ["<<GL[i]<<", "<<GU[i]<<"]\n";
 vector<int> par(keep.size()); iota(par.begin(),par.end(),0); function<int(int)> fd=[&](int x){return par[x]==x?x:par[x]=fd(par[x]);}; auto un=[&](int a,int b){a=fd(a);b=fd(b);if(a!=b)par[b]=a;}; double eps=1e-15;
 for(size_t a=0;a<keep.size();a++)for(size_t b=a+1;b<keep.size();b++){bool adj=true;for(int d=0;d<3;d++)if(keep[a].u[d]+eps<keep[b].l[d]||keep[b].u[d]+eps<keep[a].l[d]){adj=false;break;}if(adj)un(a,b);} set<int> cc; for(size_t a=0;a<keep.size();a++)cc.insert(fd(a)); cout<<"components "<<cc.size()<<"\n";
}

