#include <bits/stdc++.h>
extern "C" {
typedef double fftw_complex[2];
typedef struct fftw_plan_s *fftw_plan;
double *fftw_alloc_real(size_t);
void fftw_free(void *);
fftw_plan fftw_plan_dft_r2c(int,const int*,double*,fftw_complex*,unsigned);
fftw_plan fftw_plan_dft_c2r(int,const int*,fftw_complex*,double*,unsigned);
void fftw_execute(const fftw_plan);
void fftw_destroy_plan(fftw_plan);
int fftw_init_threads(void);
void fftw_plan_with_nthreads(int);
void fftw_cleanup_threads(void);
}
#define FFTW_ESTIMATE (1U<<6)
using namespace std;
static const int D=8, N=9, PAD=10;
static const long long NGRID = 43046721LL; // 9^8
static const long long NREAL = 4782969LL*10; // 9^7*10
static const long long NCPLX = 4782969LL*5;  // 9^7*5
static const double V[8][7]={
{4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343},
{-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338},
{0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338},
{5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343},
{-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343},
{0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338},
{-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338},
{5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338}};
struct S16 { double lz,d1,lev[4],d2; bool ok=false; };
struct S32 { double logY,dY,lev[4]; };
using Arr=array<int,8>;
static inline uint32_t enc4(const Arr&a){uint32_t x=0;for(int i=0;i<8;i++)x|=(uint32_t)a[i]<<(4*i);return x;}
static Arr decHex(const string&s){unsigned long long x=stoull(s,nullptr,16);Arr a{};for(int i=0;i<8;i++)a[i]=(x>>(4*i))&15;return a;}
static double Gm[8][8];
static double normc(const Arr&a){double t[7]={};for(int i=0;i<8;i++)for(int k=0;k<7;k++)t[k]+=a[i]*V[i][k];double s=0;for(double x:t)s+=x*x;return s;}
static long long ridx(const Arr&a){long long q=0;for(int i=0;i<7;i++)q=q*9+a[i];return q*PAD+a[7];}
static long long flat9(const Arr&a){long long q=0;for(int i=0;i<8;i++)q=q*9+a[i];return q;}
static Arr unflat9(long long q){Arr a{};for(int i=7;i>=0;i--){a[i]=q%9;q/=9;}return a;}
static long long flat5(const Arr&a){long long q=0;for(int i=0;i<8;i++)q=q*5+a[i];return q;}

int main(int argc,char**argv){
  if(argc<2){cerr<<"usage: m64_depth_fftw level16_file [alpha]\n";return 2;}
  string path=argv[1]; double alpha=argc>2?stod(argv[2]):0.8245;
  for(int i=0;i<8;i++)for(int j=0;j<8;j++){double s=0;for(int k=0;k<7;k++)s+=V[i][k]*V[j][k];Gm[i][j]=s;}
  vector<pair<Arr,S16>> rows; rows.reserve(220000);
  unordered_map<uint32_t,S16> q2,q4; q2.reserve(50000); q4.reserve(100);
  ifstream f(path); string line;
  while(getline(f,line)){
    istringstream ss(line); string hs; ss>>hs; if(!ss) continue; Arr c=decHex(hs);
    vector<string> tok; string x; while(ss>>x)tok.push_back(x); if(tok.size()!=21){cerr<<"bad tokens "<<tok.size()<<"\n";return 3;}
    auto parse=[&](int off)->S16{S16 s; if(tok[off]=="nan")return s; s.ok=true;s.lz=stod(tok[off]);s.d1=stod(tok[off+1]);for(int k=0;k<4;k++)s.lev[k]=stod(tok[off+2+k]);s.d2=stod(tok[off+6]);return s;};
    S16 s1=parse(0),s2=parse(7),s4=parse(14); rows.push_back({c,s1}); if(s2.ok)q2[enc4(c)]=s2;if(s4.ok)q4[enc4(c)]=s4;
  }
  cerr<<"rows "<<rows.size()<<" q2 "<<q2.size()<<" q4 "<<q4.size()<<"\n";
  double M=-1e300;
  for(auto &rr:rows){double n=normc(rr.first);M=max(M,rr.second.lz-alpha*n/64.0);} cerr<<setprecision(17)<<"M "<<M<<"\n";
  fftw_init_threads(); fftw_plan_with_nthreads(4);
  double *base=fftw_alloc_real(NREAL), *wrk=fftw_alloc_real(NREAL); if(!base||!wrk){cerr<<"alloc fail\n";return 4;}
  fill(base,base+NREAL,0.0); fill(wrk,wrk+NREAL,0.0);
  for(auto &rr:rows){auto c=rr.first;double n=normc(c);double ly=rr.second.lz-alpha*n/64.0; base[ridx(c)]=exp(ly-M);}
  int dims[8]={9,9,9,9,9,9,9,9};
  fftw_plan pF0=fftw_plan_dft_r2c(8,dims,base,(fftw_complex*)base,FFTW_ESTIMATE);
  fftw_plan pFw=fftw_plan_dft_r2c(8,dims,wrk,(fftw_complex*)wrk,FFTW_ESTIMATE);
  fftw_plan pBw=fftw_plan_dft_c2r(8,dims,(fftw_complex*)wrk,wrk,FFTW_ESTIMATE);
  fftw_execute(pF0); fftw_complex *F0=(fftw_complex*)base;
  // enumerate target total32 compositions and maps
  vector<Arr> states; states.reserve(2400000); vector<long long> posFlat; posFlat.reserve(2400000);
  for(long long q=0;q<NGRID;q++){long long t=q;int sum=0;Arr a{};for(int i=7;i>=0;i--){a[i]=t%9;t/=9;sum+=a[i];}if(sum==32){states.push_back(a);posFlat.push_back(q);}}
  size_t NS=states.size(); cerr<<"M32 raw states "<<NS<<"\n";
  vector<double>S0(NS),S1(NS),Ssq(NS); vector<array<double,4>>Lnum(NS);
  const double invN=1.0/(double)NGRID;
  auto extract=[&](vector<double>&dst){for(size_t z=0;z<NS;z++)dst[z]=wrk[ridx(states[z])]*invN;};
  auto fill_field=[&](int kind){fill(wrk,wrk+NREAL,0.0);for(auto &rr:rows){Arr c=rr.first;double n=normc(c),ly=rr.second.lz-alpha*n/64.0,z=exp(ly-M),d=rr.second.d1-n/64.0,val=0; if(kind==0)val=z*d; else if(kind==1)val=z*d*d; else val=z*rr.second.lev[kind-2];wrk[ridx(c)]=val;}};
  auto multF0=[&](double factor){fftw_execute(pFw);fftw_complex*Fw=(fftw_complex*)wrk;for(long long i=0;i<NCPLX;i++){double ar=Fw[i][0],ai=Fw[i][1],br=F0[i][0],bi=F0[i][1];Fw[i][0]=factor*(ar*br-ai*bi);Fw[i][1]=factor*(ar*bi+ai*br);}fftw_execute(pBw);};
  // S0 = conv(A,A)
  for(long long i=0;i<NCPLX;i++){double ar=F0[i][0],ai=F0[i][1];((fftw_complex*)wrk)[i][0]=ar*ar-ai*ai;((fftw_complex*)wrk)[i][1]=2*ar*ai;}
  fftw_execute(pBw); extract(S0); cerr<<"S0 done\n";
  // S1=2 conv(A,Ad)
  fill_field(0);multF0(2.0);extract(S1);cerr<<"S1 done\n";
  // Ssq=2 conv(A,Ad2)
  fill_field(1);multF0(2.0);extract(Ssq);cerr<<"AAd2 done\n";
  // +2 conv(Ad,Ad) : FFT Ad, square
  fill_field(0);fftw_execute(pFw);{fftw_complex*Fw=(fftw_complex*)wrk;for(long long i=0;i<NCPLX;i++){double ar=Fw[i][0],ai=Fw[i][1];Fw[i][0]=2*(ar*ar-ai*ai);Fw[i][1]=4*ar*ai;}}fftw_execute(pBw);for(size_t z=0;z<NS;z++)Ssq[z]+=wrk[ridx(states[z])]*invN;cerr<<"AdAd done\n";
  // levels
  for(int k=0;k<4;k++){fill_field(2+k);multF0(2.0);for(size_t z=0;z<NS;z++)Lnum[z][k]=wrk[ridx(states[z])]*invN;cerr<<"L"<<k<<" done\n";}
  // Build M32 q1 stats and map raw flat9 -> state index by combinatorial search: use unordered map only for complements
  vector<S32> s32(NS); vector<long long> f9(NS); unordered_map<long long,size_t> mp;mp.reserve(NS*1.3);
  for(size_t z=0;z<NS;z++){Arr c=states[z];f9[z]=flat9(c);mp[f9[z]]=z;double n=normc(c);double Eq=0,Dq=0,Lq[4]={};bool even=true;Arr h{};for(int i=0;i<8;i++){if(c[i]%2)even=false;h[i]=c[i]/2;}if(even){auto it=q2.find(enc4(h));if(it==q2.end()){cerr<<"missing q2\n";return 5;}auto s=it->second;double ly=s.lz-alpha*normc(h)/32.0;Eq=exp(ly-2*M);Dq=s.d1-normc(h)/32.0;for(int k=0;k<4;k++)Lq[k]=s.lev[k];}
    double den=S0[z]+Eq; double mean=(S1[z]+Eq*Dq)/den; double root=(Ssq[z]+Eq*Dq*Dq)/den-mean*mean; s32[z].logY=log(den)+2*M-log(2.0)+alpha*n/256.0; s32[z].dY=mean+n/256.0; s32[z].lev[0]=root; for(int k=0;k<3;k++)s32[z].lev[k+1]=(Lnum[z][k]+Eq*Lq[k])/den;
  }
  cerr<<"M32 stats built\n";
  // M32 q2 homogeneous c=(4,...,4): exact ordered convolution plus q4 correction.
  // This is Z = 1/2 [sum_ordered Y2(a)Y2(b) + Y4(2,...,2)].
  vector<double> blw,bd; vector<array<double,4>>blev; Arr root4;root4.fill(4);
  long long total5=1;for(int i=0;i<8;i++)total5*=5;
  for(long long code=0;code<total5;code++){
    long long t=code;Arr a{},b{};int sm=0;
    for(int i=7;i>=0;i--){a[i]=t%5;t/=5;sm+=a[i];b[i]=4-a[i];}
    if(sm!=16)continue;
    auto ia=q2.find(enc4(a)), ib=q2.find(enc4(b));
    if(ia==q2.end()||ib==q2.end()){cerr<<"q2 child missing\n";return 6;}
    auto A=ia->second,B=ib->second;double na=normc(a),nb=normc(b);
    blw.push_back(A.lz-alpha*na/32+B.lz-alpha*nb/32);
    bd.push_back((A.d1-na/32)+(B.d1-nb/32));
    array<double,4> l{};for(int k=0;k<4;k++)l[k]=A.lev[k]+B.lev[k];blev.push_back(l);
  }
  Arr aeq;aeq.fill(2);auto it4=q4.find(enc4(aeq));if(it4==q4.end()){cerr<<"q4 missing\n";return 7;}
  auto C4=it4->second;double neq=normc(aeq);
  blw.push_back(C4.lz-alpha*neq/16);bd.push_back(C4.d1-neq/16);array<double,4>l4{};for(int k=0;k<4;k++)l4[k]=C4.lev[k];blev.push_back(l4);
  double mx=*max_element(blw.begin(),blw.end()),sw=0;for(double x:blw)sw+=exp(x-mx);double lm=mx+log(sw),md=0;for(size_t i=0;i<blw.size();i++)md+=exp(blw[i]-lm)*bd[i];
  double vr=0;array<double,4>levq{};for(size_t i=0;i<blw.size();i++){double p=exp(blw[i]-lm);vr+=p*(bd[i]-md)*(bd[i]-md);for(int k=0;k<3;k++)levq[k+1]+=p*blev[i][k];}levq[0]=vr;
  double nroot=normc(root4);double logYq=lm-log(2.0)+alpha*nroot/128.0,dYq=md+nroot/128.0;cerr<<"M32 q2 built ordered branches "<<blw.size()<<" logYq "<<logYq<<"\n";
  // M64 root mixture over unordered root splits. s32 is for ordered child partition Z representation; root recursion's unordered sum equals half ordered + equal correction.
  double maxlw=-1e300; vector<double> lw(NS); vector<size_t> cp(NS);
  Arr full;full.fill(8);
  for(size_t z=0;z<NS;z++){Arr b{};for(int i=0;i<8;i++)b[i]=8-states[z][i];auto it=mp.find(flat9(b));if(it==mp.end()){cerr<<"comp missing\n";return 8;}cp[z]=it->second;lw[z]=s32[z].logY+s32[cp[z]].logY;maxlw=max(maxlw,lw[z]);}
  // ordered root sum /2, plus equal q2 correction /2. normalization common /2 cancels probabilities.
  double sumw=0;for(double x:lw)sumw+=exp(x-maxlw);double eqw=exp(logYq-maxlw);double den=sumw+eqw;double mu=0;for(size_t z=0;z<NS;z++)mu+=exp(lw[z]-maxlw)/den*(s32[z].dY+s32[cp[z]].dY);mu+=eqw/den*dYq;
  double rootvar=0;array<double,5> lev64{};for(size_t z=0;z<NS;z++){double p=exp(lw[z]-maxlw)/den,d=s32[z].dY+s32[cp[z]].dY;rootvar+=p*(d-mu)*(d-mu);for(int k=0;k<4;k++)lev64[k+1]+=p*(s32[z].lev[k]+s32[cp[z]].lev[k]);}double pe=eqw/den;rootvar+=pe*(dYq-mu)*(dYq-mu);for(int k=0;k<4;k++)lev64[k+1]+=pe*levq[k];lev64[0]=rootvar;
  cout<<setprecision(17)<<"alpha "<<alpha<<" M64_peq "<<pe<<" d1 "<<mu<<"\n";double sumch=0;for(int k=0;k<5;k++){double ch=alpha*alpha*lev64[k];sumch+=ch;cout<<"depth "<<k<<" CH "<<ch<<"\n";}cout<<"sumCH "<<sumch<<"\n";
  // validate M32 homogeneous q1 state
  Arr h4;h4.fill(4);auto it=mp.find(flat9(h4));if(it!=mp.end()){auto s=s32[it->second];double sum=0;cout<<"M32_hom logY "<<s.logY<<" dY "<<s.dY;for(int k=0;k<4;k++){cout<<" lev"<<k<<" "<<s.lev[k];sum+=s.lev[k];}cout<<" CHsum "<<alpha*alpha*sum<<"\n";}
  fftw_destroy_plan(pF0);fftw_destroy_plan(pFw);fftw_destroy_plan(pBw);fftw_free(base);fftw_free(wrk);fftw_cleanup_threads();
}
