#include <bits/stdc++.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
extern "C" {
typedef double fftw_complex[2]; typedef struct fftw_plan_s *fftw_plan;
double *fftw_alloc_real(size_t); void fftw_free(void*);
fftw_plan fftw_plan_dft_r2c(int,const int*,double*,fftw_complex*,unsigned);
fftw_plan fftw_plan_dft_c2r(int,const int*,fftw_complex*,double*,unsigned);
void fftw_execute(const fftw_plan); void fftw_destroy_plan(fftw_plan);
int fftw_init_threads(void); void fftw_plan_with_nthreads(int); void fftw_cleanup_threads(void);
}
#define FFTW_ESTIMATE (1U<<6)
using namespace std;
static const int N=17,PAD=18,D=7;
static const long long N3=4913LL,N6=24137569LL,NGRID=410338673LL,NREAL=N6*PAD,NCPLX=N6*9;
static const double V[8][7]={
{4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343},
{-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338},
{0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338},
{5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343},
{-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343},
{0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338},
{-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338},
{5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338}};
struct BaseRow{array<int,8> c; long long ridx; double lz,d1,d2,norm;};
struct Triple{array<double,7> v; int sum; bool even; int half9;};
static vector<Triple> TT;
static inline array<int,8> dec5(const string&s){unsigned long long x=stoull(s,nullptr,16);array<int,8>a{};for(int i=0;i<8;i++)a[i]=(x>>(5*i))&31;return a;}
static inline double normc(const array<int,8>&a){double t[7]={};for(int i=0;i<8;i++)for(int k=0;k<7;k++)t[k]+=a[i]*V[i][k];double z=0;for(double x:t)z+=x*x;return z;}
static inline long long ridx7(const array<int,8>&a){long long q=0;for(int i=0;i<6;i++)q=q*17+a[i];return q*PAD+a[6];}
void write_full(const string&p,const double*a,size_t n){FILE*f=fopen(p.c_str(),"wb");if(!f){perror("fopen write");exit(2);}size_t w=fwrite(a,sizeof(double),n,f);if(w!=n){cerr<<"short write "<<p<<"\n";exit(2);}fclose(f);}
void read_full(const string&p,double*a,size_t n){FILE*f=fopen(p.c_str(),"rb");if(!f){perror("fopen read");exit(2);}size_t r=fread(a,sizeof(double),n,f);if(r!=n){cerr<<"short read "<<p<<" "<<r<<"/"<<n<<"\n";exit(2);}fclose(f);}
void complex_combine(double*w,const string&f0,const string&f1,int mode){ // mode1:2F1F0; mode2:2F2F0+2F1^2, w starts F1 or F2
 FILE*a=fopen(f0.c_str(),"rb"); if(!a){perror("f0");exit(2);} FILE*b=nullptr;if(mode==2){b=fopen(f1.c_str(),"rb");if(!b){perror("f1");exit(2);}}
 const size_t CH=1<<20; vector<double>A(2*CH),B(mode==2?2*CH:0); fftw_complex*W=(fftw_complex*)w; long long off=0;
 while(off<NCPLX){size_t n=min<long long>(CH,NCPLX-off); if(fread(A.data(),sizeof(double),2*n,a)!=2*n){cerr<<"read f0 fail\n";exit(2);} if(mode==2&&fread(B.data(),sizeof(double),2*n,b)!=2*n){cerr<<"read f1 fail\n";exit(2);} for(size_t j=0;j<n;j++){double ar=A[2*j],ai=A[2*j+1], xr=W[off+j][0],xi=W[off+j][1]; double rr=2*(ar*xr-ai*xi), ii=2*(ar*xi+ai*xr); if(mode==2){double br=B[2*j],bi=B[2*j+1];rr+=2*(br*br-bi*bi);ii+=4*br*bi;} W[off+j][0]=rr;W[off+j][1]=ii;} off+=n; }
 fclose(a);if(b)fclose(b);
}
void prep_triples(){TT.resize(N3);for(int id=0;id<N3;id++){int x=id,d[3];for(int k=2;k>=0;k--){d[k]=x%17;x/=17;}Triple t{};t.sum=d[0]+d[1]+d[2];t.even=((d[0]|d[1]|d[2])%2==0);t.half9=(d[0]/2)*81+(d[1]/2)*9+d[2]/2;for(int k=0;k<7;k++)t.v[k]=d[0]*V[0][k]+d[1]*V[1][k]+d[2]*V[2][k];TT[id]=t;}}
void build_parent(const string&conv0,const string&conv1,double*w,const string&diagpref,double alpha,int parentM,double Schild,const string&outpref){
 vector<double>d0(4782969),d1(4782969),d2(4782969); read_full(diagpref+"g0.raw",d0.data(),d0.size());read_full(diagpref+"g1.raw",d1.data(),d1.size());read_full(diagpref+"g2.raw",d2.data(),d2.size());
 FILE*f0=fopen(conv0.c_str(),"rb"),*f1=fopen(conv1.c_str(),"rb"); if(!f0||!f1){perror("conv open");exit(2);} FILE*o0=fopen((outpref+"g0.raw").c_str(),"wb"),*o1=fopen((outpref+"g1.raw").c_str(),"wb"),*o2=fopen((outpref+"g2.raw").c_str(),"wb");if(!o0||!o1||!o2){perror("out");exit(2);} 
 const long long BR=32768; vector<double>b0(BR*PAD),b1(BR*PAD),z0(BR*PAD),z1(BR*PAD),z2(BR*PAD); double inv=1.0/(double)NGRID, ds=exp(-2*Schild), ap=4.0/(parentM*(double)parentM), Dv[7];for(int k=0;k<7;k++)Dv[k]=V[6][k]-V[7][k];double DD=0;for(double x:Dv)DD+=x*x;
 long long done=0;auto st=chrono::steady_clock::now();while(done<N6){long long nr=min<long long>(BR,N6-done),cnt=nr*PAD;if(fread(b0.data(),sizeof(double),cnt,f0)!=cnt||fread(b1.data(),sizeof(double),cnt,f1)!=cnt){cerr<<"conv read short\n";exit(2);} fill(z0.begin(),z0.begin()+cnt,0);fill(z1.begin(),z1.begin()+cnt,0);fill(z2.begin(),z2.begin()+cnt,0);
  for(long long rr=0;rr<nr;rr++){long long r=done+rr;int ia=r/N3,ib=r%N3;auto&A=TT[ia];auto&B=TT[ib];int sum6=A.sum+B.sum;int rem=parentM-sum6;int lo=max(0,rem-16),hi=min(16,rem);if(lo>hi)continue;double bv[7],n0=0,bd=0;for(int k=0;k<7;k++){bv[k]=A.v[k]+B.v[k]+rem*V[7][k];n0+=bv[k]*bv[k];bd+=bv[k]*Dv[k];}bool even6=A.even&&B.even; long long halfbase=((long long)A.half9*729+B.half9)*9; for(int j=lo;j<=hi;j++){int last=rem-j;long long p=rr*PAD+j;double n=n0+2*j*bd+j*(double)j*DD;double t0=b0[p]*inv,t1=b1[p]*inv,t2=w[(done+rr)*PAD+j]*inv;if(even6 && ((j|last)&1)==0){long long hidx=halfbase+j/2;t0+=d0[hidx]*ds;t1+=d1[hidx]*ds;t2+=d2[hidx]*ds;}double an=ap*n,fac=.5*exp(alpha*an);z0[p]=fac*t0;z1[p]=fac*(t1+an*t0);z2[p]=fac*(t2+2*an*t1+an*an*t0);} }
  fwrite(z0.data(),sizeof(double),cnt,o0);fwrite(z1.data(),sizeof(double),cnt,o1);fwrite(z2.data(),sizeof(double),cnt,o2);done+=nr;if(done%(BR*128)==0)cerr<<" buildM"<<parentM<<" rows "<<done<<"/"<<N6<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"\n";
 }
 fclose(f0);fclose(f1);fclose(o0);fclose(o1);fclose(o2);cerr<<"build parent M"<<parentM<<" complete\n";
}
int main(int argc,char**argv){if(argc<6){cerr<<"usage base alpha q2prefix workprefix q2meta\n";return 2;}string base=argv[1],q2p=argv[3],wp=argv[4],meta=argv[5];double alpha=stod(argv[2]);prep_triples();
 vector<BaseRow> rows; rows.reserve(250000); ifstream in(base);string line;double M=-1e300;while(getline(in,line)){istringstream ss(line);string hs,a,b,c;ss>>hs>>a>>b>>c;if(!ss)continue;auto cc=dec5(hs);if(a=="nan")continue;BaseRow r{cc,ridx7(cc),stod(a),stod(b),stod(c),normc(cc)};double lg=r.lz-alpha*(4.0/256.0)*r.norm;M=max(M,lg);rows.push_back(r);}cerr<<setprecision(17)<<"rows "<<rows.size()<<" S16 "<<M<<"\n";
 double*w=fftw_alloc_real(NREAL);if(!w){cerr<<"alloc fail\n";return 3;}int dims[7]={17,17,17,17,17,17,17};fftw_init_threads();fftw_plan_with_nthreads(5);fftw_plan pf=fftw_plan_dft_r2c(7,dims,w,(fftw_complex*)w,FFTW_ESTIMATE),pb=fftw_plan_dft_c2r(7,dims,(fftw_complex*)w,w,FFTW_ESTIMATE);
 auto fillbase=[&](int fld){memset(w,0,NREAL*sizeof(double));double a16=4.0/256.0;for(auto&r:rows){double lg=r.lz-alpha*a16*r.norm,z=exp(lg-M),l1=r.d1-a16*r.norm,val=fld==0?z:(fld==1?z*l1:z*(r.d2+l1*l1));w[r.ridx]=val;}};
 auto level=[&](int lev,int parentM,double Schild,const string&diagpref,const string&childpref,const string&outpref){string tag=wp+"L"+to_string(lev)+"_",F0=tag+"F0.bin",F1=tag+"F1.bin",C0=tag+"C0.bin",C1=tag+"C1.bin";auto loadfld=[&](int fld){if(lev==1)fillbase(fld);else read_full(childpref+"g"+to_string(fld)+".raw",w,NREAL);};auto tm=[&](){return chrono::steady_clock::now();};
  loadfld(0);auto t=tm();fftw_execute(pf);cerr<<"L"<<lev<<" F0 fft "<<chrono::duration<double>(tm()-t).count()<<"\n";write_full(F0,w,NREAL);{fftw_complex*W=(fftw_complex*)w;for(long long i=0;i<NCPLX;i++){double ar=W[i][0],ai=W[i][1];W[i][0]=ar*ar-ai*ai;W[i][1]=2*ar*ai;}}t=tm();fftw_execute(pb);cerr<<"L"<<lev<<" C0 ifft "<<chrono::duration<double>(tm()-t).count()<<"\n";write_full(C0,w,NREAL);
  loadfld(1);t=tm();fftw_execute(pf);cerr<<"L"<<lev<<" F1 fft "<<chrono::duration<double>(tm()-t).count()<<"\n";write_full(F1,w,NREAL);complex_combine(w,F0,"",1);t=tm();fftw_execute(pb);cerr<<"L"<<lev<<" C1 ifft "<<chrono::duration<double>(tm()-t).count()<<"\n";write_full(C1,w,NREAL);
  loadfld(2);t=tm();fftw_execute(pf);cerr<<"L"<<lev<<" F2 fft "<<chrono::duration<double>(tm()-t).count()<<"\n";complex_combine(w,F0,F1,2);t=tm();fftw_execute(pb);cerr<<"L"<<lev<<" C2 ifft "<<chrono::duration<double>(tm()-t).count()<<"\n";if(lev==2&&!childpref.empty()){remove((childpref+"g0.raw").c_str());remove((childpref+"g1.raw").c_str());remove((childpref+"g2.raw").c_str());}
  build_parent(C0,C1,w,diagpref,alpha,parentM,Schild,outpref);remove(F0.c_str());remove(F1.c_str());remove(C0.c_str());remove(C1.c_str()); };
 string p32=wp+"M32_",p64=wp+"M64_";level(1,32,M,q2p+"q2m16_","",p32);double S32=2*M;level(2,64,S32,q2p+"q2m32_",p32,p64);double S64=2*S32;
 fftw_destroy_plan(pf);fftw_destroy_plan(pb);fftw_free(w);fftw_cleanup_threads();
 // parse q2m64 homogeneous from meta (very simple: find array after key)
 ifstream mf(meta);string ms((istreambuf_iterator<char>(mf)),{});auto pos=ms.find("\"q2m64_hom\"");pos=ms.find('[',pos);auto end=ms.find(']',pos);string arr=ms.substr(pos+1,end-pos-1);replace(arr.begin(),arr.end(),',',' ');istringstream as(arr);double Q0,Q1,Q2;as>>Q0>>Q1>>Q2;double eqscale=exp(-2*S64);Q0*=eqscale;Q1*=eqscale;Q2*=eqscale;
 // mmap final three files
 vector<int>fd(3);vector<double*>g(3);for(int k=0;k<3;k++){string p=p64+"g"+to_string(k)+".raw";fd[k]=open(p.c_str(),O_RDONLY);if(fd[k]<0){perror("open g64");return 4;}g[k]=(double*)mmap(nullptr,NREAL*sizeof(double),PROT_READ,MAP_SHARED,fd[k],0);if(g[k]==MAP_FAILED){perror("mmap");return 4;}}
 long double o0=0,o1=0,o2=0,ed=0,ed2=0,eint=0; long long half=(NGRID-1)/2;auto st=chrono::steady_clock::now();for(long long q=0;q<=half;q++){long long qc=NGRID-1-q;long long r=q/17,j=q%17,rc=qc/17,jc=qc%17;long long i=r*PAD+j,ic=rc*PAD+jc;double a0=g[0][i],b0=g[0][ic];if(a0==0||b0==0)continue;double a1=g[1][i],b1=g[1][ic],a2=g[2][i],b2=g[2][ic];long double mult=(q==qc)?1.0L:2.0L;long double ww=mult*a0*b0;double da=a1/a0,db=b1/b0,d=da+db;double la2=a2/a0-da*da,lb2=b2/b0-db*db;o0+=ww;o1+=ww*d;o2+=ww*(d*d+la2+lb2);ed+=ww*d;ed2+=ww*d*d;eint+=ww*(la2+lb2);if(q>0&&q%50000000LL==0)cerr<<"root q "<<q<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"\n";}
 // equal branch q2
 double deq=Q1/Q0,leq=Q2/Q0-deq*deq;o0+=Q0;o1+=Q1;o2+=Q2;ed+=Q0*deq;ed2+=Q0*deq*deq;eint+=Q0*leq;
 long double mean=o1/o0,total=o2/o0-mean*mean,root=ed2/o0-(ed/o0)*(ed/o0),internal=eint/o0;long double CH=alpha*alpha*total;cerr<<setprecision(17);cout<<setprecision(17)<<"alpha "<<alpha<<"\nS16 "<<M<<" S64 "<<S64<<"\npeq "<<(double)(Q0/o0)<<"\nCH "<<(double)CH<<"\nroot_CH "<<(double)(alpha*alpha*root)<<"\ninternal_CH "<<(double)(alpha*alpha*internal)<<"\nroot_fraction "<<(double)(root/total)<<"\nroot_deficit "<<(double)(1-root/total)<<"\nlogder2 "<<(double)total<<"\n";
 for(int k=0;k<3;k++){munmap(g[k],NREAL*sizeof(double));close(fd[k]);}
}
