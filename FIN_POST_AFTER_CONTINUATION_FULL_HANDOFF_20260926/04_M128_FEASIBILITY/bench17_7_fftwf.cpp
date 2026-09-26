#include <bits/stdc++.h>
extern "C" {
typedef float fftwf_complex[2]; typedef struct fftwf_plan_s *fftwf_plan;
float *fftwf_alloc_real(size_t); void fftwf_free(void*);
fftwf_plan fftwf_plan_dft_r2c(int,const int*,float*,fftwf_complex*,unsigned);
void fftwf_execute(const fftwf_plan); void fftwf_destroy_plan(fftwf_plan);
int fftwf_init_threads(void); void fftwf_plan_with_nthreads(int); void fftwf_cleanup_threads(void);
}
#define FFTW_ESTIMATE (1U<<6)
int main(){const long long N6=24137569LL; const long long NREAL=N6*18; std::cerr<<"NREAL "<<NREAL<<" bytes "<<NREAL*4.0/1e9<<"\n"; float* a=fftwf_alloc_real(NREAL); if(!a){std::cerr<<"alloc fail\n";return 2;} memset(a,0,NREAL*sizeof(float)); a[0]=1; int dims[7]={17,17,17,17,17,17,17}; fftwf_init_threads(); fftwf_plan_with_nthreads(5); auto t=std::chrono::steady_clock::now(); auto p=fftwf_plan_dft_r2c(7,dims,a,(fftwf_complex*)a,FFTW_ESTIMATE); std::cerr<<"plan "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<"\n"; t=std::chrono::steady_clock::now(); fftwf_execute(p); std::cerr<<"fft "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<" first "<<a[0]<<"\n"; fftwf_destroy_plan(p);fftwf_free(a);fftwf_cleanup_threads();}
