#include <bits/stdc++.h>
extern "C" {
typedef double fftw_complex[2]; typedef struct fftw_plan_s *fftw_plan;
double *fftw_alloc_real(size_t); void fftw_free(void*);
fftw_plan fftw_plan_dft_r2c(int,const int*,double*,fftw_complex*,unsigned);
void fftw_execute(const fftw_plan); void fftw_destroy_plan(fftw_plan);
int fftw_init_threads(void); void fftw_plan_with_nthreads(int); void fftw_cleanup_threads(void);
}
#define FFTW_ESTIMATE (1U<<6)
int main(){const long long N6=24137569LL,NREAL=N6*18;std::cerr<<"bytes "<<NREAL*8.0/1e9<<"\n";double*a=fftw_alloc_real(NREAL);if(!a){std::cerr<<"alloc fail\n";return 2;}memset(a,0,NREAL*sizeof(double));a[0]=1;int d[7]={17,17,17,17,17,17,17};fftw_init_threads();fftw_plan_with_nthreads(5);auto p=fftw_plan_dft_r2c(7,d,a,(fftw_complex*)a,FFTW_ESTIMATE);auto t=std::chrono::steady_clock::now();fftw_execute(p);std::cerr<<"fft "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<"\n";fftw_destroy_plan(p);fftw_free(a);fftw_cleanup_threads();}
