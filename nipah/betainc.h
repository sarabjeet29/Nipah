#include <cstdlib>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_gamma.h>
#include <complex>
#include <ctime>
#include <vector>

#define MAXIT 100000000
#define EPS 1e-15 
#define FPMIN 1.0e-30
using namespace std;

typedef complex<double> cdouble;

#define REAL(z,i) ((z)[i][0])
#define IMAG(z,i) ((z)[i][1])

cdouble lngamma( cdouble );
cdouble gamma( cdouble );
cdouble betai( cdouble, cdouble, cdouble );
cdouble upperbetai( cdouble, cdouble, cdouble );
cdouble betacf( cdouble, cdouble, cdouble );
void nrerror(const char[]);
