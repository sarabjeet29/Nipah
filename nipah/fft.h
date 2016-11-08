#include <fftw3.h>
#include "pgf.h"

void getfftInput1D_Gxy(cdouble, cdouble, cdouble, double, int, fftw_complex* );
void getfftInput2D_Hxzuw(cdouble, cdouble, cdouble, cdouble, cdouble, cdouble, double,
                         int, fftw_complex*);
void runfft1D( int, fftw_complex*, fftw_complex* );
void runfft2D( int, int, fftw_complex*, fftw_complex*);
