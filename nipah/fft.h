#include <fftw3.h>
#include "pgf.h"

void getfftInput1D_Gxy( cdouble, cdouble, cdouble, double, int, cdouble* );
void runfft1D( int, cdouble*, cdouble* );
void runfft2D( int, int, cdouble **, cdouble **);
