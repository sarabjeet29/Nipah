#include <fftw3.h>
#include "pgf.h"

void getfftInput1D_Gxy( cdouble, cdouble, cdouble, double, int, vector<cdouble>& );
void getfftInput2D_Hxzuw( cdouble, cdouble, cdouble, cdouble, cdouble, cdouble, double, int, vector<vector<cdouble>>&);
void runfft1D( int, vector<cdouble>&, vector<cdouble>& );
void runfft2D( int, int, vector<vector<cdouble>>&, vector<vector<cdouble>>&);
void runfft2DflatInput(int M, int N, vector<cdouble>& fftInput, vector<vector<double>>& fftOutput);
void runfft2DflatIO(int M, int N, vector<cdouble>& fftInput, vector<cdouble>& fftOutput);