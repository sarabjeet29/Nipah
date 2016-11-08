#include "fft.h"

void getfftInput1D_Gxy(cdouble x, cdouble alpha, cdouble nu, double R, int N,
                       fftw_complex *fftInput )
{
    cdouble tempval;
    for( int i = 0; i < N; i++ )
    {
        tempval = Gxy(x, cdouble(R,0) * exp(cdouble(0,2*PI*i/N)), alpha, nu) / cdouble(N,0);
        REAL(fftInput, i) = real(tempval);
        IMAG(fftInput, i) = imag(tempval);
    }
};

void getfftInput2D_Hxzuw( cdouble x, cdouble u, cdouble alpha, cdouble nu, cdouble p, cdouble q, double R, int N, fftw_complex *fftInput)
{
    int index;
    cdouble tempval;
    for( int i = 0; i < N; i++ )
    {
        for(int j = 0; j < N; j++)
        {
            index = i * N + j;
            tempval = Hxzuw(x, cdouble(R,0) * exp(cdouble(0,2*PI*i/N)), u, cdouble(R,0) * exp(cdouble(0,2*PI*j/N)), alpha, nu, p, q) / cdouble(N*N,0);
            REAL(fftInput, index) = real(tempval);
            IMAG(fftInput, index) = imag(tempval);
        };
    };
};

void runfft1D( int N, fftw_complex *fftInput, fftw_complex *fftOutput )
{
	fftw_plan plan;
	
//	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
	plan = fftw_plan_dft_1d(N, fftInput, fftOutput, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
    
//	for( int j = 0; j < N; j++ )
//		fftOutput.push_back(cdouble((REAL(out,j))/N, IMAG(out,j)/N));

	fftw_destroy_plan(plan);
//	fftw_free(in); fftw_free(out);
};

void runfft2D( int M, int N, fftw_complex *fftInput, fftw_complex *fftOutput)
 {
     fftw_plan plan;
//     int size = M * N, index1D;
     
//     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
//     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);

     plan = fftw_plan_dft_2d(M, N, fftInput, fftOutput, FFTW_FORWARD, FFTW_ESTIMATE);
     fftw_execute(plan); /* repeat as needed */

//     for( int i = 0; i < M; i++)
//     {
//         for( int j = 0; j < N; j++)
//         {
//             index1D = i * N + j;
//             row.push_back(cdouble(REAL(out,index1D) / size, IMAG(out,index1D) / size));
//         };
//         fftOutput.push_back(row);
//         row.clear();
//     };
     
 	 fftw_destroy_plan(plan);
// 	 fftw_free(in); fftw_free(out);
 };
