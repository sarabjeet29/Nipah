#include "fft.h"

void getfftInput1D_Gxy( cdouble x, cdouble alpha, cdouble nu, double R, int N, cdouble fftInput[] )
{
	for( int i = 0; i < N; i++ )
	{
		fftInput[i] = Gxy(x, cdouble(R,0) * exp(cdouble(0,2*PI*i/N)), alpha, nu);
	};
};

void runfft1D( int N, cdouble fftInput[], cdouble fftOutput[] )
{
	fftw_complex *in, *out;
	fftw_plan plan;
	
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
	for( int i = 0; i < N; i++ )
	{
		REAL(in,i) = real( fftInput[i] );
		IMAG(in,i) = imag( fftInput[i] );
	};		
	
	plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
    
	for( int j = 0; j < N; j++ )
		fftOutput[j] = cdouble( (REAL(out,j))/N, IMAG(out,j)/N );

	fftw_destroy_plan(plan);
	fftw_free(in); fftw_free(out);
};

void getfftInput2D(int M, int N, cdouble * rawInput[], cdouble fftInput[])
{
    for( int i = 0; i < M; i++)
    {
        for(int j = 0; j < N; j++)
        {
            fftInput[i*N+j] = rawInput[i][j];
        };
    };
};

void runfft2D( int M, int N, cdouble * rawInput[], cdouble * fftOutput[] )
 {
     fftw_complex *in, *out;
     fftw_plan plan;
     int size = M * N, index1D;
     
     cdouble * fftInput;
     fftInput = new cdouble[size];
     getfftInput2D(M, N, rawInput, fftInput);

     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);

     for( int i = 0; i < size; i++ )
     {
         REAL(in,i) = real( fftInput[i] );
         IMAG(in,i) = imag( fftInput[i] );
     };

     plan = fftw_plan_dft_2d(M, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
     fftw_execute(plan); /* repeat as needed */

     for( int i = 0; i < M; i++)
     {
         for( int j = 0; j < N; j++)
         {
             index1D = i * N + j;
             fftOutput[i][j] = cdouble(REAL(out,index1D) / size, IMAG(out,index1D) / size);
         };
     };
     
 	 fftw_destroy_plan(plan);
 	 fftw_free(in); fftw_free(out);
 };
