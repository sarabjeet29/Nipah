#include "fft.h"

void getfftInput1D_Gxy( cdouble x, cdouble alpha, cdouble nu, double R, int N,
                       vector<cdouble>& fftInput )
{
	for( int i = 0; i < N; i++ )
		fftInput.push_back(Gxy(x, cdouble(R,0) * exp(cdouble(0,2*PI*i/N)), alpha, nu));
};

void getfftInput2D_Hxzuw( cdouble x, cdouble u, cdouble alpha, cdouble nu, cdouble p, cdouble q, double R, int N, vector<vector<cdouble>>& fftInput)
{
    vector<cdouble> row;
    for( int i = 0; i < N; i++ )
    {
        for(int j = 0; j < N; j++)
            row.push_back(Hxzuw(x, cdouble(R,0) * exp(cdouble(0,2*PI*i/N)), u, cdouble(R,0) * exp(cdouble(0,2*PI*j/N)), alpha, nu, p, q));

        fftInput.push_back(row);
        row.clear();
    };
};


void runfft1D( int N, vector<cdouble>& fftInput, vector<cdouble>& fftOutput )
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
		fftOutput.push_back(cdouble((REAL(out,j))/N, IMAG(out,j)/N));

	fftw_destroy_plan(plan);
	fftw_free(in); fftw_free(out);
};

void getfftInput2D(int M, int N, vector<vector<cdouble>>& rawInput,
                   vector<cdouble>& fftInput)
{
    for( int i = 0; i < M; i++)
    {
        for(int j = 0; j < N; j++)
            fftInput.push_back(rawInput[i][j]);
    };
};

void runfft2D( int M, int N, vector<vector<cdouble>>& rawInput,
               vector<vector<cdouble>>& fftOutput)
 {
     fftw_complex *in, *out;
     fftw_plan plan;
     int size = M * N, index1D;
     
     vector<cdouble> fftInput, row;
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
             row.push_back(cdouble(REAL(out,index1D) / size, IMAG(out,index1D) / size));
         };
         fftOutput.push_back(row);
         row.clear();
     };
     
 	 fftw_destroy_plan(plan);
 	 fftw_free(in); fftw_free(out);
 };

void runfft2DflatInput(int M, int N, vector<cdouble>& fftInput, vector<vector<double>>& fftOutput)
{
    vector<double> row;
    fftw_complex *in, *out;
    fftw_plan plan;
    int size = M * N, index1D;
    
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
            row.push_back(REAL(out,index1D) / size);
        };
        fftOutput.push_back(row);
        row.clear();
    };
    
    fftw_destroy_plan(plan);
    fftw_free(in); fftw_free(out);
};

void runfft2DflatIO(int M, int N, vector<cdouble>& fftInput, vector<cdouble>& fftOutput)
{
    fftw_complex *in, *out;
    fftw_plan plan;
    int size = M * N;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
    
    for( int i = 0; i < size; i++ )
    {
        REAL(in,i) = real( fftInput[i] );
        IMAG(in,i) = imag( fftInput[i] );
    };
    
    plan = fftw_plan_dft_2d(M, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan); /* repeat as needed */
    
    for( int i = 0; i < size; i++)
        fftOutput.push_back(cdouble(REAL(out,i) / size, IMAG(out,i) / size));
    
    fftw_destroy_plan(plan);
    fftw_free(in); fftw_free(out);
};

