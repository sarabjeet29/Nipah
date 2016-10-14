#include "calc.h"

void getDistOutbSizes( double alpha, double nu, vector<double>& distOutbSizes, int N, double R )
{
	cdouble *fftInput, *fftOutput;
	fftInput = new cdouble[N];
    fftOutput = new cdouble[N];
	getfftInput1D_Gxy( cdouble(1,0), cdouble(alpha,0), cdouble(nu,0), R, N, fftInput );
	runfft1D( N, fftInput, fftOutput );
    
    for(int i = 0; i < N; i++)
        distOutbSizes.push_back(real(fftOutput[i]));
    
    delete fftInput;
    delete fftOutput;
};

void getDistSumPrimary(double alpha, double nu, vector<int>& outbsPrimaryInfo,
                       vector<double>& distOutbSizes, vector<double>& distSumOutbSizesPrimary, int N, int M, double R)
{
    cdouble *outerfftInput, *outerfftOutput, *innerfftInput, *innerfftOutput;
    outerfftInput = new cdouble[M];
    outerfftOutput = new cdouble[M];
   
    for( int i = 0; i < M; i++ )
    {
        innerfftInput = new cdouble[N];
        innerfftOutput = new cdouble[N];

        getfftInput1D_Gxy( cdouble(R,0) * exp(cdouble(0,2*PI*i/M)), cdouble(alpha,0), cdouble(nu,0), R, N, innerfftInput);
        runfft1D( N, innerfftInput, innerfftOutput );
        
        outerfftInput[i] = cdouble(1,0);
        for(int j : outbsPrimaryInfo)
            outerfftInput[i] *= innerfftOutput[j] / cdouble(distOutbSizes[j],0);

        delete innerfftInput;
        delete innerfftOutput;
    };
    runfft1D( M, outerfftInput, outerfftOutput);

    for(int i = 0; i < M; i++)
        distSumOutbSizesPrimary.push_back(real(outerfftOutput[i]));
    
    delete outerfftInput;
    delete outerfftOutput;
};

double getProbObs(double alpha, double nu, vector<int>& outbsPrimaryInfo,
                  vector<int>& outbsNoPrimaryInfo, int sumPrimaryOutbs, int N, int M,
                  double R )
{
    vector<double> distOutbSizes, distSumOutbSizesPrimary;
    getDistOutbSizes(alpha, nu, distOutbSizes, N, R);
    getDistSumPrimary(alpha, nu, outbsPrimaryInfo, distOutbSizes, distSumOutbSizesPrimary, N, M, R);

    double prob = distSumOutbSizesPrimary[sumPrimaryOutbs];
    
    for(int i : outbsNoPrimaryInfo)
        prob *= distOutbSizes[i];
    
    for(int j : outbsPrimaryInfo)
        prob *= distOutbSizes[j];
    
    return prob;
};




