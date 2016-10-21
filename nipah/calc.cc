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
                       vector<double>& distOutbSizes, int N, int M, double R,
                       vector<double>& distSumOutbSizesPrimary)
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
    getDistSumPrimary(alpha, nu, outbsPrimaryInfo, distOutbSizes, N, M, R, distSumOutbSizesPrimary);

    double prob = distSumOutbSizesPrimary[sumPrimaryOutbs];
    
    for(int i : outbsNoPrimaryInfo)
        prob *= distOutbSizes[i];
    
    for(int j : outbsPrimaryInfo)
        prob *= distOutbSizes[j];
    
    return prob;
};

void getPmPrimaryInfo(double alpha, double nu, vector<int>& outbsPrimaryInfo, int sumPrimaryOutbs, int targetOutbSize, int N, int M, double R, vector<double>& Pmn)
{
    int lenPmn = targetOutbSize + 1;
    vector<double> distOutbSizes, distSumPrimary;
    getDistOutbSizes(alpha, nu, distOutbSizes, N, R);
    getDistSumPrimary(alpha, nu, outbsPrimaryInfo, distOutbSizes, N, M, R, distSumPrimary);
    
    cdouble **fftInputL1, **fftOutputL1, *fftInputL2, *fftOutputL2;
    fftInputL1 = new cdouble*[lenPmn];
    fftOutputL1 = new cdouble*[lenPmn];
    fftInputL2 = new cdouble[M];
    fftOutputL2 = new cdouble[M];
    
    for( int k = 0; k < lenPmn; k++)
    {
        fftInputL1[k] = new cdouble[M];
        fftOutputL1[k] = new cdouble[M];
        for( int i = 0; i < M; i++ )
        {
            fftInputL2 = new cdouble[N];
            fftOutputL2 = new cdouble[N];
        
            getfftInput1D_Gxy( cdouble(R,0) * exp(cdouble(0,2*PI*i/M)),
                               cdouble(alpha,0), cdouble(nu,0), R, N, fftInputL2);
            runfft1D( N, fftInputL2, fftOutputL2 );
        
            fftInputL1[k][i] = cdouble(1,0);
            for(int s : outbsPrimaryInfo)
            {
                if( s != targetOutbSize )
                fftInputL1[k][i] *= fftOutputL2[s] / distOutbSizes[s];
            };
            
            if( targetOutbSize == 12) //12 is repeated so adding it back
                fftInputL1[k][i] *= fftOutputL2[12] / distOutbSizes[12];
            
            delete fftInputL2;
            delete fftOutputL2;
            
            fftInputL2 = new cdouble[N];
            fftOutputL2 = new cdouble[N];

            getfftInput1D_Gxy( cdouble(R,0) * exp(cdouble(0,2*PI*i/M)) *
                              cdouble(R,0) * exp(cdouble(0,2*PI*k/(lenPmn))),
                              cdouble(alpha,0), cdouble(nu,0), R, N, fftInputL2);
            runfft1D( N, fftInputL2, fftOutputL2 );
            
            fftInputL1[k][i] *= fftOutputL2[targetOutbSize] / distOutbSizes[targetOutbSize];
        };
    };
    runfft2D( lenPmn, M, fftInputL1, fftOutputL1);
    
    for(int i = 0; i < lenPmn; i++)
        Pmn.push_back(real(fftOutputL1[i][sumPrimaryOutbs]) / distSumPrimary[sumPrimaryOutbs]);
    
    delete fftInputL1;
    delete fftOutputL1;
};

void getPmNoPrimaryInfo(double alpha, double nu, int N, double R,
                        vector< vector<double> >& Pmn)
{
    vector<double> distOutbSizes, temprow;
    getDistOutbSizes(alpha, nu, distOutbSizes, N, R);
    
    cdouble **fftInput, **fftOutput;
    fftInput = new cdouble*[N];
    fftOutput = new cdouble*[N];
    for( int i = 0; i < N; i++ )
    {
        fftInput[i] = new cdouble[N];
        fftOutput[i] = new cdouble[N];
        for( int j = 0; j < N; j++ )
        {
            fftInput[i][j] = Gxy(cdouble(R,0) * exp(cdouble(0,2*PI*j/N)), cdouble(R,0) * exp(cdouble(0,2*PI*i/N)), alpha, nu);
        };
    };
    runfft2D( N, N, fftInput, fftOutput);

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            temprow.push_back(real(fftOutput[i][j])/distOutbSizes[i]);
        };
        Pmn.push_back(temprow);
        temprow.clear();
    };
};


