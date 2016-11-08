#include "calc.h"
#include <iostream>

void getDistOutbs(double alpha, double nu, vector<double>& distOutbs, int N, double R )
{
	fftw_complex *fftInput, *fftOutput;
    fftInput  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	getfftInput1D_Gxy( cdouble(1,0), cdouble(alpha,0), cdouble(nu,0), R, N, fftInput );
	runfft1D( N, fftInput, fftOutput );
    
    for(int i = 0; i < N; i++)
        distOutbs.push_back(REAL(fftOutput, i));
    
    fftw_free(fftInput); fftw_free(fftOutput);
};

void getDistSumPrimaryOutbs(double alpha, double nu, vector<int>& outbsPrimaryInfo,
                            vector<double>& distOutbs, int N, int M, double R,
                            vector<double>& distSumOutbSizesPrimary)
{
    cdouble tempval;
    fftw_complex *outerfftInput, *outerfftOutput, *innerfftInput, *innerfftOutput;

    outerfftInput  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    outerfftOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    
    for( int i = 0; i < M; i++ )
    {
        innerfftInput  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        innerfftOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

        getfftInput1D_Gxy( cdouble(R,0) * exp(cdouble(0,2*PI*i/M)), cdouble(alpha,0), cdouble(nu,0), R, N, innerfftInput);
        runfft1D( N, innerfftInput, innerfftOutput );
        
        tempval = cdouble(1,0);
        for(int j : outbsPrimaryInfo)
            tempval *= cdouble(REAL(innerfftOutput, j), IMAG(innerfftOutput, j)) / cdouble(distOutbs[j],0);

        REAL(outerfftInput,i) = real(tempval);
        IMAG(outerfftInput,i) = imag(tempval);
        
        fftw_free(innerfftInput); fftw_free(innerfftOutput);
    };
    runfft1D( M, outerfftInput, outerfftOutput);

    for(int i = 0; i < M; i++)
        distSumOutbSizesPrimary.push_back(REAL(outerfftOutput, i));
    
    fftw_free(outerfftInput); fftw_free(outerfftOutput);
};

double getProbObsOutbs(double alpha, double nu, vector<int>& outbsPrimaryInfo,
                       vector<int>& outbsNoPrimaryInfo, int sumPrimaryOutbs, int N,
                       int M, double R )
{
    vector<double> distOutbs, distSumOutbSizesPrimary;
    getDistOutbs(alpha, nu, distOutbs, N, R);
    getDistSumPrimaryOutbs(alpha, nu, outbsPrimaryInfo, distOutbs, N, M, R,
                           distSumOutbSizesPrimary);

    double prob = distSumOutbSizesPrimary[sumPrimaryOutbs];
    
    for(int i : outbsNoPrimaryInfo)
        prob *= distOutbs[i];
    
    for(int j : outbsPrimaryInfo)
        prob *= distOutbs[j];
    
    return prob;
};

void getPmPrimaryInfo(double alpha, double nu, vector<int>& outbsPrimaryInfo, int sumPrimaryOutbs, int targetOutbSize, int N, int M, double R, vector<double>& Pmn)
{
    int lenPmn = targetOutbSize + 1;
    vector<double> distOutbs, distSumPrimary;
    getDistOutbs(alpha, nu, distOutbs, N, R);
    getDistSumPrimaryOutbs(alpha, nu, outbsPrimaryInfo, distOutbs, N, M, R, distSumPrimary);
    
    fftw_complex *outerfftInput, *outerfftOutput, *innerfftInput, *innerfftOutput;
    cdouble ctempVal;
    
    outerfftInput  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * lenPmn);
    outerfftOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * lenPmn);
    
    for( int i = 0; i < M; i++ )
    {
        innerfftInput  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        innerfftOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        
        getfftInput1D_Gxy(cdouble(R,0) * exp(cdouble(0,2*PI*i/M)),
                          cdouble(alpha,0), cdouble(nu,0), R, N, innerfftInput);
        runfft1D( N, innerfftInput, innerfftOutput );
        
        ctempVal = cdouble(1,0);
        for(int s : outbsPrimaryInfo)
        {
            if( s != targetOutbSize )
                ctempVal *= cdouble(REAL(innerfftOutput,s), IMAG(innerfftOutput,s)) / distOutbs[s];
        };
            
        if(targetOutbSize == 12) //12 is repeated so adding it back
            ctempVal *= cdouble(REAL(innerfftOutput,12), IMAG(innerfftOutput,12)) / distOutbs[12];
            
        fftw_free(innerfftInput); fftw_free(innerfftOutput);

        for( int k = 0; k < lenPmn; k++)
        {
            innerfftInput  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
            innerfftOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
            getfftInput1D_Gxy(cdouble(R,0) * exp(cdouble(0,2*PI*i/M)) *
                              cdouble(R,0) * exp(cdouble(0,2*PI*k/(lenPmn))),
                              cdouble(alpha,0), cdouble(nu,0), R, N, innerfftInput);

            runfft1D(N, innerfftInput, innerfftOutput);
            
            ctempVal *= cdouble(REAL(innerfftOutput, targetOutbSize),
                                IMAG(innerfftOutput, targetOutbSize)) / distOutbs[targetOutbSize];
            
            REAL(outerfftInput, i * lenPmn + k) = real(ctempVal);
            IMAG(outerfftInput, i * lenPmn + k) = imag(ctempVal);
            
            fftw_free(innerfftInput); fftw_free(innerfftOutput);
        };
    };
    runfft2D(M, lenPmn, outerfftInput, outerfftOutput);
    
    for(int j = 0; j < lenPmn; j++)
        Pmn.push_back(REAL(outerfftOutput, sumPrimaryOutbs * lenPmn + j) /
                      distSumPrimary[sumPrimaryOutbs]);

    fftw_free(outerfftInput); fftw_free(outerfftOutput);
};

void getPmNoPrimaryInfo(double alpha, double nu, int N, double R, vector<vector<double>>& Pmn)
{
    vector<double> distOutbs, temprow;
    getDistOutbs(alpha, nu, distOutbs, N, R);
    
    cdouble ctempval;
    fftw_complex *fftInput, *fftOutput;
    fftInput  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    fftOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);

    for( int i = 0; i < N; i++ )
    {
        for( int j = 0; j < N; j++ )
        {
            ctempval = Gxy(cdouble(R,0) * exp(cdouble(0,2*PI*j/N)), cdouble(R,0) *
                           exp(cdouble(0,2*PI*i/N)), alpha, nu);
            REAL(fftInput, i * N + j) = real(ctempval);
            IMAG(fftInput, i * N + j) = imag(ctempval);
        };
    };
    runfft2D( N, N, fftInput, fftOutput);

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            temprow.push_back(REAL(fftOutput, i * N + j) / distOutbs[i]);

        Pmn.push_back(temprow);
        temprow.clear();
    };

    fftw_free(fftInput); fftw_free(fftOutput);
};

void getDistOutbsAndDeaths(double alpha, double nu, double p, double q,
                           vector<vector<double>>& distOutbsAndDeaths, int N, double R )
{
    vector<double> temprow;
    fftw_complex *fftInput, *fftOutput;
    fftInput  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    fftOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);

    getfftInput2D_Hxzuw( cdouble(1,0), cdouble(1,0), cdouble(alpha,0), cdouble(nu,0), cdouble(p,0), cdouble(q,0), R, N, fftInput);
    
    runfft2D( N, N, fftInput, fftOutput);
    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            temprow.push_back(REAL(fftOutput, i * N + j));

        distOutbsAndDeaths.push_back(temprow);
        temprow.clear();
    };

    fftw_free(fftInput); fftw_free(fftOutput);
};

void getDistSumPrimaryOutbsAndDeaths(double alpha, double nu, double p, double q,
                             vector<int>& outbsPrimaryInfo,
                             vector<int>& deathsPrimaryInfo,
                             vector<vector<double>>& distOutbsAndDeaths,
                             int N, int M, double R,
                             vector<vector<double>>& distSumPrimaryOutbsAndDeaths)
{
    int outbSize, nDeaths;
    cdouble ctempval;
    vector<double> tempRow;
    vector<cdouble> ctempRow;
    fftw_complex *outerfftInput, *outerfftOutput, *innerfftInput, *innerfftOutput;
    
    outerfftInput  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * M);
    outerfftOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * M);

    for( int i = 0; i < M; i++ )
    {
        cout << i <<" ";
        for( int j = 0; j < M; j++ )
        {
            innerfftInput  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
            innerfftOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);

            getfftInput2D_Hxzuw( cdouble(R,0) * exp(cdouble(0,2*PI*i/M)), cdouble(R,0) * exp(cdouble(0,2*PI*j/M)), alpha, nu, p, q, R, N, innerfftInput );

            runfft2D( N, N, innerfftInput, innerfftOutput );
            
            ctempval = cdouble(1,0);
            for(int k = 0; k < outbsPrimaryInfo.size(); k++)
            {
                outbSize = outbsPrimaryInfo[k];
                nDeaths = deathsPrimaryInfo[k];
                ctempval *= cdouble(REAL(innerfftOutput, outbSize * N + nDeaths),
                                    IMAG(innerfftOutput, outbSize * N + nDeaths)) /distOutbsAndDeaths[outbSize][nDeaths];
            };
            fftw_free(innerfftInput); fftw_free(innerfftOutput);
            REAL(outerfftInput, i * M + j) = real(ctempval);
            REAL(outerfftInput, i * M + j) = imag(ctempval);
        };
    };
    cout << "\n";
    runfft2D( M, M, outerfftInput, outerfftOutput);
    
    for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < M; j++)
            tempRow.push_back(REAL(outerfftOutput, i * M + j));
        
        distSumPrimaryOutbsAndDeaths.push_back(tempRow);
        tempRow.clear();
    };
};

double getProbObsOubsAndDeaths(double alpha, double nu, double p, double q,
                               vector<int>& outbsPrimaryInfo,
                               vector<int>& outbsNoPrimaryInfo,
                               vector<int>& deathsPrimaryInfo,
                               vector<int>& deathsNoPrimaryInfo, int sumPrimaryOutbs,
                               int sumPrimaryDeaths, int N, int M, double R )
{
    vector<vector<double>> distOutbsAndDeaths, distSumPrimaryOutbsAndDeaths;
    getDistOutbsAndDeaths(alpha, nu, p, q,distOutbsAndDeaths, N, R );
    getDistSumPrimaryOutbsAndDeaths(alpha, nu, p, q, outbsPrimaryInfo, deathsPrimaryInfo, distOutbsAndDeaths, N, M, R, distSumPrimaryOutbsAndDeaths);
    
    double prob = distSumPrimaryOutbsAndDeaths[sumPrimaryOutbs][sumPrimaryDeaths];
    
    for( int i = 0; i < outbsNoPrimaryInfo.size(); i++ )
        prob *= distOutbsAndDeaths[outbsNoPrimaryInfo[i]][deathsNoPrimaryInfo[i]];

    for( int j = 0; j < outbsPrimaryInfo.size(); j++ )
        prob *= distOutbsAndDeaths[outbsPrimaryInfo[j]][deathsPrimaryInfo[j]];

    return prob;
};




