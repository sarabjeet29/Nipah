#include "calc.h"
#include <iostream>

void getDistOutbs(int N, vector<cdouble>& cachedGxy, vector<double>& distOutbs)
{
    vector<cdouble> fftInput(cachedGxy.begin(), cachedGxy.begin() + N);
    vector<cdouble> fftOutput;
	runfft1D( N, fftInput, fftOutput );
    
    for(int i = 0; i < N; i++)
        distOutbs.push_back(real(fftOutput[i]));
};

void getDistSumPrimaryOutbs(vector<int>& outbsPrimaryInfo, vector<double>& distOutbs, int N, int M,
                            vector<cdouble>& cachedGxy, vector<double>& distSumOutbSizesPrimary)
{
    cdouble tempval;
    vector<cdouble> outerfftInput, outerfftOutput, innerfftOutput;
    for( int i = 0; i < M; i++ )
    {
        vector<cdouble> innerfftInput(cachedGxy.begin() + i * N, cachedGxy.begin() + (i + 1) * N);
        runfft1D( N, innerfftInput, innerfftOutput );
        tempval = cdouble(1,0);
        for(int j : outbsPrimaryInfo)
            tempval *= innerfftOutput[j] / cdouble(distOutbs[j],0);

        outerfftInput.push_back(tempval);
        innerfftInput.clear();
        innerfftOutput.clear();
    };
    runfft1D( M, outerfftInput, outerfftOutput);

    for(int i = 0; i < M; i++)
        distSumOutbSizesPrimary.push_back(real(outerfftOutput[i]));
};

double getProbObsOutbs(double alpha, double nu, vector<int>& outbsPrimaryInfo,
                       vector<int>& outbsNoPrimaryInfo, int sumPrimaryOutbs, int N,
                       int M, vector<cdouble>& cachedGxy)
{
    vector<double> distOutbs, distSumOutbSizesPrimary;
    getDistOutbs(N, cachedGxy, distOutbs);
    getDistSumPrimaryOutbs(outbsPrimaryInfo, distOutbs, N, M, cachedGxy, distSumOutbSizesPrimary);

    double prob = distSumOutbSizesPrimary[sumPrimaryOutbs];
    
    for(int i : outbsNoPrimaryInfo)
        prob *= distOutbs[i];
    
    for(int j : outbsPrimaryInfo)
        prob *= distOutbs[j];
    
    return prob;
};

void getPmPrimaryInfo(vector<int>& outbsPrimaryInfo, int sumPrimaryOutbs, int targetOutbSize, int N, int M, vector<cdouble>& cachedGxy, vector<double>& Pmn)
{
    int lenPmn = M;
    vector<double> distOutbs, distSumPrimary;
    getDistOutbs(N, cachedGxy, distOutbs);
    getDistSumPrimaryOutbs(outbsPrimaryInfo, distOutbs, N, M, cachedGxy, distSumPrimary);
    
    vector<vector<cdouble>> fftInputL1, fftOutputL1;
    vector<cdouble> fftOutputL2, ctempRow;
    cdouble ctempVal;
    int scaledIndex, combinedIndex;
    
    for( int i = 0; i < M; i++ )
    {
        vector<cdouble> fftInputL2(cachedGxy.begin() + i * N, cachedGxy.begin() + (i + 1) * N);
        runfft1D( N, fftInputL2, fftOutputL2 );
        
        ctempVal = cdouble(1,0);
        for(int s : outbsPrimaryInfo)
        {
            if( s != targetOutbSize )
                ctempVal *= fftOutputL2[s] / distOutbs[s];
        };
            
        if(targetOutbSize == 12) //12 is repeated so adding it back
            ctempVal *= fftOutputL2[12] / distOutbs[12];
            
        fftInputL2.clear();
        fftOutputL2.clear();

        for( int k = 0; k < lenPmn; k++)
        {
            scaledIndex = k * M / lenPmn;
            if( i + scaledIndex < M)
                combinedIndex = i + scaledIndex;
            else
                combinedIndex = i + scaledIndex - M;

            vector<cdouble> fftInputL2(cachedGxy.begin() + combinedIndex * N,
                                       cachedGxy.begin() + (combinedIndex + 1) * N);
            runfft1D( N, fftInputL2, fftOutputL2 );
            ctempRow.push_back(ctempVal * fftOutputL2[targetOutbSize] / distOutbs[targetOutbSize]);
            fftInputL2.clear();
            fftOutputL2.clear();
        };
        
        fftInputL1.push_back(ctempRow);
        ctempRow.clear();
    };
    runfft2D( M, lenPmn, fftInputL1, fftOutputL1);
    
    for(int j = 0; j < targetOutbSize + 1; j++)
        Pmn.push_back(real(fftOutputL1[sumPrimaryOutbs][j]) / distSumPrimary[sumPrimaryOutbs]);
};

void getPmNoPrimaryInfo(int M, int N, vector<cdouble> cachedGxy, vector<vector<double>>& Pmn)
{
    vector<double> distOutbs, temprow;
    getDistOutbs(N, cachedGxy, distOutbs);
    vector<vector<cdouble>> fftInput, fftOutput;
    for( int i = 0; i < M; i++ )
    {
        vector<cdouble> ctemprow(cachedGxy.begin() + i * N, cachedGxy.begin() + (i + 1) * N);
        fftInput.push_back(ctemprow);
        ctemprow.clear();
    };
    runfft2D( M, N, fftInput, fftOutput);

    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < min(M, j + 1); i++)
            temprow.push_back(real(fftOutput[i][j])/distOutbs[j]);

        Pmn.push_back(temprow);
        temprow.clear();
    };
};

void getDistOutbsAndDeaths(double alpha, double nu, double p, double q,
                           vector<vector<double>>& distOutbsAndDeaths, int N,
                           double R )
{
    vector<double> temprow;
    vector<vector<cdouble>> fftInput, fftOutput;

    getfftInput2D_Hxzuw( cdouble(1,0), cdouble(1,0), cdouble(alpha,0), cdouble(nu,0), cdouble(p,0), cdouble(q,0), R, N, fftInput);
    
    runfft2D( N, N, fftInput, fftOutput);
    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            temprow.push_back(real(fftOutput[i][j]));

        distOutbsAndDeaths.push_back(temprow);
        temprow.clear();
    };
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
    vector<vector<cdouble>> outerfftInput, outerfftOutput, innerfftInput,innerfftOutput;

    for( int i = 0; i < M; i++ )
    {
        cout << i <<" ";
        for( int j = 0; j < M; j++ )
        {
            getfftInput2D_Hxzuw( cdouble(R,0) * exp(cdouble(0,2*PI*i/M)), cdouble(R,0) * exp(cdouble(0,2*PI*j/M)), alpha, nu, p, q, R, N, innerfftInput );

            runfft2D( N, N, innerfftInput, innerfftOutput );
            
            ctempval = cdouble(1,0);
            for(int k = 0; k < outbsPrimaryInfo.size(); k++)
            {
                outbSize = outbsPrimaryInfo[k];
                nDeaths = deathsPrimaryInfo[k];
                ctempval *= innerfftOutput[outbSize][nDeaths] /distOutbsAndDeaths[outbSize][nDeaths];
            };
            ctempRow.push_back(ctempval);
            innerfftInput.clear();
            innerfftOutput.clear();
        };
        outerfftInput.push_back(ctempRow);
        ctempRow.clear();
    };
    cout << "\n";
    runfft2D( M, M, outerfftInput, outerfftOutput);
    
    for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < M; j++)
            tempRow.push_back(real(outerfftOutput[i][j]));
        
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




