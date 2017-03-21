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

void getDistSumPrimaryOutbs(vector<int>& outbsPrimaryInfo, vector<double>& distOutbs, int M, int N,
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
                       vector<int>& outbsNoPrimaryInfo, int sumPrimaryOutbs, int M,
                       int N, vector<cdouble>& cachedGxy)
{
    vector<double> distOutbs, distSumOutbSizesPrimary;
    getDistOutbs(N, cachedGxy, distOutbs);
    getDistSumPrimaryOutbs(outbsPrimaryInfo, distOutbs, M, N, cachedGxy, distSumOutbSizesPrimary);

    double prob = distSumOutbSizesPrimary[sumPrimaryOutbs];
    
    for(int i : outbsNoPrimaryInfo)
        prob *= distOutbs[i];
    
    for(int j : outbsPrimaryInfo)
        prob *= distOutbs[j];
    
    return prob;
};

void getPmPrimaryInfo(vector<int>& outbsPrimaryInfo, int sumPrimaryOutbs, int targetOutbSize, int M, int N, vector<cdouble>& cachedGxy, vector<double>& Pmn)
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

void getDistOutbsAndDeaths(int M, int N, vector<cdouble> cachedHxzuw,
                           vector<vector<double>>& distOutbsAndDeaths)
{
    vector<cdouble> fftInput(cachedHxzuw.begin(), cachedHxzuw.begin() + M*M);
    vector<cdouble> fftOutput;
    runfft2DflatInput(N, N, fftInput, distOutbsAndDeaths);
};

void getDistSumPrimaryOutbsAndDeaths(vector<int>& outbsPrimaryInfo, vector<int>& deathsPrimaryInfo,
                                     vector<vector<double>>& distOutbsAndDeaths, int M, int N,
                                     vector<cdouble>& cachedHxzuw,
                                     vector<vector<double>>& distSumPrimaryOutbsAndDeaths)
{
    int outbSize, nDeaths;
    cdouble ctempval;
    vector<double> tempRow;
    vector<cdouble> outerfftInput, outerfftOutput, innerfftInput,innerfftOutput;
    int startIndex, endIndex, tempIndex;
    
    for( int i = 0; i < M; i++ )
    {
        cout << i <<" ";
        for( int j = 0; j < M; j++ )
        {
            startIndex = i * M * N * N + j * N * N;
            endIndex = (i * M + j + 1) * N * N;
            vector<cdouble> innerfftInput(cachedHxzuw.begin() + startIndex,
                                          cachedHxzuw.begin() + endIndex);

            runfft2DflatIO( N, N, innerfftInput, innerfftOutput );
            
            ctempval = cdouble(1,0);
            for(int k = 0; k < outbsPrimaryInfo.size(); k++)
            {
                outbSize = outbsPrimaryInfo[k];
                nDeaths = deathsPrimaryInfo[k];
                tempIndex = outbSize * M + nDeaths;
                ctempval *= innerfftOutput[tempIndex] /distOutbsAndDeaths[outbSize][nDeaths];
            };
            outerfftInput.push_back(ctempval);
            innerfftInput.clear();
            innerfftOutput.clear();
        };
    };
    cout << "\n";
    runfft2DflatInput( M, M, outerfftInput, distSumPrimaryOutbsAndDeaths);
};

double getProbObsOutbsAndDeaths(vector<int>& outbsPrimaryInfo,
                               vector<int>& outbsNoPrimaryInfo,
                               vector<int>& deathsPrimaryInfo,
                               vector<int>& deathsNoPrimaryInfo, int sumPrimaryOutbs,
                               int sumPrimaryDeaths, int M, int N,
                               vector<cdouble> cachedHxzuw)
{
    vector<vector<double>> distOutbsAndDeaths, distSumPrimaryOutbsAndDeaths;
    getDistOutbsAndDeaths(M, N, cachedHxzuw, distOutbsAndDeaths);
    getDistSumPrimaryOutbsAndDeaths(outbsPrimaryInfo, deathsPrimaryInfo, distOutbsAndDeaths, M, N,
                                    cachedHxzuw, distSumPrimaryOutbsAndDeaths);
    
    double prob = distSumPrimaryOutbsAndDeaths[sumPrimaryOutbs][sumPrimaryDeaths];
    
    for( int i = 0; i < outbsNoPrimaryInfo.size(); i++ )
        prob *= distOutbsAndDeaths[outbsNoPrimaryInfo[i]][deathsNoPrimaryInfo[i]];

    for( int j = 0; j < outbsPrimaryInfo.size(); j++ )
        prob *= distOutbsAndDeaths[outbsPrimaryInfo[j]][deathsPrimaryInfo[j]];

    return prob;
};




