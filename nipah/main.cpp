//
//  main.cpp
//  nipah
//
//  Created by Sarabjeet Singh on 9/19/16.
//  Copyright © 2016 Sarabjeet Singh. All rights reserved.
//

#include "calc.h"
#include <fstream>

int main( void )
{
    ofstream outputFile;
    double alpha = 0.9, nu = 0.900001, R = 1.0;
    
    static const int arr1[] = {5,4,7,3,1,16,44,12,66};
    vector<int> outbsNoPrimaryInfo (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]));

    static const int arr2[] = {5,4,5,0,1,14,40,10,45};
    vector<int> deathsNoPrimaryInfo (arr2,  arr2 + sizeof(arr2) / sizeof(arr2[0]));
    
    static const int arr3[] = {13,12,31,36,12,7,8,3};
    vector<int> outbsPrimaryInfo (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]));
    
    static const int arr4[] = {9,8,23,27,11,3,5,1};
    vector<int> deathsPrimaryInfo (arr4, arr4 + sizeof(arr4) / sizeof(arr4[0]));
    
    int N = 1000, M = 100;
    double *distOutbSizes, *distSumOutbSizesPrimary;
    distOutbSizes = new double[N];
    distSumOutbSizesPrimary = new double[M];
    getDistOutbSizes(alpha, nu, distOutbSizes, N, R);
    getDistSumPrimary(alpha, nu, distOutbSizes, outbsPrimaryInfo, distSumOutbSizesPrimary, N, M, R);
    
    outputFile.open("output.txt");
    for( int i = 0; i < N; i++ )
        outputFile << distOutbSizes[i]<<"\n";
    outputFile.close();
    
    outputFile.open("output2.txt");
    for( int i = 0; i < M; i++ )
        outputFile << distSumOutbSizesPrimary[i]<<"\n";
    outputFile.close();
    
    delete distOutbSizes;
    delete distSumOutbSizesPrimary;
    
    return 0;
};
