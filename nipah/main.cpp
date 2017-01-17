//
//  main.cpp
//  nipah
//
//  Created by Sarabjeet Singh on 9/19/16.
//  Copyright © 2016 Sarabjeet Singh. All rights reserved.
//

#include "cachedata.h"

int main( void )
{
    int N, M, sumPrimaryOutbs, sumPrimaryDeaths;
    double R, alpha_min, alpha_max, alpha_step, nu_min, nu_max, nu_step, p_min, p_max, p_step, q_min, q_max, q_step;
    
    N = 100;
    M = 100;
    R = 1.0;

    sumPrimaryOutbs = 60;
    sumPrimaryDeaths = 45;
    
    alpha_min = 0.38;
    alpha_max = 0.66;
    alpha_step = 0.01;
    
    nu_min = 0.6;
    nu_max = 2.0;
    nu_step = 0.05;
    
    p_min = 0.55;
    p_max = 0.9;
    p_step = 0.01;
    
    q_min = 0.6;
    q_max = 0.95;
    q_step = 0.01;

    vector<double> distOutbs, distSumPrimaryOutbs;
    vector<vector<double>> distOutbsAndDeaths, distSumPrimaryOutbsAndDeaths;

    static const int arr1[] = {5,4,7,3,1,16,44,12,66};
    vector<int> outbsNoPrimaryInfo (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]));

    static const int arr2[] = {5,4,5,0,1,14,40,10,45};
    vector<int> deathsNoPrimaryInfo (arr2,  arr2 + sizeof(arr2) / sizeof(arr2[0]));
    
    static const int arr3[] = {13,12,31,36,12,7,8,3};
    vector<int> outbsPrimaryInfo (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]));
    
    static const int arr4[] = {9,8,23,27,11,3,5,1};
    vector<int> deathsPrimaryInfo (arr4, arr4 + sizeof(arr4) / sizeof(arr4[0]));

    vector<double> alpha_vec, nu_vec, p_vec, q_vec;
    
    for( double alpha = alpha_min; alpha < alpha_max + alpha_step; alpha += alpha_step)
        alpha_vec.push_back(alpha);

    for( double nu = nu_min; nu < nu_max + nu_step; nu += nu_step)
        nu_vec.push_back(nu);

    for( double p = p_min; p < p_max + p_step; p += p_step)
        p_vec.push_back(p);

    for( double q = q_min; q < q_max + q_step; q += q_step)
        q_vec.push_back(q);
    
    cache_Gxy(alpha_vec, nu_vec, R, M, N);
  
    save_prob_alpha_nu(alpha_vec, nu_vec, outbsPrimaryInfo, outbsNoPrimaryInfo,
                       sumPrimaryOutbs, N, M);

    save_PmnPrimaryInfo(alpha_vec, nu_vec, outbsPrimaryInfo, sumPrimaryOutbs, M, N);

    save_PmnNoPrimaryInfo(alpha_vec, nu_vec, outbsNoPrimaryInfo, M, N);
    
    save_prob_p_q(alpha_vec, nu_vec, p_vec, q_vec, outbsPrimaryInfo, outbsNoPrimaryInfo,
                  deathsPrimaryInfo, deathsNoPrimaryInfo, sumPrimaryOutbs, sumPrimaryDeaths,
                  N, M, R );
    
    cache_Hxzuw(alpha_vec, nu_vec, p_vec, q_vec, R, M, N);

    return 0;
};
