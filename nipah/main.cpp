//
//  main.cpp
//  nipah
//
//  Created by Sarabjeet Singh on 9/19/16.
//  Copyright © 2016 Sarabjeet Singh. All rights reserved.
//

#include "calc.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/stat.h>

void save_prob_alpha_nu(vector<double>& alpha_vec, vector<double>& nu_vec,
                        vector<int>& outbsPrimaryInfo, vector<int>& outbsNoPrimaryInfo,
                        int sumPrimaryOutbs, int N, int M, double R )
{
    ofstream outputFile;
    double pij;
    vector< vector<double> > probObs;
    
    outputFile.open("./results/temp_prob_alpha_nu.txt");
    for(double alpha : alpha_vec)
    {
        vector<double> probRow;
        for( double nu : nu_vec)
        {
            cout << "alpha: " << alpha <<", nu: "<<nu<<"\n";
            if( abs(round(nu / alpha ) - nu / alpha) < 1e-12 )
            {
                cout <<"adjusted for alpha: "<<alpha<<", nu: "<<nu<< "\n";
                nu += 1e-6;
            };
            pij = getProbObs(alpha, nu, outbsPrimaryInfo, outbsNoPrimaryInfo, sumPrimaryOutbs, N, M, R );
            probRow.push_back(pij);
            outputFile << pij << ' ';
        };
        probObs.push_back(probRow);
        outputFile << "\n";
    };
    outputFile.close();
};

void save_PmnPrimaryInfo(vector<double> alpha_vec, vector<double> nu_vec,
                         vector<int> outbsPrimaryInfo, int sumPrimaryOutbs,
                         int N, int M, double R )
{
    ostringstream dirStringStream, fileStringStream;
    ofstream myfile;
    vector<double> Pmn;
    const char * filepath;
    const char * dirpath;
    double alpha = 0.7, nu = 2.0;

    dirStringStream << "./results/alpha=" << int(alpha*100)/100.0 << ",nu=" << int(nu*100)/100.0;
    dirpath = dirStringStream.str().c_str();
    mkdir(dirpath,0777);
    
    dirStringStream << "/p";
    dirpath = dirStringStream.str().c_str();
    mkdir(dirpath,0777);

    for(int s: outbsPrimaryInfo)
    {
        fileStringStream << dirStringStream.str().c_str() << "/" << s <<".txt";
        filepath = fileStringStream.str().c_str();
        fileStringStream.str(string());
        myfile.open ( filepath );
        getPmPrimaryInfo(alpha, nu, outbsPrimaryInfo, sumPrimaryOutbs, s, N, M, R, Pmn);
        for( double pij : Pmn)
            myfile << pij <<'\n';
        Pmn.clear();
        myfile.close();
    };
};

void save_PmnNoPrimaryInfo(vector<double> alpha_vec, vector<double> nu_vec,
                         vector<int> outbsNoPrimaryInfo, int N, double R )
{
    ostringstream dirStringStream, fileStringStream;
    ofstream myfile;
    const char * filepath;
    const char * dirpath;
    vector<vector<double>> Pmn;
    
    double alpha = 0.52, nu = 1.25;
    
    dirStringStream << "./results/alpha=" << int(alpha*100)/100.0 << ",nu=" << int(nu*100)/100.0;
    dirpath = dirStringStream.str().c_str();
    mkdir(dirpath,0777);

    dirStringStream << "/np";
    dirpath = dirStringStream.str().c_str();
    mkdir(dirpath,0777);

    getPmNoPrimaryInfo(alpha, nu, N, R, Pmn);
    for( int s: outbsNoPrimaryInfo)
    {
        fileStringStream << dirStringStream.str().c_str() << "/" << s <<".txt";
        filepath = fileStringStream.str().c_str();
        fileStringStream.str(string());
        myfile.open ( filepath );

        for( double pij : Pmn[s])
            myfile << pij <<'\n';

        myfile.close();
    };
};

int main( void )
{
    int N, M, sumPrimaryOutbs;
    double R, alpha_min, alpha_max, alpha_step, nu_min, nu_max, nu_step;
    
    N = 300;
    M = 100;
    R = 1.0;

    sumPrimaryOutbs = 60;
    
    alpha_min = 0.38;
    alpha_max = 0.66;
    alpha_step = 0.01;
    
    nu_min = 0.6;
    nu_max = 2.0;
    nu_step = 0.05;

    vector<double> distOutbSizes, distSumOutbSizesPrimary;

    static const int arr1[] = {5,4,7,3,1,16,44,12,66};
    vector<int> outbsNoPrimaryInfo (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]));

    static const int arr2[] = {5,4,5,0,1,14,40,10,45};
    vector<int> deathsNoPrimaryInfo (arr2,  arr2 + sizeof(arr2) / sizeof(arr2[0]));
    
    static const int arr3[] = {13,12,31,36,12,7,8,3};
    vector<int> outbsPrimaryInfo (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]));
    
    static const int arr4[] = {9,8,23,27,11,3,5,1};
    vector<int> deathsPrimaryInfo (arr4, arr4 + sizeof(arr4) / sizeof(arr4[0]));

    vector<double> alpha_vec, nu_vec;
    
    for( double alpha = alpha_min; alpha < alpha_max + alpha_step; alpha += alpha_step)
        alpha_vec.push_back(alpha);

    for( double nu = nu_min; nu < nu_max + nu_step; nu += nu_step)
        nu_vec.push_back(nu);
    
    save_prob_alpha_nu(alpha_vec, nu_vec, outbsPrimaryInfo, outbsNoPrimaryInfo,
                       sumPrimaryOutbs, N, M, R );
    
    save_PmnPrimaryInfo(alpha_vec, nu_vec, outbsPrimaryInfo, sumPrimaryOutbs,
                        N, M, R );

    save_PmnNoPrimaryInfo(alpha_vec, nu_vec, outbsNoPrimaryInfo, N, R);
    return 0;
};
