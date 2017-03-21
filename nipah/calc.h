#include "fft.h"

void getDistOutbs(int N, vector<cdouble>& cachedGxy, vector<double>& distOutbs);

double getProbObsOutbs(double alpha, double nu, vector<int>& outbsPrimaryInfo,
                       vector<int>& outbsNoPrimaryInfo, int sumPrimaryOutbs, int M,
                       int N, vector<cdouble>& cachedGxy);

void getPmPrimaryInfo(vector<int>& outbsPrimaryInfo, int sumPrimaryOutbs, int targetOutbSize, int M, int N, vector<cdouble>& cachedGxy, vector<double>& Pmn);

void getPmNoPrimaryInfo(int M, int N, vector<cdouble> cachedGxy, vector<vector<double>>& Pmn);

void getDistOutbsAndDeaths(double, double, double, double, vector<vector<double>>&, int,double);

void getDistSumPrimaryOutbsAndDeaths(double, double, double, double, vector<int>&, vector<int>&, vector<vector<double>>&, int, int, double, vector<vector<double>>&);

//double getProbObsOutbsAndDeaths(vector<int>& outbsPrimaryInfo,
//                               vector<int>& outbsNoPrimaryInfo,
//                               vector<int>& deathsPrimaryInfo,
//                               vector<int>& deathsNoPrimaryInfo, int sumPrimaryOutbs,
//                               int sumPrimaryDeaths, int M, int N,
//                               vector<cdouble> cachedHxzuw);




