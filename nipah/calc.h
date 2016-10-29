#include "fft.h"

void getDistOutbs( double, double, vector<double>&, int, double );
void getDistSumPrimaryOutbs(double, double, vector<int>&, vector<double>&,
                       int, int, double, vector<double>&);
double getProbObsOutbs(double, double, vector<int>&, vector<int>&, int, int, int, double);
void getPmPrimaryInfo(double, double, vector<int>&, int, int, int, int, double,
                      vector<double>&);
void getPmNoPrimaryInfo(double, double, int, double, vector<vector<double>>&);
void getDistOutbsAndDeaths(double, double, double, double, vector<vector<double>>&, int,double);
void getDistSumPrimaryOutbsAndDeaths(double, double, double, double, vector<int>&, vector<int>&, vector<vector<double>>&, int, int, double, vector<vector<double>>&);
double getProbObsOubsAndDeaths(double, double, double, double, vector<int>&,
                               vector<int>&, vector<int>&, vector<int>&, int, int,
                               int, int, double);






