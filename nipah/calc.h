#include "fft.h"

void getDistOutbSizes( double, double, vector<double>&, int, double );
void getDistSumPrimary(double, double, vector<int>&, vector<double>&,
                       int, int, double, vector<double>&);
double getProbObs(double, double, vector<int>&, vector<int>&, int, int, int, double);
void getPmPrimaryInfo(double, double, vector<int>&, int, int, int, int, double,
                      vector<double>&);
void getPmNoPrimaryInfo(double, double, int, double, vector<vector<double>>&);



