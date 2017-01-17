#include "calc.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/stat.h>
#include <zlib.h>

void cache_Gxy(vector<double>& alpha_vec, vector<double>& nu_vec, double R, int M, int N);
void cache_Hxzuw(vector<double>& alpha_vec, vector<double>& nu_vec, vector<double>& p_vec,
                 vector<double> q_vec, double R, double M, double N);

void save_prob_alpha_nu(vector<double>& alpha_vec, vector<double>& nu_vec,
                        vector<int>& outbsPrimaryInfo, vector<int>& outbsNoPrimaryInfo,
                        int sumPrimaryOutbs, int N, int M);

void save_PmnPrimaryInfo(vector<double> alpha_vec, vector<double> nu_vec,
                         vector<int> outbsPrimaryInfo, int sumPrimaryOutbs,
                         int M, int N);

void save_PmnNoPrimaryInfo(vector<double> alpha_vec, vector<double> nu_vec,
                           vector<int> outbsNoPrimaryInfo, int M, int N);

void save_prob_p_q(vector<double>& alpha_vec, vector<double>& nu_vec,
                   vector<double>& p_vec, vector<double>& q_vec,
                   vector<int>& outbsPrimaryInfo, vector<int>& outbsNoPrimaryInfo,
                   vector<int>& deathsPrimaryInfo, vector<int>& deathsNoPrimaryInfo,
                   int sumPrimaryOutbs, int sumPrimaryDeaths, int N, int M, double R );