#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/stat.h>
#include <zlib.h>
#include "betainc.h"
#include <iterator>

void read_cached_Gxy(double alpha, double nu, int M, int N, vector<cdouble>& Gxy);
void read_cached_Hzw(double alpha, double nu, double p, double q, int N,
                     vector<cdouble>& cachedHzw);
void read_cached_Hxzuw(double alpha, double nu, double p, double q, int M, int N,
                       vector<cdouble>& Hxzuw);
