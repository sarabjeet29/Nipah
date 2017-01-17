#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/stat.h>
#include "fft.h"

void read_cached_Gxy(double alpha, double nu, int M, int N, vector<cdouble>& Gxy);