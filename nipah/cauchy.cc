//
//  cauchy.cpp
//  nipah
//
//  Created by Sarabjeet Singh on 5/25/17.
//  Copyright © 2017 Sarabjeet Singh. All rights reserved.
//

#include "cauchy.h"

int integrandFunc(unsigned ndim, const double *x, void *fdata, unsigned fdim,
                        double *fval)
{
    cdouble *params = (cdouble *) fdata;
    cdouble alpha, nu, p, q, outbsize, ndeaths, param_primaryoutbsize, param_primarydeaths;
    alpha = params[0];
    nu = params[1];
    p = params[2];
    q = params[3];
    outbsize = params[4];
    ndeaths = params[5];
    param_primaryoutbsize = params[6];
    param_primarydeaths = params[7];
    cdouble funcval = exp(-outbsize * cdouble(0, x[0]) - ndeaths * cdouble(0, x[1])) *
    Hxzuw(param_primaryoutbsize, exp(cdouble(0, x[0])), param_primarydeaths, exp(cdouble(0, x[1])), alpha, nu, p, q) / cdouble(4 * pow(PI, 2), 0);

    fval[0] = real(funcval);
    fval[1] = imag(funcval);
    return 0; // success
};

cdouble cauchyintegral(cdouble x, cdouble u, int outbsize, int ndeaths, cdouble alpha, cdouble nu, cdouble p, cdouble q)
{
    double xmin[2] = {0, 0}, xmax[2] = {2*PI, 2*PI};
    double tol = 1e-5, fval[2], err[2];
    unsigned dim = 2, fdim = 2;
    
    cdouble fdata[8] = {alpha, nu, p, q, cdouble(outbsize,0), cdouble(ndeaths,0), x, u};
    hcubature(fdim, integrandFunc, fdata, dim, xmin, xmax, 0, 0, tol, ERROR_INDIVIDUAL, fval,             err);
    return cdouble(fval[0], fval[1]);
};

//x: primary outbreak size
//z: total outbreak size
//u: primary deaths
//w: total deaths

//hcubature_v(unsigned fdim, int f, void *fdata,
//            unsigned dim, const double *xmin, const double *xmax,
//            unsigned maxEval, double reqAbsError, double reqRelError,
//            error_norm norm, double *val, double *err);



//double xmin[3] = {-2,-2,-2}, xmax[3] = {2,2,2}, sigma = 0.5, val, err;
//hcubature(1, f, &sigma, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
//printf("Computed integral = %0.10g +/- %g\n", val, err);
