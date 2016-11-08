#include "pgf.h"
#include <iostream>

cdouble lambda0( cdouble y, cdouble alpha )
{
	return (alpha + cdouble(1,0) - sqrt( pow(alpha+cdouble(1,0),2) - 	
	cdouble(4,0)*alpha*y))/(cdouble(2,0)*alpha);
};

cdouble Gxy( cdouble x, cdouble y, cdouble alpha, cdouble nu )
{
    cout << "x: " << x << ", y: " << y << '\n';
    if( x == cdouble(1,0) && y == cdouble(1,0) )
        return cdouble(1,0);
    
    else if( x == cdouble(1,0) && alpha == cdouble(1,0) && nu == cdouble(1,0))
    {
        return cdouble(1,0) + sqrt(y) / log((sqrt(cdouble(1,0) - y)) / (cdouble(1,0) + sqrt(y)));
    }
    else
    {
        cdouble L0, L1, L0byL1, a, b, z, incbeta;
        L0 = lambda0( y, alpha );
        L1 = cdouble(1,0) + cdouble(1,0)/alpha - L0;
        L0byL1 = L0 / L1;
        a = cdouble(1,0) - nu*x/alpha;
        b = nu/alpha*((cdouble(1,0) - L0*x)/(L1 - L0));
        z = cdouble(1,0) - L0byL1;
        incbeta = upperbetai( a, b, z );
        return cdouble(1,0) - alpha*L1*pow(z,a)*pow(L0byL1,b)/(nu*incbeta);
    };
};

cdouble Hxzuw( cdouble x, cdouble z, cdouble u, cdouble w, cdouble alpha, cdouble nu, cdouble p, cdouble q)
{
	return Gxy( x * (cdouble(1,0) - p + p*u*w) / (cdouble(1,0) - q + q*w), (cdouble(1,0) - q + q*w) * z, 
	alpha, nu );
};



