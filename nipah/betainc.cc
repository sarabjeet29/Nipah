#include "betainc.h"

void nrerror(const char error_text[])
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
 	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
};

cdouble lngamma( cdouble val )
{
	double real_val = real( val );
	double imag_val = imag( val );
	gsl_sf_result result_logradius, result_phi;
    gsl_sf_lngamma_complex_e( real_val, imag_val, &result_logradius, &result_phi);
    return cdouble ( result_logradius.val, result_phi.val );
};

cdouble gamma( cdouble val )
{
	cdouble mylngamma = lngamma( val );
	double radius = exp( real( mylngamma ) ); 
	double theta  = imag( mylngamma );
	cdouble mygamma = polar( radius, theta );
	return mygamma;
};

cdouble betai( cdouble a, cdouble b, cdouble x )
{
	cdouble bt, mybeta;
	
	if (x == cdouble (0,0) || x == cdouble (1,1)) 
		bt = cdouble (0,0);
	else
		bt=exp(a*log(x)+b*log(cdouble(1,0) -x)); 

	if (real(x) < real((a+cdouble(1,0))/(a+b+cdouble(2,0)))) 
		return bt*betacf(a,b,x)/a;
	else
	{
		mybeta = exp( lngamma(a) + lngamma(b) - lngamma(a+b) );
		return mybeta - bt*betacf(b,a,1.0-x)/b; 
	};	
};

cdouble upperbetai( cdouble a, cdouble b, cdouble x )
{
	if( real(a) < -1 )
		return ((a+b)*upperbetai(a+cdouble(1,0),b,x) - pow(x,a)*pow((cdouble(1,0)-x),b))/a;

 	if( real(b) < 0 )
 		return (a+b)*upperbetai(a,b+cdouble(1,0),x)/b + exp(a*log(x)+b*log(cdouble(1,0) -x))/b;

	cdouble bt, mybeta;
	mybeta = exp( lngamma(a) + lngamma(b) - lngamma(a+b) );
		
	if (x == cdouble (0,0) || x == cdouble (1,1)) 
		bt = cdouble (0,0);
	else
		bt=exp(a*log(x)+b*log(cdouble(1,0) -x)); 

	if (real(x) < real((a+cdouble(1,0))/(a+b+cdouble(2,0)))) 
		return mybeta - bt*betacf(a,b,x)/a;
	else
		return bt*betacf(b,a,1.0-x)/b; 
};

cdouble betacf( cdouble a, cdouble b, cdouble x )
{
	int n;
	cdouble aa,c,d,del,h,qab,qam,qap,m,m2;

	qab = a + b;
	qap = a + cdouble(1,0);
	qam = a - cdouble(1,0);
	c = cdouble(1,0);
	d = cdouble(1,0) - qab*x/qap;
	if ( abs(d) < FPMIN) 
		d=FPMIN; 
	d = cdouble(1,0) / d;
	h = d;
	for (n = 1; n <= MAXIT; n++) 
	{
		m  = cdouble(n,0);
		m2 = cdouble(2,0)*m;
		aa = m*(b-m)*x/((qam+m2)*(a+m2)); 
		d = cdouble(1,0) + aa*d;
		if (abs(d) < FPMIN) 
			d=FPMIN; 
		c = cdouble(1,0) + aa/c;
		if (abs(c) < FPMIN) 
			c=FPMIN;
		d = cdouble(1,0)/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d = cdouble(1,0)+aa*d;
		if (abs(d) < FPMIN) 
			d=FPMIN; 
		c = cdouble(1,0) + aa/c;
		if (abs(c) < FPMIN) 
			c=FPMIN; 
		d = cdouble(1,0)/d;
		del = d*c;
		h *= del;
		if (abs(del-1.0) < EPS) 
			break;
	};
	
	if (n > MAXIT)
		nrerror("a or b too big, or MAXIT too small in betacf"); 

	return h;
};
