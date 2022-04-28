//Tomas Co 10/07/2014

#include "polynomials.h"
//	x^3 + a x^2 + b x + c

std::vector<float> third_degree_polynomial::tomas_co()
{
	float p = (3 * coefs[1] - pow(coefs[2], 2))/3;
	float q = (2 * pow(coefs[2], 3) - 9 * coefs[2] * coefs[1] + 27 * coefs[0]) / 27;
	std::vector<float> est_roots(3);
	float A = 2 * sqrt(-p / 3); //p is less than 0
	float phi = acos((3 * q) / (A * p));
	float B = -coefs[2] / 3;
	for (int i = 0; i < 3; ++i)
		est_roots[i] = A * cos((phi + 2*i*M_PI)/3) + B;
	
	return est_roots;
}