//Tomas Co 10/07/2014

#include "polynomials.h"
//	x^3 + a x^2 + b x + c

template<typename T>
std::vector<T> tomas_co(third_degree_polynomial<T> P)
{
	float one_third = 1.0 / 3.0; float one_twentyseventh = 1.0 / 27.0;

	std::vector<T> coefs = P.get_coefs();
	float p = (3 * coefs[1] - coefs[2] * coefs[2]) * one_third;
	float q = (2 * pow(coefs[2], 3) - 9 * coefs[2] * coefs[1] + 27 * coefs[0]) * one_twentyseventh;
	std::vector<T> est_roots(3);
	float A = 2 * sqrt(-p * one_third); //p is less than 0
	float phi = acos((3 * q) / (A * p));
	float B = -coefs[2] * one_third;
	for (int i = 0; i < 3; ++i)
		est_roots[i] = A * cos((phi + 2*i*M_PI) * one_third) + B;
	
	return est_roots;
}