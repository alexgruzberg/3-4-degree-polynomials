//A simplified expression for the solution of cubic polynomial equations using function evaluation-Tiruneh-2020//
// https://arxiv.org/abs/2002.06976 //

#include "polynomials.h"
//	x^3 + a x^2 + b x + c

template<typename T>
std::vector<T> tiruneh(third_degree_polynomial<T> P)
{
	float one_third = 1.0 / 3.0; float half = 1.0 / 2.0;

	std::vector<T> coefs = P.get_coefs();
	float z = -coefs[2] * one_third;
	float Q = (-coefs[2] * coefs[2] * one_third + coefs[1]) * one_third;
	float R = -(pow(z, 3) + coefs[2] * z * z + coefs[1] * z + coefs[0]) * half;
	float theta = acos(R / sqrt(pow(-Q, 3)));
	std::vector<T> est_roots(3);
	float two_sqrt_q = 2 * sqrt(-Q);
	est_roots[0] = (two_sqrt_q * cos(theta * one_third) + z);
	est_roots[1] = (two_sqrt_q * cos((theta + 2 * M_PI) * one_third) + z);
	est_roots[2] = (two_sqrt_q * cos((theta + 4 * M_PI) * one_third) + z);
	return est_roots;
}