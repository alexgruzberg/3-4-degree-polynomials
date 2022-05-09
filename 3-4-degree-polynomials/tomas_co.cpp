//Tomas Co 10/07/2014
//---ONLY FOR REAL ROOTS---//

#include "polynomials.h"
//	x^3 + a x^2 + b x + c

template<typename T>
std::vector<T> tomas_co(third_degree_polynomial<T> P)
{
	float one_third = 1.0 / 3.0; float one_twentyseventh = 1.0 / 27.0;

	std::vector<float> coefs = P.get_coefs();
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

template<typename T>
std::vector<std::complex<T>> tomas_co(third_degree_polynomial<std::complex<T>> P)
{
	float one_third = 1.0 / 3.0; float one_twentyseventh = 1.0 / 27.0; float half = 1.0 / 2.0;

	std::vector<float> coefs = P.get_coefs();
	float p = (3 * coefs[1] - coefs[2] * coefs[2]) * one_third;
	float q = (2 * pow(coefs[2], 3) - 9 * coefs[2] * coefs[1] + 27 * coefs[0]) * one_twentyseventh;
	std::vector<std::complex<T>> est_roots(3);
	if (p < 0)
	{
		float A = 2 * sqrt(-p * one_third); //p is less than 0
		float phi = acos((3 * q) / (A * p));
		float B = -coefs[2] * one_third;
		for (int i = 0; i < 3; ++i)
			est_roots[i] = A * cos((phi + 2 * i * M_PI) * one_third) + B;
	}
	else if (p > 0)
	{
		float A = 2 * sqrt(p * one_third);
		float phi = asinh(3 * q / (A * p));
		est_roots[0] = est_roots[1] = std::complex<T>(0,0);
		est_roots[2] = std::complex<T>(-3 * A * half * sinh(phi * one_third),0);
	}
	return est_roots;
}