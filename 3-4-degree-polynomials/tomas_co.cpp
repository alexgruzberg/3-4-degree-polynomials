//Tomas Co 10/07/2014
//--------------------------//
//---ONLY FOR REAL ROOTS---//
//------------------------//
//https://www.academia.edu/23308251/Real_Roots_of_Cubic_Equation

#include "polynomials.h"
//	x^3 + a x^2 + b x + c

template<typename T>
std::vector<T> tomas_co(third_degree_polynomial<T> P)
{
	float one_third = 1.0 / 3.0; float one_twentyseventh = one_third * one_third * one_third;

	std::vector<float> coefs = P.get_coefs();
	float d = coefs[0], c = coefs[1], b = coefs[2];

	float bb = b * b;
	float p = (3 * c - bb) * one_third;
	float q = (b * (2 * bb - 9 * c) + 27 * d) * one_twentyseventh;

	std::vector<T> est_roots(3);

	float p_one_third = p * one_third;
	if (p_one_third > 0)
		throw sqrt_of_negative_number();

	float A = 2 * sqrt(-p_one_third); //p is less than 0

	float acos_arg = (3 * q) / (A * p);
	if (isinf(acos_arg) || isnan(acos_arg))
		throw division_by_zero();
	if (acos_arg > 1 || acos_arg < -1)
		throw invalid_types_of_complex();

	float phi = acos(acos_arg);
	float B = -b * one_third;

	for (int i = 0; i < 3; ++i)
		est_roots[i] = A * cos((phi + 2*i*M_PI) * one_third) + B;
	
	return est_roots;
}

template<typename T>
std::vector<std::complex<T>> tomas_co(third_degree_polynomial<std::complex<T>> P)
{
	float one_third = 1.0 / 3.0; float one_twentyseventh = one_third * one_third * one_third;

	std::vector<float> coefs = P.get_coefs();
	float d = coefs[0], c = coefs[1], b = coefs[2];

	float bb = b * b;
	float p = (3 * c - bb) * one_third;
	float q = (b * (2 * bb - 9 * c) + 27 * d) * one_twentyseventh;

	std::vector<std::complex<T>> est_roots(3);
	if (p < 0)
	{
		float p_one_third = p * one_third;

		float A = 2 * sqrt(-p_one_third);

		float acos_arg = (3 * q) / (A * p);

		if (isinf(acos_arg) || isnan(acos_arg))
			throw division_by_zero();
		if (acos_arg > 1 || acos_arg < -1)
			throw invalid_types_of_complex();

		float phi = acos(acos_arg);
		float B = -b * one_third;
		for (int i = 0; i < 3; ++i)
			est_roots[i] = A * cos((phi + 2 * i * M_PI) * one_third) + B;
	}
	else if (p > 0)
	{
		float A = 2 * sqrt(p * one_third);

		float asinh_arg = (3 * q) / (A * p);
		if (isinf(asinh_arg) || isnan(asinh_arg))
			throw division_by_zero();

		float phi = asinh(asinh_arg);
		est_roots[0] = est_roots[1] = std::complex<T>(0,0);
		est_roots[2] = std::complex<T>(-3 * A * 0.5 * sinh(phi * one_third),0);
	}
	return est_roots;
}