//Tomas Co 10/07/2014
//https://www.academia.edu/23308251/Real_Roots_of_Cubic_Equation

//-----------------------------------------------------------//
//---Modified by Maria Ermakova to work with complex roots---//
//-----------------------------------------------------------//

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

	float A = 2 * sqrt(-p_one_third);

	float acos_arg = (3 * q) / (A * p);
	if (isinf(acos_arg) || isnan(acos_arg))
		throw division_by_zero();
	if (acos_arg > 1 || acos_arg < -1)
		throw invalid_types_of_complex();

	float third_phi = acos(acos_arg) * one_third;
	float B = -b * one_third;
	float two_third_pi = 2 * M_PI * one_third;

	for (int i = 0; i < 3; ++i)
		est_roots[i] = A * cos(third_phi + i * two_third_pi) + B;
	
	sort(est_roots.begin(), est_roots.end());

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

	float B = -b * one_third;

	std::vector<std::complex<T>> est_roots(3);
	if (p < 0)
	{
		float p_one_third = p * one_third;

		float A = 2 * sqrt(-p_one_third);

		float C = (3 * q) / (A * p);

		if (isinf(C) || isnan(C))
			throw division_by_zero();

		float two_third_pi = 2 * M_PI * one_third;

		if (abs(C) <= 1)	//all roots should be real
		{
			float third_phi = acos(C) * one_third;
			for (int i = 0; i < 3; ++i)
				est_roots[i] = A * cos(third_phi + i * two_third_pi) + B;
		}
		else if (C < -1)	//the polynomial has complex conjugate
		{
			float _log = -log(abs(C + sqrt(C * C - 1)));
			if (isinf(_log) || isnan(_log))
				throw nan_value();

			std::complex<T> phi(M_PI, _log);
			std::complex<float> third_phi = phi * one_third;
			for (int i = 0; i < 3; ++i)
				est_roots[i > 1 ? 0 : i + 1] = A * cos(third_phi + i * two_third_pi) + B;
		}
		else if (C > 1)		//the polynomial has complex conjugate
		{
			std::complex<T> phi(0, log(C + sqrt(C * C - 1)));
			std::complex<float>  third_phi = phi * one_third;
			for (int i = 0; i < 3; ++i)
				est_roots[2 - i] = A * cos(third_phi + i * two_third_pi) + B;
		}
	}
	else if (p > 0)		//the polynomial has complex conjugate
	{
		float A = 2 * sqrt(p * one_third);

		float asinh_arg = (3 * q) / (A * p);
		if (isinf(asinh_arg) || isnan(asinh_arg))
			throw division_by_zero();

		float phi = asinh(asinh_arg);
		phi = phi * one_third;
		std::complex<T> two_third_pi_i(0, 2 * M_PI * one_third);
		for (int i = 0; i < 3; ++i)
			est_roots[2-i] = -A * sinh(phi + std::complex<T>(i,0) * two_third_pi_i) + B;
	}
	return est_roots;
}