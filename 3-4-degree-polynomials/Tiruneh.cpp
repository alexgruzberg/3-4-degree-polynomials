//A simplified expression for the solution of cubic polynomial equations using function evaluation-Tiruneh-2020//
// https://arxiv.org/abs/2002.06976 //

#include "polynomials.h"
//	x^3 + a x^2 + b x + c

template<typename T>
std::vector<T> tiruneh(third_degree_polynomial<T> P)
{
	float one_third = 1.0 / 3.0;

	std::vector<T> coefs = P.get_coefs();
	float d = coefs[0], c = coefs[1], b = coefs[2];

	float z = -b * one_third;
	float Q = (-b * b * one_third + c) * one_third;
	float R = -(z * (z * (z + b) + c) + d) * 0.5;

	if (Q == 0 && R == 0)
	{
		float x = -b * one_third;
		std::vector<T> est_roots = { x, x, x };
		return est_roots;
	}

	float minus_Q_cube = -Q * Q * Q;

	if (isnan(minus_Q_cube))
		throw division_by_zero();
	if (minus_Q_cube < 0)
		throw sqrt_of_negative_number();

	float theta = acos(R / sqrt(minus_Q_cube));
	std::vector<T> est_roots(3);

	float two_sqrt_q = 2 * sqrt(-Q);
	est_roots[0] = (two_sqrt_q * cos(theta * one_third) + z);
	est_roots[1] = (two_sqrt_q * cos((theta + 2 * M_PI) * one_third) + z);
	est_roots[2] = (two_sqrt_q * cos((theta + 4 * M_PI) * one_third) + z);

	return est_roots;
}

template<typename T>
std::vector<std::complex<T>> tiruneh(third_degree_polynomial<std::complex<T>> P)
{
	float one_third = 1.0 / 3.0;

	std::vector<float> coefs = P.get_coefs();
	float d = coefs[0], c = coefs[1], b = coefs[2];

	float z = -b * one_third;
	float Q = (-b * b * one_third + c) * one_third;
	float R = -( z * ( z * (z + b) + c) + d) * 0.5;

	std::vector<std::complex<T>> est_roots(3);

	float Q_cube = Q * Q * Q;
	float D = Q_cube + R * R;
	if (D < 0)
	{
		if (isnan(Q_cube))
			throw division_by_zero();
		if (Q_cube > 0)
			throw sqrt_of_negative_number();

		float theta = acos(R / sqrt(-Q_cube));
		float two_sqrt_q = 2 * sqrt(-Q);
		est_roots[0] = (two_sqrt_q * cos(theta * one_third) + z);
		est_roots[1] = (two_sqrt_q * cos((theta + 2 * M_PI) * one_third) + z);
		est_roots[2] = (two_sqrt_q * cos((theta + 4 * M_PI) * one_third) + z);
	}
	else if (D > 0)
	{


		float B = cbrt(R + sqrt(D)) + cbrt(R - sqrt(D));
		float BB_plus_four_Q = B * B + 4 * Q;
		if (BB_plus_four_Q < 0 || D < 0)
			throw sqrt_of_negative_number();

		float root_real = z - 0.5 * B, root_imag = sqrt(3) * 0.5 * sqrt(BB_plus_four_Q);
		est_roots[0] = std::complex<T>(root_real, root_imag);
		est_roots[1] = std::complex<T>(root_real, -root_imag);
		est_roots[2] = z + B;
	}
	else
	{
		float x = -b * one_third;
		std::vector<std::complex<T>> est_roots = { x, x, x };
		return est_roots;
	}

	return est_roots;
}