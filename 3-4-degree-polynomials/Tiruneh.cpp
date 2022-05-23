//A simplified expression for the solution of cubic polynomial equations using function evaluation-Tiruneh-2020//
// https://www.researchgate.net/publication/339325448_A_simplified_expression_for_the_solution_of_cubic_polynomial_equations_using_function_evaluation //

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

	if (Q > 0)
		throw sqrt_of_negative_number();

	float sqrt_minus_cube = sqrt(-Q);

	float arg = R / (-Q * sqrt_minus_cube);

	if (isinf(arg) || isnan(arg))
		throw division_by_zero();
	if (arg < -1 || arg > 1)
		throw arccos_out_of_range();

	float theta = acos(arg);
	std::vector<T> est_roots(3);

	float two_sqrt_q = 2 * sqrt_minus_cube;
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

	float D = Q * Q * Q + R * R;
	if (D < 0)
	{
		if (Q > 0)
			throw sqrt_of_negative_number();

		float sqrt_minus_cube = sqrt(-Q);

		float arg = R / (-Q * sqrt_minus_cube);

		if (isinf(arg) || isnan(arg))
			throw division_by_zero();
		if (arg < -1 || arg > 1)
			throw arccos_out_of_range();

		float theta = acos(arg);
		float two_sqrt_q = 2 * sqrt_minus_cube;
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