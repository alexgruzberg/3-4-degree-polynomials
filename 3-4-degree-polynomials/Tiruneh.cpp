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
	if (Q == 0 && R == 0)
	{
		std::vector<T> est_roots = { -coefs[2] * one_third, -coefs[2] * one_third, -coefs[2] * one_third };
		return est_roots;
	}
	float theta = acos(R / sqrt(pow(-Q, 3)));
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
	float one_third = 1.0 / 3.0; float half = 1.0 / 2.0;

	std::vector<float> coefs = P.get_coefs();
	float z = -coefs[2] * one_third;
	float Q = (-coefs[2] * coefs[2] * one_third + coefs[1]) * one_third;
	float R = -(pow(z, 3) + coefs[2] * z * z + coefs[1] * z + coefs[0]) * half;
	std::vector<std::complex<T>> est_roots(3);
	float D = Q * Q * Q + R * R;
	if (D < 0)
	{
		float theta = acos(R / sqrt(pow(-Q, 3)));
		float two_sqrt_q = 2 * sqrt(-Q);
		est_roots[0] = (two_sqrt_q * cos(theta * one_third) + z);
		est_roots[1] = (two_sqrt_q * cos((theta + 2 * M_PI) * one_third) + z);
		est_roots[2] = (two_sqrt_q * cos((theta + 4 * M_PI) * one_third) + z);
	}
	else if (D > 0)
	{
		float B = cbrt(R + sqrt(R * R + Q * Q * Q)) + cbrt(R - sqrt(R * R + Q * Q * Q));
		float root_real = z - half * B, root_imag = sqrt(3) * half * sqrt(B * B + 4 * Q);
		est_roots[0] = std::complex<T>(root_real, root_imag);
		est_roots[1] = std::complex<T>(root_real, -root_imag);
		est_roots[2] = z + B;
	}
	else
	{
		std::vector<std::complex<T>> est_roots = { -coefs[2] * one_third, -coefs[2] * one_third, -coefs[2] * one_third };
		return est_roots;
	}
	return est_roots;
}