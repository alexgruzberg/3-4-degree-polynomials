//CARDON’S METHOD TO SOLVE A CUBIC EQUATION
#include "polynomials.h"
#include <complex>
using namespace std::complex_literals;
//	x^3 + 3 a x^2 + 3 b x + c

template<typename T>
std::vector<T> cardon(third_degree_polynomial<T> P)
{
	float one_third = 1.0 / 3.0;

	std::vector<float> coefs = P.get_coefs();
	float d = coefs[0], c = coefs[1], b = coefs[2];

	b = b * one_third; c = c * one_third;
	float bb = b * b;
	float H = c - bb; float G = b * (2 * bb - 3 * c) + d;
	float delta = G * G + 4 * H * H * H;

	std::complex<float> w(-0.5, sqrt(3) * 0.5);
	std::vector<T> est_roots(3);

	if (delta > 0)	//corresponding roots are all real and different
	{
		std::complex<float> H_c = (H, 0);
		float re_phi_c = cbrt(G + sqrt(delta)) / cbrt(2);
		std::complex<float> phi_c(re_phi_c, 0);

		float H_div_re_phi_c = H / re_phi_c;
		if (isnan(H_div_re_phi_c))
			throw division_by_zero();
		std::complex<float> H_div_phi = H_c / phi_c;
		if (isnan(abs(H_div_phi)))
			throw division_by_zero();

		est_roots[0] = -re_phi_c + H_div_re_phi_c;
		est_roots[1] = (w * (-phi_c + w * H_div_phi)).real();
		est_roots[2] = (w * (-w * phi_c + H_div_phi)).real();
	}
	else if (delta == 0)	//corresponding roots are all real and 2 of them are the same
	{
		float phi = cbrt(G * 0.5);
		est_roots[0] = -2 * phi;
		est_roots[1] = est_roots[2] = phi;
	}
	else    //there might be a pair of complex conjugates
	{
		std::complex<float> H_c(H, 0);
		float imag_delta = sqrt(abs(delta)) * 0.5;
		std::complex<float> phi_c(G * 0.5, imag_delta);
		phi_c = cbrt(phi_c);
		est_roots[0] = (-phi_c + H_c / phi_c).real();
		est_roots[1] = ( w * (-phi_c + w * H_c / phi_c)).real();
		est_roots[2] = ( w * (-w * phi_c + H_c / phi_c)).real();
	}
	for (int i = 0; i < 3; ++i)
		est_roots[i] = est_roots[i] - b;
	return est_roots;
}

template <typename T>
std::vector<std::complex<T>> cardon(third_degree_polynomial<std::complex<T>> P)
{
	float one_third = 1.0 / 3.0;

	std::vector<float> coefs = P.get_coefs();
	float d = coefs[0], c = coefs[1], b = coefs[2];

	b = b * one_third; c = c * one_third;
	float bb = b * b;
	float H = c - bb; float G = b * (2 * bb - 3 * c) + d;
	float delta = G * G + 4 * H * H * H;

	std::complex<float> w(-0.5, sqrt(3) * 0.5);
	std::vector<std::complex<T>> est_roots(3);

	if (delta > 0)	//corresponding roots are all real and different
	{
		std::complex<float> H_c(H, 0);
		float re_phi_c = cbrt(G + sqrt(delta)) / cbrt(2);
		std::complex<float> phi_c(re_phi_c, 0);

		std::complex<float> H_div_phi = H_c / phi_c;
		if (isnan(abs(H_div_phi)))
			throw division_by_zero();

		est_roots[0] = w * (-w * phi_c + H_div_phi);
		est_roots[1] = w * (-phi_c + w * H_div_phi);
		est_roots[2] = -phi_c + H_div_phi;
	}
	else if (delta == 0)	//corresponding roots are all real and 2 of them are the same
	{
		float phi = cbrt((G + sqrt(delta)) * 0.5);
		est_roots[0] = -2 * phi;
		est_roots[1] = est_roots[2] = phi;
	}
	else    //there might be a pair of complex conjugates
	{
		std::complex<float> H_c(H, 0);
		float imag_delta = sqrt(abs(delta)) * 0.5;
		std::complex<float> phi_c(G * 0.5, imag_delta);
		phi_c = pow(phi_c, one_third);		//there is no cbrt<complex<T>>

		std::complex<float> H_div_phi = H_c / phi_c;
		if (isnan(abs(H_div_phi)))
			throw division_by_zero();

		est_roots[0] = w * (-w * phi_c + H_div_phi);
		est_roots[1] = w * (-phi_c + w * H_div_phi);
		est_roots[2] = -phi_c + H_div_phi;
	}
	for (int i = 0; i < 3; ++i)
		est_roots[i] = est_roots[i] - b;
	return est_roots;
}