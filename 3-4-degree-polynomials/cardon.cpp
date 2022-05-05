//CARDON’S METHOD TO SOLVE A CUBIC EQUATION

#include "polynomials.h"
#include <complex>
using namespace std::complex_literals;
//	x^3 + 3 a x^2 + 3 b x + c

template<typename T>
std::vector<T> third_degree_polynomial<T>::cardon()
{
	float a = this->coefs[2] / 3; float b = this->coefs[1] / 3;
	float H = b - pow(a,2); float G = 2 * pow(a, 3) - 3 * a * b + this->coefs[0];
	float delta = pow(G, 2) + 4 * pow(H, 3);
	std::complex<double> w(-0.5, sqrt(3) / 2);
	std::vector<float> est_roots(3);
	if (delta > 0)	//corresponding roots are all different
	{
		std::complex<double> H_c = (H, 0);
		double re_phi_c = pow((G + sqrt(delta)), (float)1 / (float)3) / pow((float)2, (float)1 / (float)3);
		std::complex<double> phi_c((double)re_phi_c, 0);
		est_roots[0] = -phi_c.real() + H / phi_c.real();
		est_roots[1] = (-w * phi_c + pow(w, 2) * H_c/ phi_c).real();
		est_roots[2] = (-pow(w,2) * phi_c + w * H_c / phi_c).real();
	}
	else if (delta == 0)
	{
		float phi = pow((G + sqrt(delta)) / 2, 1 / 3);
		est_roots[0] = -2 * phi;
		est_roots[1] = est_roots[2] = phi;
	}
	else //delta < 0
	{
		std::complex<double> H_c(H, 0);
		float imag_delta = sqrt(abs(delta)) / 2;
		std::complex<double> phi_c(G / 2, imag_delta);
		phi_c = pow(phi_c, 1.0/3);
		est_roots[0] = (-phi_c + H_c / phi_c).real();
		est_roots[1] = (-w * phi_c + pow(w, 2) * H_c / phi_c).real();
		est_roots[2] = (-pow(w, 2) * phi_c + w * H_c / phi_c).real();
	}
	for (int i = 0; i < 3; ++i)
		est_roots[i] = est_roots[i] - a;
	return est_roots;
}