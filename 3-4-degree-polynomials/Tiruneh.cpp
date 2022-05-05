//A simplified expression for the solution of cubic polynomial equations using function evaluation-Tiruneh-2020//
// https://arxiv.org/abs/2002.06976 //

#include "polynomials.h"
//	x^3 + a x^2 + b x + c

template<typename T>
std::vector<T> third_degree_polynomial<T>::tiruneh()
{
	float z = -this->coefs[2]/3;
	float Q = (-pow(this->coefs[2],2) / 3 + this->coefs[1]) / 3;
	float R = -(pow(z, 3) + this->coefs[2] * pow(z, 2) + this->coefs[1] * z + this->coefs[0]) / 2;
	float theta = acos(R / sqrt(pow(-Q, 3)));
	std::vector<T> est_roots(3);
	est_roots[0] = (2 * sqrt(-Q) * cos(theta / 3) + z);
	est_roots[1] = (2 * sqrt(-Q) * cos((theta + 2 * M_PI) / 3) + z);
	est_roots[2] = (2 * sqrt(-Q) * cos((theta + 4 * M_PI) / 3) + z);
	return est_roots;
}