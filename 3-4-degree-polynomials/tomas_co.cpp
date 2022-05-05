//Tomas Co 10/07/2014

#include "polynomials.h"
//	x^3 + a x^2 + b x + c

template<typename T>
std::vector<T> third_degree_polynomial<T>::tomas_co()
{
	float p = (3 * this->coefs[1] - pow(this->coefs[2], 2))/3;
	float q = (2 * pow(this->coefs[2], 3) - 9 * this->coefs[2] * this->coefs[1] + 27 * this->coefs[0]) / 27;
	std::vector<T> est_roots(3);
	float A = 2 * sqrt(-p / 3); //p is less than 0
	float phi = acos((3 * q) / (A * p));
	float B = -this->coefs[2] / 3;
	for (int i = 0; i < 3; ++i)
		est_roots[i] = A * cos((phi + 2*i*M_PI)/3) + B;
	
	return est_roots;
}