//A simplified expression for the solution of cubic polynomial equations using function evaluation-Tiruneh-2020//
#include "polynomials.h"

//	x^3 + a x^2 + b x + c

std::vector<float> tiruneh(third_degree_polynomial p)
{
	std::vector<float> coefs = p.get_coefs();
	float z = -coefs[2]/3; //z = -b/3
	float Q = (-pow(coefs[2],2) / 3 + coefs[2]) / 3;
	float R = -(pow(z, 3) + coefs[2] * pow(z, 2) + coefs[1] * z + coefs[0]) / 2;
	
	std::vector<float> roots;
	return roots;
}