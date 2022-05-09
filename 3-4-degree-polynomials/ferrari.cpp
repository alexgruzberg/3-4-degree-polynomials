// Ferrari's solution for quartic equations //
// https://math.stackexchange.com/questions/785/is-there-a-general-formula-for-solving-4th-degree-equations-quartic //
// https://en.wikipedia.org/wiki/Quartic_function#Solution_methods //

#include "polynomials.h";
#include "cardon.cpp";

/*template<typename T>
std::vector<T> ferrari(fourth_degree_polynomial<T> P)
{
	std::complex<float> half(0.5, 0); std::complex<float> quarter(0.25, 0);

	std::vector<T> coefs = P.get_coefs();
	std::complex<float> p((8 * coefs[2] - 3 * coefs[3] * coefs[3]) * 0.125, 0);
	std::complex<float> q((coefs[3] * coefs[3] * coefs[3] - 4 * coefs[3] * coefs[2] + 8 * coefs[1]) * 0.125, 0);
	std::complex<float> r((-3 * pow(coefs[3], 4) + 256 * coefs[0] - 64 * coefs[3] * coefs[1] + 16 * coefs[3] * coefs[3] * coefs[2]) * 0.00390625, 0);
	// depressed polynomial y^4 + p * y^2 + q * y + r = 0, x = y - b / 4;

	std::vector<T> roots(4);

	if (q == std::complex<float>(0, 0))
	{
		// y^4 + p * y^2 + r = 0
		std::complex<float> D1 = p * p - std::complex<float>(4,0) * r;
		roots[1] = -sqrt((-p + sqrt(D1)) * half) - coefs[3] * quarter;
		roots[2] = -sqrt((-p - sqrt(D1)) * half) - coefs[3] * quarter;
		roots[3] = sqrt((-p + sqrt(D1)) * half) - coefs[3] * quarter;
		roots[4] = sqrt((-p - sqrt(D1)) * half) - coefs[3] * quarter;

		return roots;
	}

	// resolvent cubic : 8 * m^3 + 8 * p * m^2 + (2 * p^2 - 8 * r) * m - q^2 = 0 //

	third_degree_polynomial<T> Resolvent = third_degree_polynomial<T>(1, p, p * p * 0.25 - r, q * q * -0.125);
	std::vector<T> Resolvent_roots = cardon(Resolvent);
	T m = Resolvent_roots[0];

	roots[1] = -coefs[3] * 0.25 - (sqrt(2 * m) + sqrt(-(2 * p + 2 * m - sqrt(2) * q / sqrt(m)))) * 0.5;
	roots[2] = -coefs[3] * 0.25 - (sqrt(2 * m) - sqrt(-(2 * p + 2 * m - sqrt(2) * q / sqrt(m)))) * 0.5;
	roots[3] = -coefs[3] * 0.25 + (sqrt(2 * m) - sqrt(-(2 * p + 2 * m + sqrt(2) * q / sqrt(m)))) * 0.5;
	roots[4] = -coefs[3] * 0.25 + (sqrt(2 * m) + sqrt(-(2 * p + 2 * m + sqrt(2) * q / sqrt(m)))) * 0.5;

	return roots;
}*/