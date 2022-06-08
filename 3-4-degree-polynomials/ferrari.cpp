// Ferrari's solution for quartic equations //
// https://math.stackexchange.com/questions/785/is-there-a-general-formula-for-solving-4th-degree-equations-quartic //
// https://en.wikipedia.org/wiki/Quartic_function#Solution_methods //

#include "polynomials.h";
#include "cardon.cpp";
#include <complex>
using namespace std::complex_literals;

template<typename T>
std::vector<T> ferrari(fourth_degree_polynomial<T> P)
{
	std::vector<float> coefs = P.get_coefs();

	std::complex<float> c_half(0.5, 0); std::complex<float> c_two(2, 0);
	
	float b = coefs[3], c = coefs[2], d = coefs[1], e = coefs[0];

	std::complex<float> p((8 * c - 3 * b * b) * 0.125, 0);
	std::complex<float> q((b * (b * b - 4 * c) + 8 * d) * 0.125, 0);
	std::complex<float> r((b * (-3 * b * b * b - 16 * (4 * d - b * c)) + 256 * e) * 0.00390625, 0);
	// depressed polynomial y^4 + p * y^2 + q * y + r = 0, x = y - b / 4;

	std::vector<T> roots(4);

	std::complex<float> fourth_coef_three = -b * c_half * c_half;

	if (q.real() == 0)
	{
		// y^4 + p * y^2 + r = 0
		std::complex<float> D1 = p * p - c_two * c_two * r;
		std::complex<float> sqrt_d_one = sqrt(D1);
		roots[0] = (-sqrt((-p + sqrt(sqrt_d_one)) * c_half) - fourth_coef_three).real();
		roots[1] = (-sqrt((-p - sqrt(sqrt_d_one)) * c_half) - fourth_coef_three).real();
		roots[2] = (sqrt((-p + sqrt(sqrt_d_one)) * c_half) - fourth_coef_three).real();
		roots[3] = (sqrt((-p - sqrt(sqrt_d_one)) * c_half) - fourth_coef_three).real();

		return roots;
	}

	// resolvent cubic : 8 * m^3 + 8 * p * m^2 + (2 * p^2 - 8 * r) * m - q^2 = 0 //

	third_degree_polynomial<std::complex<float>> Resolvent(p.real(), p.real() * p.real() * 0.25 - r.real(), q.real() * q.real() * -0.125);
	std::vector<std::complex<float>> Resolvent_roots = cardon(Resolvent);
	std::complex<float> m = Resolvent_roots[0];
	std::complex<float> sqrt_two_m = sqrt(c_two * m); std::complex < float> two_p_two_m = c_two * (p + m);
	std::complex<float> sqrt_two_q_sqrt_m = sqrt(c_two) * q / sqrt(m);

	roots[0] = (fourth_coef_three - (sqrt_two_m + sqrt(-(two_p_two_m - sqrt_two_q_sqrt_m))) * c_half).real();
	roots[1] = (fourth_coef_three - (sqrt_two_m - sqrt(-(two_p_two_m - sqrt_two_q_sqrt_m))) * c_half).real();
	roots[2] = (fourth_coef_three + (sqrt_two_m - sqrt(-(two_p_two_m + sqrt_two_q_sqrt_m))) * c_half).real();
	roots[3] = (fourth_coef_three + (sqrt_two_m + sqrt(-(two_p_two_m + sqrt_two_q_sqrt_m))) * c_half).real();

	sort(roots.begin(), roots.end());

	return roots;
}

template<typename T>
std::vector<std::complex<T>> ferrari(fourth_degree_polynomial<std::complex<T>> P)
{
	std::vector<float> coefs = P.get_coefs();

	std::complex<float> c_half(0.5, 0); std::complex<float> c_two(2, 0);

	float b = coefs[3], c = coefs[2], d = coefs[1], e = coefs[0];

	std::complex<float> p((8 * c - 3 * b * b) * 0.125, 0);
	std::complex<float> q((b * (b * b - 4 * c) + 8 * d) * 0.125, 0);
	std::complex<float> r((b * (-3 * b * b * b - 16 * (4 * d - b * c)) + 256 * e) * 0.00390625, 0);
	// depressed polynomial y^4 + p * y^2 + q * y + r = 0, x = y - b / 4;

	std::vector<std::complex<T>> roots(4);

	std::complex<float> fourth_coef_three = -b * c_half * c_half;

	if (q.real() == 0)
	{
		// y^4 + p * y^2 + r = 0
		std::complex<float> D1 = p * p - c_two * c_two * r;
		std::complex<float> sqrt_d_one = sqrt(D1);
		roots[0] = -sqrt((-p + sqrt_d_one) * c_half) - fourth_coef_three;
		roots[1] = -sqrt((-p - sqrt_d_one) * c_half) - fourth_coef_three;
		roots[2] = sqrt((-p + sqrt_d_one) * c_half) - fourth_coef_three;
		roots[3] = sqrt((-p - sqrt_d_one) * c_half) - fourth_coef_three;

		return roots;
	}

	// resolvent cubic : 8 * m^3 + 8 * p * m^2 + (2 * p^2 - 8 * r) * m - q^2 = 0 //

	third_degree_polynomial<std::complex<float>> Resolvent(p.real(), p.real() * p.real() * 0.25 - r.real(), q.real() * q.real() * -0.125);
	std::vector<std::complex<float>> Resolvent_roots = cardon(Resolvent);
	std::complex<float> m = Resolvent_roots[0];
	std::complex<float> sqrt_two_m = sqrt(c_two * m); std::complex<float> two_p_two_m = c_two * (p + m);
	std::complex<float> sqrt_two_q_sqrt_m = sqrt(c_two) * q / sqrt(m);

	roots[0] = fourth_coef_three - (sqrt_two_m + sqrt(-(two_p_two_m - sqrt_two_q_sqrt_m))) * c_half;
	roots[1] = fourth_coef_three - (sqrt_two_m - sqrt(-(two_p_two_m - sqrt_two_q_sqrt_m))) * c_half;
	roots[2] = fourth_coef_three + (sqrt_two_m - sqrt(-(two_p_two_m + sqrt_two_q_sqrt_m))) * c_half;
	roots[3] = fourth_coef_three + (sqrt_two_m + sqrt(-(two_p_two_m + sqrt_two_q_sqrt_m))) * c_half;


	sort(roots.begin(), roots.end(), [](std::complex<T> a, std::complex<T> b) {return (abs(a) > abs(b)); });
	for (int i = 0; i < 2; ++i)
	{
		if (roots[2*i].imag() < -0.000001)
			swap(roots[2 * i], roots[2 * i + 1]);
	}

	return roots;
}