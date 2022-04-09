#include "polynomials.h"
#include <algorithm>

///                                 ///   
///         COMPLEX NUMBERS         ///
///                                 ///   

void print(complex z)
{
    if (z.imag() == 0)
        std::cout << z.real();
    else if (z.imag() > 0)
        std::cout << z.real() << "+" << z.imag() << "i";
    else // z.imag() < 0
        std::cout << z.real() << z.imag() << "i";
}

///                             ///
///         POLYNOMIAL          ///
///                             ///

polynomial::polynomial() : coefs(0), roots(0)
{
}

polynomial::polynomial(int degree) : coefs(degree, 0), roots(degree, 0)
{
}

polynomial::~polynomial()
{
}



complex polynomial::error_est_abs(complex_vector est_roots)
{
    //if (est_roots.size()!=roots.size())
    //  exception
    //else
    complex est = 0;
    for (int i = 0; i < roots.size(); ++i)
        est += abs(est_roots[i] - roots[i]);
    return est;
}

void polynomial::info()
{
    std::cout << "~~~~~~~~~~~~~~~~" << std::endl << "Degree of the polynomial: " << roots.size();
    std::cout << std::endl << "The polynomial: x^" << coefs.size();
    for (int i = coefs.size() - 1; i > 0; --i)
    {
        std::cout << " + ("; print(coefs[i]); std::cout << ")";
        std::cout << "*x^" << i;
    }
    std::cout << " + "; print(coefs[0]);
    std::cout << std::endl << "Roots of the polynomial: ";
    for (auto v : roots)
    {
        print(v); std::cout << ", ";
    }
    std::cout << std::endl << "~~~~~~~~~~~~~~~~";
}

complex_vector polynomial::get_coefs()
{
    return coefs;
}

complex_vector polynomial::get_roots()
{
    return roots;
}

void polynomial::count_coefs()
{
}



///                                       ///
///          3 DEGREE POLYNOMIAL          ///
///                                       ///



third_degree_polynomial::third_degree_polynomial() : polynomial(3)
{
}

third_degree_polynomial::third_degree_polynomial(complex a, complex b, complex c) : polynomial(3)
{
    roots[0] = a; roots[1] = b; roots[2] = c;
    count_coefs();
}

third_degree_polynomial::~third_degree_polynomial()
{
}

void third_degree_polynomial::count_coefs()
{
    coefs[0] = -(roots[0] * roots[1] * roots[2]);
    coefs[1] = roots[0] * roots[1] + roots[0] * roots[2] + roots[1] * roots[2];
    coefs[2] = -(roots[0] + roots[1] + roots[2]);
}



///                                       ///
///          4 DEGREE POLYNOMIAL          ///
///                                       ///



fourth_degree_polynomial::fourth_degree_polynomial() : polynomial(4)
{
}

fourth_degree_polynomial::fourth_degree_polynomial(complex a, complex b, complex c, complex d) : polynomial(4)
{
    roots[0] = a; roots[1] = b; roots[2] = c; roots[3] = d;
    count_coefs();
}

fourth_degree_polynomial::~fourth_degree_polynomial()
{
}

void fourth_degree_polynomial::count_coefs()
{
    coefs[0] = roots[0] * roots[1] * roots[2] * roots[3];
    coefs[1] = -(roots[0] * roots[1] * roots[2] + roots[0] * roots[1] * roots[3] + roots[0] * roots[2] * roots[3] + roots[1] * roots[2] * roots[3]);
    coefs[2] = roots[0] * roots[1] + roots[0] * roots[2] + roots[0] * roots[3] + roots[1] * roots[2] + roots[1] * roots[3] + roots[2] * roots[3];
    coefs[3] = -(roots[0] + roots[1] + roots[2] + roots[3]);
}