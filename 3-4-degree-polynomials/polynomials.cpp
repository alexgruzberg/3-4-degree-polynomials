#include "polynomials.h"

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



float polynomial::error_est_sum(std::vector<float> est_roots)
{
    //if (est_roots.size()!=roots.size())
    //  exception
    //else
    float est = 0;
    for (int i = 0; i < roots.size(); ++i)
        est += abs(est_roots[i] - roots[i]);
    return est;
}

float polynomial::error_est_max(std::vector<float> est_roots)
{
    //if (est_roots.size()!=roots.size())
    //  exception
    //else
    float est = 0;
    for (int i = 0; i < roots.size(); ++i)
        est += abs(est_roots[i] - roots[i]) / std::max(est_roots[i], roots[i]);
    return est;
}

void polynomial::info()
{
    std::cout << "~~~~~~~~~~~~~~~~" << std::endl << "Degree of the polynomial: " << roots.size();
    std::cout << std::endl << "The polynomial: x^" << coefs.size();
    for (int i = coefs.size() - 1; i > 0; --i)
    {
        if (coefs[i] < 0) std::cout << " - " << std::abs(coefs[i]);
        else std::cout << " + " << coefs[i];
        std::cout << "*x^" << i;
    }
    if (coefs[0] < 0) std::cout << " - " << std::abs(coefs[0]);
    else std::cout << " + " << coefs[0];
    std::cout << std::endl << "Roots of the polynomial: ";
    for (auto v : roots)
    {
        std::cout << v << ", ";
    }
    std::cout << std::endl << "~~~~~~~~~~~~~~~~" << std::endl;
}

std::vector<float> polynomial::get_coefs()
{
    return coefs;
}

std::vector<float> polynomial::get_roots()
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

third_degree_polynomial::third_degree_polynomial(float a, float b, float c) : polynomial(3)
{
    roots[0] = a; roots[1] = b; roots[2] = c;
    count_coefs();
}

third_degree_polynomial::third_degree_polynomial(std::vector<float> rec_roots) : polynomial(3)
{
    for (int i = 0; i < 3; ++i) roots[i] = rec_roots[i];
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

fourth_degree_polynomial::fourth_degree_polynomial(float a, float b, float c, float d) : polynomial(4)
{
    roots[0] = a; roots[1] = b; roots[2] = c; roots[3] = d;
    count_coefs();
}

fourth_degree_polynomial::fourth_degree_polynomial(std::vector<float> rec_roots) : polynomial(4)
{
    for (int i = 0; i < 4; ++i) roots[i] = rec_roots[i];
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