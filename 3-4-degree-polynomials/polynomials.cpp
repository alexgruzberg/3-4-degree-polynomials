#pragma once
#include "polynomials.h"    

///                             ///
///         POLYNOMIAL          ///
///                             ///

template<typename T>
polynomial<T>::polynomial() : coefs(0), roots(0)
{
}

template<typename T>
polynomial<T>::polynomial(int degree) : coefs(degree, 0), roots(degree, 0)
{
}

template<typename T>
polynomial<T>::~polynomial()
{
}



template<typename T>
float polynomial<T>::error_est_sum(std::vector<T> est_roots)
{
    //if (est_roots.size()!=roots.size())
    //  exception
    //else
    float est = 0;
    for (int i = 0; i < this->roots.size(); ++i)
        est += abs(est_roots[i] - this->roots[i]);
    return est;
}

template<typename T>
float polynomial<T>::error_est_max(std::vector<T> est_roots)
{
    //if (est_roots.size()!=roots.size())
    //  exception
    //else
    float est = 0;
    for (int i = 0; i < this->roots.size(); ++i)
        est += abs(est_roots[i] - this->roots[i]) / std::max(est_roots[i], this->roots[i]);
    return est;
}

template<typename T>
void polynomial<T>::info()
{
    std::cout << "~~~~~~~~~~~~~~~~" << std::endl << "Degree of the polynomial: " << this->roots.size();
    std::cout << std::endl << "The polynomial: x^" << this->coefs.size();
    for (int i = this->coefs.size() - 1; i > 0; --i)
    {
        if (coefs[i] < 0) std::cout << " - " << std::abs(this->coefs[i]);
        else std::cout << " + " << this->coefs[i];
        std::cout << "*x^" << i;
    }
    if (this->coefs[0] < 0) std::cout << " - " << std::abs(this->coefs[0]);
    else std::cout << " + " << this->coefs[0];
    std::cout << std::endl << "Roots of the polynomial: ";
    for (auto v : this->roots)
    {
        std::cout << v << ", ";
    }
    std::cout << std::endl << "~~~~~~~~~~~~~~~~" << std::endl;
}

template<typename T>
std::vector<T> polynomial<T>::get_coefs()
{
    return this->coefs;
}

template<typename T>
std::vector<T> polynomial<T>::get_roots()
{
    return roots;
}

template<typename T>
void polynomial<T>::count_coefs()
{
}



///                                       ///
///          3 DEGREE POLYNOMIAL          ///
///                                       ///



template<typename T>
third_degree_polynomial<T>::third_degree_polynomial() : polynomial<T>(3)
{
}

template<typename T>
third_degree_polynomial<T>::third_degree_polynomial(float a, float b, float c) : polynomial<T>(3)
{
    this->coefs[2] = a; this->coefs[1] = b; this->coefs[0] = c;
}

template<typename T>
third_degree_polynomial<T>::third_degree_polynomial(std::vector<T> rec_roots) : polynomial<T>(3)
{
    for (int i = 0; i < 3; ++i) this->roots[i] = rec_roots[i];
    std::sort(this->roots.begin(), this->roots.end());
    count_coefs();
}

template<typename T>
third_degree_polynomial<T>::~third_degree_polynomial()
{
}

template<typename T>
void third_degree_polynomial<T>::count_coefs()
{
    this->coefs[0] = -(this->roots[0] * this->roots[1] * this->roots[2]);
    this->coefs[1] = this->roots[0] * this->roots[1] + this->roots[0] * this->roots[2] + this->roots[1] * this->roots[2];
    this->coefs[2] = -(this->roots[0] + this->roots[1] + this->roots[2]);
}



///                                       ///
///          4 DEGREE POLYNOMIAL          ///
///                                       ///



template<typename T>
fourth_degree_polynomial<T>::fourth_degree_polynomial() : polynomial<T>(4)
{
}

template<typename T>
fourth_degree_polynomial<T>::fourth_degree_polynomial(float a, float b, float c, float d) : polynomial<T>(4)
{
    this->coefs[3] = a; this->coefs[2] = b; this->coefs[1] = c; this->coefs[0] = d;
}

template<typename T>
fourth_degree_polynomial<T>::fourth_degree_polynomial(std::vector<T> rec_roots) : polynomial<T>(4)
{
    for (int i = 0; i < 4; ++i) this->roots[i] = rec_roots[i];
    std::sort(this->roots.begin(), this->roots.end());
    count_coefs();
}

template<typename T>
fourth_degree_polynomial<T>::~fourth_degree_polynomial()
{
}

template<typename T>
void fourth_degree_polynomial<T>::count_coefs()
{
    this->coefs[0] = this->roots[0] * this->roots[1] * this->roots[2] * this->roots[3];
    this->coefs[1] = -(this->roots[0] * this->roots[1] * this->roots[2] + this->roots[0] * this->roots[1] * this->roots[3] + this->roots[0] * this->roots[2] * 
                        this->roots[3] + this->roots[1] * this->roots[2] * this->roots[3]);
    this->coefs[2] = this->roots[0] * this->roots[1] + this->roots[0] * this->roots[2] + this->roots[0] * this->roots[3] + this->roots[1] * this->roots[2] + 
                        this->roots[1] * this->roots[3] + this->roots[2] * this->roots[3];
    this->coefs[3] = -(this->roots[0] + this->roots[1] + this->roots[2] + this->roots[3]);
}
