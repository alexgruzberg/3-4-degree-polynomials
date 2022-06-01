#pragma once
#define _USE_MATH_DEFINES

#include "exceptions.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <complex>

template<typename T>
class polynomial //abstract class for polynomial
{
public:
    polynomial();
    polynomial(int degree);
    virtual ~polynomial() = 0;

    void info();    //useful info about polynomial
    std::vector<float> get_coefs();
    std::vector<T> get_roots();
    float error_est_sum(std::vector<T> est_roots);  //error estimation
    float error_est_max(std::vector<T> est_roots);
protected:
    virtual void count_coefs(); //counts coefficients of the polynomial using received roots
    std::vector<float> coefs;   //coefficients of a polynomial
    std::vector<T> roots;
};

template <typename T>
polynomial<T>::polynomial() : coefs(0), roots(0) {}

template <typename T>
polynomial<T>::polynomial(int degree) : coefs(degree, 0), roots(degree, 0) {}

template <typename T>
polynomial<T>::~polynomial() {}

template <typename T>
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

template <typename T>
std::vector<float> polynomial<T>::get_coefs()
{
    return this->coefs;
}

template <typename T>
std::vector<T> polynomial<T>::get_roots()
{
    return this->roots;
}

template <typename T>
float polynomial<T>::error_est_sum(std::vector<T> est_roots)
{
    if (est_roots.size() != roots.size())
        throw unexcpected_number_of_roots();
    float est = 0;
    for (int i = 0; i < this->roots.size(); ++i)
        est += abs(est_roots[i] - this->roots[i]);
    return est;
}

template <typename T>
float polynomial<T>::error_est_max(std::vector<T> est_roots)
{
    if (est_roots.size() != roots.size())
        throw unexcpected_number_of_roots();
    float est = 0;
    for (int i = 0; i < this->roots.size(); ++i)
        est += abs(est_roots[i] - this->roots[i]) / std::max(abs(est_roots[i]), abs(this->roots[i]));
    return est;
}

template <typename T>
void polynomial<T>::count_coefs() {}




//----------------------------------------------------------------------------------------------------------//




template<typename T>
class third_degree_polynomial : public polynomial<T>        //third degree polynomial that has no complex roots
    ///
    ///         x^3 + coefs[2]*x^2 + coefs[1] * x + coefs[0]
    ///
{
public:
    third_degree_polynomial();
    third_degree_polynomial(float a, float b, float c); //constructor that receives coefficients

    third_degree_polynomial(std::vector<T> rec_roots); //constructor that receives roots

    ~third_degree_polynomial();
protected:
    void count_coefs();
};

template <typename T>
third_degree_polynomial<T>::third_degree_polynomial() : polynomial<T>(3) {}

template <typename T>
third_degree_polynomial<T>::third_degree_polynomial(float a, float b, float c) : polynomial<T>(3) //constructor that receives coefficients
{
    this->coefs[2] = a; this->coefs[1] = b; this->coefs[0] = c;
}

template <typename T>
third_degree_polynomial<T>::third_degree_polynomial(std::vector<T> rec_roots) : polynomial<T>(3) //constructor that receives roots
{
    for (int i = 0; i < 3; ++i) this->roots[i] = rec_roots[i];
    std::sort(this->roots.begin(), this->roots.end());
    count_coefs();
}

template <typename T>
third_degree_polynomial<T>::~third_degree_polynomial() {}

template <typename T>
void third_degree_polynomial<T>::count_coefs()
{
    this->coefs[0] = -(this->roots[0] * this->roots[1] * this->roots[2]);
    this->coefs[1] = this->roots[0] * this->roots[1] + this->roots[0] * this->roots[2] + this->roots[1] * this->roots[2];
    this->coefs[2] = -(this->roots[0] + this->roots[1] + this->roots[2]);
}

template<typename T>
class third_degree_polynomial<std::complex<T>> : public polynomial<std::complex<T>>     //third degree polynomial that can have complex roots
    ///
    ///         x^3 + coefs[2]*x^2 + coefs[1] * x + coefs[0]
    ///
{
public:
    third_degree_polynomial();
    third_degree_polynomial(float a, float b, float c); //constructor that receives coefficients

    third_degree_polynomial(std::vector<std::complex<T>> rec_roots); //constructor that receives roots

    ~third_degree_polynomial();
protected:
    void count_coefs();
};

template <typename T>
third_degree_polynomial<std::complex<T>>::third_degree_polynomial() : polynomial<std::complex<T>>(3) {}

template <typename T>
third_degree_polynomial<std::complex<T>>::third_degree_polynomial(float a, float b, float c) : polynomial<std::complex<T>>(3) //constructor that receives coefficients
{
    this->coefs[2] = a; this->coefs[1] = b; this->coefs[0] = c;
}

template <typename T>
third_degree_polynomial<std::complex<T>>::third_degree_polynomial(std::vector<std::complex<T >> rec_roots) : polynomial<std::complex<T>>(3) //constructor that receives roots
{
    for (int i = 0; i < 3; ++i) this->roots[i] = rec_roots[i];
    count_coefs();
}

template <typename T>
third_degree_polynomial<std::complex<T>>::~third_degree_polynomial() {}

template <typename T>
void third_degree_polynomial<std::complex<T>>::count_coefs()
{
     this->coefs[0] = (-(this->roots[0] * this->roots[1] * this->roots[2])).real();
     this->coefs[1] = (this->roots[0] * this->roots[1] + this->roots[0] * this->roots[2] + this->roots[1] * this->roots[2]).real();
     this->coefs[2] = ((-(this->roots[0] + this->roots[1] + this->roots[2]))).real();
}





//----------------------------------------------------------------------------------------------------------//




template<typename T>
class fourth_degree_polynomial : public polynomial<T>       //fourth degree polynomial without complex roots
    ///
    ///         x^4 + coefs[3]*x^3 + coefs[2]*x^2 + coefs[1] * x + coefs[0]
    ///
{
public:
    fourth_degree_polynomial();
    fourth_degree_polynomial(float a, float b, float c, float d); //constructor that receives coefficients
    fourth_degree_polynomial(std::vector<T> rec_roots); //constructor that receives roots
    ~fourth_degree_polynomial();
protected:
    void count_coefs();
};

template <typename T>
fourth_degree_polynomial<T>::fourth_degree_polynomial() : polynomial<T>(4) {}

template<typename T>
fourth_degree_polynomial<T>::fourth_degree_polynomial(float a, float b, float c, float d) : polynomial<T>(4)
{
    this->coefs[3] = a; this->coefs[2] = b; this->coefs[1] = c; this->coefs[0] = d;
}

template <typename T>
fourth_degree_polynomial<T>::fourth_degree_polynomial(std::vector<T> rec_roots) : polynomial<T>(4) //constructor that receives roots
{
    for (int i = 0; i < 4; ++i) this->roots[i] = rec_roots[i];
    std::sort(this->roots.begin(), this->roots.end());
    count_coefs();
}

template <typename T>
fourth_degree_polynomial<T>::~fourth_degree_polynomial() {}

template <typename T>
void fourth_degree_polynomial<T>::count_coefs()
{
    this->coefs[0] = this->roots[0] * this->roots[1] * this->roots[2] * this->roots[3];
    this->coefs[1] = -(this->roots[0] * this->roots[1] * this->roots[2] + this->roots[0] * this->roots[1] * this->roots[3] + this->roots[0] * this->roots[2] *
        this->roots[3] + this->roots[1] * this->roots[2] * this->roots[3]);
    this->coefs[2] = this->roots[0] * this->roots[1] + this->roots[0] * this->roots[2] + this->roots[0] * this->roots[3] + this->roots[1] * this->roots[2] +
        this->roots[1] * this->roots[3] + this->roots[2] * this->roots[3];
    this->coefs[3] = -(this->roots[0] + this->roots[1] + this->roots[2] + this->roots[3]);
}



template<typename T>
class fourth_degree_polynomial<std::complex<T>> : public polynomial<std::complex<T>>       //fourth degree polynomial with complex roots
    ///
    ///         x^4 + coefs[3]*x^3 + coefs[2]*x^2 + coefs[1] * x + coefs[0]
    ///
{
public:
    fourth_degree_polynomial();
    fourth_degree_polynomial(float a, float b, float c, float d); //constructor that receives coefficients
    fourth_degree_polynomial(std::vector<std::complex<T>> rec_roots); //constructor that receives roots
    ~fourth_degree_polynomial();
protected:
    void count_coefs();
};

template <typename T>
fourth_degree_polynomial<std::complex<T>>::fourth_degree_polynomial() : polynomial<std::complex<T>>(4) {}

template<typename T>
fourth_degree_polynomial<std::complex<T>>::fourth_degree_polynomial(float a, float b, float c, float d) : polynomial<std::complex<T>>(4)
{
    this->coefs[3] = a; this->coefs[2] = b; this->coefs[1] = c; this->coefs[0] = d;
}

template <typename T>
fourth_degree_polynomial<std::complex<T>>::fourth_degree_polynomial(std::vector<std::complex<T>> rec_roots) : polynomial<std::complex<T>>(4) //constructor that receives roots
{
    for (int i = 0; i < 4; ++i) this->roots[i] = rec_roots[i];
    std::sort(this->roots.begin(), this->roots.end());
    count_coefs();
}

template <typename T>
fourth_degree_polynomial<std::complex<T>>::~fourth_degree_polynomial() {}

template <typename T>
void fourth_degree_polynomial<std::complex<T>>::count_coefs()
{
    this->coefs[0] = (this->roots[0] * this->roots[1] * this->roots[2] * this->roots[3]).real();
    this->coefs[1] = (-(this->roots[0] * this->roots[1] * this->roots[2] + this->roots[0] * this->roots[1] * this->roots[3] + this->roots[0] * this->roots[2] *
        this->roots[3] + this->roots[1] * this->roots[2] * this->roots[3])).real();
    this->coefs[2] = (this->roots[0] * this->roots[1] + this->roots[0] * this->roots[2] + this->roots[0] * this->roots[3] + this->roots[1] * this->roots[2] +
        this->roots[1] * this->roots[3] + this->roots[2] * this->roots[3]).real();
    this->coefs[3] = (-(this->roots[0] + this->roots[1] + this->roots[2] + this->roots[3])).real();
}