#pragma once
#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>

template<typename T>
class polynomial //abstract class for polynomial
{
public:
    polynomial();
    polynomial(int degree);
    virtual ~polynomial() = 0;

    void info();
    std::vector<T> get_coefs();
    std::vector<T> get_roots();
    float error_est_sum(std::vector<T> est_roots);  //error estimation
    float error_est_max(std::vector<T> est_roots);
protected:
    virtual void count_coefs();
    std::vector<T> roots;
    std::vector<float> coefs;   //coefficients of a polynomial
};

template<typename T>
class third_degree_polynomial : public polynomial<T>
    ///
    ///         x^3 + coefs[2]*x^2 + coefs[1] * x + coefs[0]
    ///
{
public:
    third_degree_polynomial();
    third_degree_polynomial(float a, float b, float c); //constructor that receives coefficients
    third_degree_polynomial(std::vector<T> rec_roots); //constructor that receives roots
    ~third_degree_polynomial();

    std::vector<T> tiruneh();
    std::vector<T> cardon();
    std::vector<T> tomas_co();
protected:
    void count_coefs();
};

template<typename T>
class fourth_degree_polynomial : public polynomial<T>
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

#include "polynomials.cpp"