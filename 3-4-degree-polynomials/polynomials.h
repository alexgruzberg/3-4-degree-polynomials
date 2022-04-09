#pragma once

#include <vector>
#include <iostream>
#include "math.h"
#include <complex>

typedef std::vector<std::complex<float>> complex_vector;
typedef std::complex<float> complex;

class polynomial //abstract class for polynomial
{
public:
    polynomial();
    polynomial(int degree);
    virtual ~polynomial() = 0;

    void info();
    complex_vector get_coefs();
    complex_vector get_roots();
protected:
    virtual void count_coefs();
    complex error_est_abs(complex_vector est_roots);  //error estimation
    complex_vector roots;
    complex_vector coefs;   //coefficients of a polynomial
};

class third_degree_polynomial : public polynomial
    ///
    ///         x^3 + coefs[2]*x^2 + coefs[1] * x + coefs[0]
    ///
{
public:
    third_degree_polynomial();
    third_degree_polynomial(complex a, complex b, complex c); //constructor that receives roots
    ~third_degree_polynomial();
protected:
    void count_coefs();
};

class fourth_degree_polynomial : public polynomial
    ///
    ///         x^4 + coefs[3]*x^3 + coefs[2]*x^2 + coefs[1] * x + coefs[0]
    ///
{
public:
    fourth_degree_polynomial();
    fourth_degree_polynomial(complex a, complex b, complex c, complex d); //constructor that receives roots
    ~fourth_degree_polynomial();
protected:
    void count_coefs();
};