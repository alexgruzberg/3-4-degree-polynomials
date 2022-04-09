#pragma once

#include <vector>
#include <iostream>
#include "math.h"

class polynomial //abstract class for polynomial
{
public:
    polynomial();
    polynomial(int degree);
    virtual ~polynomial() = 0;

    void info();
    std::vector<float> get_coefs();
    std::vector<float> get_roots();
protected:
    virtual void count_coefs();
    float error_est_sum(std::vector<float> est_roots);  //error estimation
    float error_est_max(std::vector<float> est_roots);
    std::vector<float> roots;
    std::vector<float> coefs;   //coefficients of a polynomial
};

class third_degree_polynomial : public polynomial
    ///
    ///         x^3 + coefs[2]*x^2 + coefs[1] * x + coefs[0]
    ///
{
public:
    third_degree_polynomial();
    third_degree_polynomial(float a, float b, float c); //constructor that receives roots
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
    fourth_degree_polynomial(float a, float b, float c, float d); //constructor that receives roots
    ~fourth_degree_polynomial();
protected:
    void count_coefs();
};