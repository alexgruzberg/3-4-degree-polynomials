/*
Здесь представлена реализация метода [преобразования] Чирнгауза для нахождения 
корней кубического уравнения, записанного в виде x^3 + bx^2 + cx + d = 0.

Программа адаптирована под многочлены, содержащие комплексные коэффициенты,
поэтому класс third_degree_polynomial не используется ввиду отсутствия возможности работать с комплексными коэффициентами.

Даны краткие комментарии к функциям (полное описание см. в документации).
*/

#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>
#include <vector>
#include <ctime>

#include "polynomials.h"
#include "exceptions.h"

using namespace std;
//функция для нахождения всех [трёх] значений кубического корня
template<typename T>
vector<complex<T>> cubic_roots(complex<T> z)
{
    complex<double> one_third = { 1.0 / 3.0, 0 };

    complex<double> phi = arg(z);   //аргумент комплексного числа
    complex<double> module = abs(z);    //модуль комплексного числа

    complex<double> cubic_root_module = pow(module, one_third);  //кубический корень из модуля комплексного числа

    vector<complex<T>> roots;
    for (int n = 0; n <= 2; n++)    //вычисление значений корня
    {
        double p = 2 * M_PI * n;

        double alfa = real((phi + complex<double>{(p), 0})* one_third);
        complex<double> value = cubic_root_module * complex<double>{ cos(alfa), sin(alfa) };
        roots.push_back(value);
    }
    return roots;
}

//функция для нахождения корней квадратного уравнения со старшим коэффициентом a = 1
template<typename T>
vector<complex<T>> roots_of_square_poly(vector<complex<T>> s)
{
    complex<double> one_second = { 0.5, 0 };
    complex<double> D = s[0] * s[0] - complex<double>{ 4.0, 0 } * s[1]; //дискриминант
    complex<double> rootD = sqrt(D);    //корень из дискриминанта
    complex<T> x_1 = (-s[0] + rootD) * one_second;
    complex<T> x_2 = (-s[0] - rootD) * one_second;
    return { x_1, x_2 };
}

//функция для приведения общего вида многочлена 3-ей степени [x^3 + bx^2 + cx + d] к виду t^3 + pt + q
template<typename T>
vector<complex<T>> canonical_reduction(vector<complex<T>> P)
{
    complex<double> one_third = { 1.0 / 3.0, 0 };
    complex<double> one_twentyseventh = { 1.0 / 27.0, 0 };
   
    complex<double> d = P[2], c = P[1], b = P[0];   //коэффициенты

    complex<T> p = (complex<double>{ 3, 0 } * c - b * b)* one_third;    //новые коэффициенты
    complex<T> q = (complex<double>{27, 0} * d - complex<double>{9, 0} * b * c + complex<double>{2, 0} * b * b * b) * one_twentyseventh;

    return { p,q };
}

//Преобразование Чирнгауза
template<typename T>
vector<complex<T>> Tschirnhaus_transformation(vector<complex<T>> F)
{
    if ((F[1] == complex<double>{0, 0}) && (F[2] == complex<double>{0, 0}))     //быстрый подсчёт корней при c = d = 0
    {
        return { complex<T>{0, 0},complex<T>{0, 0},-F[0] };
    }
    complex<double> one_third = { 1.0 / 3.0, 0 };
    complex<double> two_third = { 2.0 / 3.0, 0 };
    complex<double> one_twentyseventh = { 1.0 / 27.0, 0 };

    vector<complex<double>> reduct = canonical_reduction(F);    //приведение к упрощённому виду [t^3 + pt + q]
    complex<double> p =  reduct[0];
    complex<double> q = reduct[1];

    complex<double> qq = q * q;
    complex<double> qqq = q * q * q;
    complex<double> pp = p * p;
    complex<double> ppp = p * p * p;
    complex<double> reverse_p = complex<double>{ 1,0 } / p;
    if ((isinf(real(reverse_p)) || isnan(real(reverse_p))) && (isinf(imag(reverse_p)) || isnan(imag(reverse_p))))   //проверка на ноль в знаменателе [если да - поиск корней для t^3 + q = 0]
    {
        //throw division_by_zero();
        complex<double> y_cubic_0 = -q;
        vector<complex<T>> x0 = cubic_roots(y_cubic_0);
        for (int l = 0; l <= 2; l++)
        {
            x0[l] = x0[l] - F[0] * one_third;
        }
        return x0;
    }
    complex<double> reverse_pp = reverse_p* reverse_p;

    vector<complex<double>> find_a = { complex<double>{ 3, 0 }*q* reverse_p, -one_third * p };  //коэффициенты для квадратного уравнения
    vector<complex<double>> a = roots_of_square_poly(find_a);   //корни квадратного уравнения
    complex<double> b = two_third * p;
    /*Далее в вычислениях участвует только один корень квадратного уравнения (пояснение см. в документации)*/
    complex<double> y_cubic_1 = complex<double>{ 9, 0 } * qqq * a[0] * reverse_pp + complex<double>{ 4, 0 }*one_third* a[0] * p * q - complex<double>{ 8, 0 }*one_twentyseventh* ppp - complex<double>{ 2, 0 } * qq;

    vector<complex<double>> y1 = cubic_roots(y_cubic_1);    //проверка кубических корней из y1 на нулевое значение
    if (y1 == vector<complex<double>>{0, 0, 0})
    {
        vector<complex<double>> sq_poly_coefs0 = { a[0], b };
        vector<complex<T>> main_roots0 = roots_of_square_poly(sq_poly_coefs0);
        main_roots0.push_back(a[0]);
        return main_roots0;
    }
    
    double eps = 1e-3;
    vector<complex<double>> six_roots;  
    for (int i = 0; i <= 2; i++)    //вычисление шести корней
    {
        complex<double> sum_b_y1 = b + y1[i];
        vector<complex<double>> sq_poly_coefs1 = { a[0], sum_b_y1 };
        vector<complex<double>> main_roots1 = roots_of_square_poly(sq_poly_coefs1);
        six_roots.push_back(main_roots1[0] - F[0] * one_third);
        six_roots.push_back(main_roots1[1] - F[0] * one_third);
    }
    vector<complex<T>> final_roots;
    complex<double> d0 = F[2], c0 = F[1], b0 = F[0];
    for (int j = 0; j <= 5; j++)    //отбрасывание посторонних корней
    {
        complex<double> value = six_roots[j];
        complex<double> x_square = value * value;
        complex<double> x_cubic = value * x_square;
        double re_root = abs(real(x_cubic + b0 * x_square + c0 * value + d0));
        double im_root = abs(imag(x_cubic + b0 * x_square + c0 * value + d0));
        if (re_root <= eps && im_root <= eps)
        {
            final_roots.push_back(value);
        }
    }
    return final_roots;
}

template<typename T>
void print(vector<complex<T>> v)    //вывод в консоль
{
    for (unsigned int k = 0; k < v.size(); k++)
    {
        cout << "x_" << k+1 << " = " << real(v[k]) << " + " << imag(v[k]) << "*i" << "\n";
    }
}
