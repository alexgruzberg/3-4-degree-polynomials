#include "polynomials.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <math.h>
#include <corecrt_math_defines.h>
#include <complex>
using namespace std;

const float one_third = 1.0f / 3.0f; /*Константа, обозначаюшая число 1\3 */

inline float algorith(const float b, const float c, const float d)/*Функция алгоритм, возвращает число float, принимает три константы float, b это коэфициент при x^2,c при x и d - это свободный член */
{/*Функция algorith нужна для вычисления вещественнеого корня в уравнении 3 степени
 Если их несколько, берется тот, который или меньше 0, или 0, или последний*/

    float e = b * one_third;
    float f = c - b * e;/*От значения данной переменной зависит сколько будет вещественных корней и формула их вычисления*/
    float g = e * (c - 2 * e * e) - d;
    /*Переменные e, g, h, j, нужны для вычисления комплесных корней в разных случаях */
    if (f == 0)
    {
        return cbrt(g) - e;/*Возвращение единственного вещественного корня*/
    }
    else
    {
        float h = sqrt(4 * abs(f) * one_third);
        float i = 4 * g / (h * h * h);
        if (f > 0)
        {
            return h * sinh(one_third * asinh(i)) - e;/*Возвращение единственного вещественного корня*/
        }
        else
        {
            if (abs(i) > 1)
            {
                return (h * signbit(i) * cosh(one_third * acosh(abs(i))) - e);/*Возвращение единственного вещественного корня*/
            }
            else/*Если переменная f меньше нуля и модель переменной i меньше 1, то у данного многочлена 3 вещественных корня*/
            {

                float j = acos(i) * one_third;
                if ((h * cos(j) - e) <= 0)
                    return  h * cos(j) - e;
                else if ((h * cos(((float((2 * M_PI) * one_third)) + j)) - e) <= 0)
                    return h * cos(((float((2 * M_PI) * one_third)) + j)) - e;
                else
                    return h * cos(((float((2 * M_PI) * one_third)) - j)) - e;

            }
        }

    }
}

void A_Salzer(const float z, const float x, const float y, const float v, vector<complex<float>>& solution)/*Функция А_Салзер ничего не возвращает, принимает 4 константы, обозначающие коэфициенты многочлена, и вектор solution, куда будут вносится ответы*/
{
    //z-коэфициент при x^3
    //x-коэфициент при x^2
    //y-коэфициент при x
    //v-свободный член
    float z2 = z * z;//Переменная, введенная чтобы не вычислять данное выражение множество раз
    float b = -x, c = z * y - 4 * v,
        d = v * (4 * x - z2) - y * y;
    /*Переменные b, c, d это коэфициенты нового многочлена 3 степени, вещественный корень которого нужно найти*/
    float x1 = algorith(b, c, d);/*Корень многочлена 3 степени*/
    complex<float> m = 0.25f * z2 - x + x1;
    m = sqrt(m);
    complex<float> n;
    if (real(m) == 0 && imag(m) == 0)/* Если переменная m равна нулю, то n вычисляется по другой формуле*/
    {
        n = sqrt(0.25f * x1 * x1 - v);
    }
    else
    {
        n = (z * x1 - 2 * y) * 0.25f / m;
    }
    if (imag(m) == 0)/* Если переменная m вещественная, то корни многочлена вычисляеются одним методоном*/
    {
        complex <float> al = 0.5f * z2 - x1 - x, be = 4 * real(n) - z * real(m);
        complex <float> gm = sqrt(al + be), dl = sqrt(al - be);
        /*Вычисление самих корней многочлена*/
        complex<float> temp = (-0.5f * z + m) * 0.5f;//Временная переменная для менее затратного вычисления
        solution[0] = temp + gm * 0.5f;
        solution[2] = temp - gm * 0.5f;
        temp = temp - m;//изменение значения временной переменной
        solution[1] = temp + dl * 0.5f;
        solution[3] = temp - dl * 0.5f;
    }
    else/* Если m комплексная, то другим*/
    {
        //Если m комплесное число, то все корни являются комплесными
        m = complex <float>(imag(m), -real(m));
        n = complex <float>(imag(n), -real(n));
        complex <float> al = 0.5f * z2 - x1 - x, be = 4.0f * n - z * m;
        complex <float> p = sqrt(al * al + be * be), gm = sqrt(0.5f * (al + p));
        complex <float> dl;
        if (real(gm) == 0 && imag(gm) == 0)/* Если переменная gm равна 0, то переменная dl вычисляется по другому, а переменная al меняет свое значение */
        {
            al = (-real(al), -imag(al));
            dl = sqrt(al);
        }
        else
        {
            dl = 0.5f * be / gm;
        }
        complex <float> t = 0.5f * (m + dl);
        /*Вычисление самих корней многочлена*/
        float temp = -0.25f * z - imag(t);//Временная переменная для менее затратного вычисления
        solution[0] = complex <float>(temp + real(gm) * 0.5f, (real(t) + imag(gm) * 0.5f));
        solution[1] = conj(solution[0]);//Вычисление комплексно сопряженного корня
        t = 0.5f * (m - dl);
        solution[2] = complex <float>(temp - real(gm) * 0.5f, real(t) - imag(gm) * 0.5f);
        solution[3] = conj(solution[2]);//Вычисление комплексно сопряженного корня

    }
}