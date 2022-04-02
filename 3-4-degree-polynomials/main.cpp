#include "polynomials.h"

int main()
{
    third_degree_polynomial P(1, 1, 1);
    P.info();
    third_degree_polynomial Q(5, 2, 3);
    Q.info();
    fourth_degree_polynomial T(1, 1, 1, 1);
    T.info();
    fourth_degree_polynomial S(4, 2, 5, 2);
    S.info();

}