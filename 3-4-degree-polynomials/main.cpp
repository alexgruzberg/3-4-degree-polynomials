#include "finding_roots.h"
#include <random>
#include <chrono>

using namespace std;

//setting the range for the RNG//
uniform_int_distribution<mt19937::result_type> udist(0, (float)1e8);
mt19937 rng;

// Template for finding roots for various methods  //
template <typename pol, typename root_type>
inline float root_finder(pol& polynomial, vector<root_type>(*method)(pol), float& sum_avg, int& exceptions)
{
    vector<root_type> estimated_roots = method(polynomial);
    float error = polynomial.error_est_sum(estimated_roots);
    sum_avg += error;
    if (isnan(error) || isinf(error))
        throw division_by_zero();
    return error;
}

// Information about all error estimates  //
template <typename T>
inline void error_estimation_info(vector<float>& error_estimations, float& sum_avg, int& exceptions, int& tests, vector<T>& worst_case, float& max)
{
    sort(error_estimations.begin(), error_estimations.end());
    cout << "Error estimation sum| x-x'|" << endl;
    if (error_estimations.size() > 0) // There are some error estimations //
    {
        cout << "Min : " << error_estimations[0] << endl;
        cout << "Max : " << max << endl;
        sum_avg /= tests; cout << "Average : " << sum_avg << endl;
        cout << "Median : " << error_estimations[(tests - exceptions) / 2] << endl;
        cout << "Roots with the worst error estimation : ";
        for (auto v : worst_case)
            cout << v << " ";
    }
    cout << "Number of exceptions : " << exceptions << endl;
    cout << endl << endl << endl;

    error_estimations.clear(); sum_avg = 0; exceptions = 0; max = 0;
}

// Method for comparing the modulus of complex numbers  //
template <typename T>
bool complex_greater(complex<T> a, complex<T> b)
{
    return abs(a) > abs(b);
}

int main()
{
    // counting the current time in ms to use it as a seed for the rng   //
    auto time = chrono::system_clock::now();
    auto time_ms = chrono::time_point_cast<chrono::milliseconds>(time);
    auto value = time_ms.time_since_epoch();
    long dur = value.count();


    // giving the seed to our RNG //
    std::mt19937::result_type const seedval = dur;
    rng.seed(seedval);



    // Checking methods on different versions of polynomials //


    // estimating the error for polynomial with 3 different real roots //

    std::cout << "Cubic polynomials with 3 different real roots within the interval [0,x]" << std::endl;
    int tests; std::cout << "Enter the number of tests: "; std::cin >> tests;
    float range; std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;
    int exceptions = 0;
    vector<float> worst_case_cubic(3); float max = 0; float cur_error;

    vector<float> random_roots_cubic(3);
    vector<float> error_est_sum; float sum_avg = 0;

    vector<vector<float>(*)(third_degree_polynomial<float>)> estimating_functions_cubic(3);
    estimating_functions_cubic[0] = &cardon; estimating_functions_cubic[1] = &tiruneh; estimating_functions_cubic[2] = &tomas_co;

    vector<string> estimating_functions_cubic_names(estimating_functions_cubic.size());
    estimating_functions_cubic_names[0] = "Cardons method to solve a cubic equation";
    estimating_functions_cubic_names[1] = "A simplified expression for the solution of cubic polynomial equations using function evaluation-Tiruneh-2020";
    estimating_functions_cubic_names[2] = "Tomas Co Real roots for cubic equation";

    for (int k = 0; k < estimating_functions_cubic.size(); k++)
    {
        cout << "Method: " << estimating_functions_cubic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            for (int j = 0; j < 3; ++j)
                random_roots_cubic[j] = range * udist(rng) / (float)1e8;
            try
            {
                third_degree_polynomial<float> P(random_roots_cubic);
                cur_error = root_finder(P, estimating_functions_cubic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of cubics //
                {
                    max = cur_error;
                    worst_case_cubic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_cubic, max);
    }





    //Cubic polynomials with a complex conjugate pair //

    std::cout << "Cubic polynomials with real root within the interval [0,x] and a complex conjugate with Re z and |Im z| within the same interval." << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    vector<complex<float>> random_roots_cubic_c(3);
    vector<complex<float>> worst_case_cubic_c(3);

    vector<vector<complex<float>>(*)(third_degree_polynomial<complex<float>>)> estimating_functions_cubic_complex(3);
    estimating_functions_cubic_complex[0] = &cardon; estimating_functions_cubic_complex[1] = &tiruneh; estimating_functions_cubic_complex[2] = &tomas_co;
    vector<string> estimating_functions_cubic_complex_names(estimating_functions_cubic.size());
    estimating_functions_cubic_complex_names[0] = "Cardon’s method to solve a cubic equation";
    estimating_functions_cubic_complex_names[1] = "A simplified expression for the solution of cubic polynomial equations using function evaluation-Tiruneh-2020";
    estimating_functions_cubic_complex_names[2] = "modified for complex roots Tomas Co Real roots for cubic equation";

    for (int k = 0; k < estimating_functions_cubic_complex.size(); k++)
    {
        cout << "Method: " << estimating_functions_cubic_complex_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            float rand_real = range * udist(rng) / (float)1e8; float rand_img = range * udist(rng) / (float)1e8;
            random_roots_cubic_c[0] = complex<float>(rand_real, rand_img);
            random_roots_cubic_c[1] = complex<float>(rand_real, -rand_img);
            random_roots_cubic_c[2] = complex<float>(range * udist(rng) / (float)1e8, 0);
            try
            {
                third_degree_polynomial<complex<float>> P(random_roots_cubic_c);
                cur_error = root_finder(P, estimating_functions_cubic_complex[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of cubics //
                {
                    max = cur_error;
                    worst_case_cubic_c = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_cubic, max);
    }









    //---Cubic polynomials with 3 same real roots within the interval [0,x]---//

    std::cout << "Cubic polynomials with 3 same real roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions_cubic.size(); k++)
    {
        cout << "Method: " << estimating_functions_cubic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots_cubic[0] = range * udist(rng) / (float)1e8;
            random_roots_cubic[2] = random_roots_cubic[1] = random_roots_cubic[0];
            try
            {
                third_degree_polynomial<float> P(random_roots_cubic); cur_error = root_finder(P, estimating_functions_cubic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of cubics //
                {
                    max = cur_error;
                    worst_case_cubic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_cubic, max);
    }









    //---Cubic polynomials with 2 same real roots---//

    std::cout << "Cubic polynomials with 2 same real roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions_cubic.size(); ++k)
    {
        cout << "Method: " << estimating_functions_cubic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots_cubic[0] = range * udist(rng) / (float)1e8;
            random_roots_cubic[2] = random_roots_cubic[1] = range * udist(rng) / (float)1e8;
            try
            {
                third_degree_polynomial<float> P(random_roots_cubic);
                cur_error = root_finder(P, estimating_functions_cubic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of cubics //
                {
                    max = cur_error;
                    worst_case_cubic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_cubic, max);
    }









    //---Cubic polynomials with 3 almost similar real roots---//

    std::cout << "Cubic polynomials with 3 almost similar (<= 0.001 difference) real roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions_cubic.size(); ++k)
    {
        cout << "Method: " << estimating_functions_cubic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots_cubic[0] = range * udist(rng) / (float)1e8;
            random_roots_cubic[1] = random_roots_cubic[0] + range * udist(rng) / (float)(1000 * 1e8);
            random_roots_cubic[2] = random_roots_cubic[0] + range * udist(rng) / (float)(1000 * 1e8);
            try
            {
                third_degree_polynomial<float> P(random_roots_cubic);
                cur_error = root_finder(P, estimating_functions_cubic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of cubics //
                {
                    max = cur_error;
                    worst_case_cubic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_cubic, max);
    }









    //---Cubic polynomials with 2 almost similar real roots---//
    std::cout << "Cubic polynomials with 2 almost similar (<= 0.001 difference) real roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions_cubic.size(); ++k)
    {
        cout << "Method: " << estimating_functions_cubic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots_cubic[0] = range * udist(rng) / (float)1e8;
            random_roots_cubic[1] = random_roots_cubic[0] + range * udist(rng) / (float)(1000 * 1e8);
            random_roots_cubic[2] = range * udist(rng) / (float)1e8;
            try
            {
                third_degree_polynomial<float> P(random_roots_cubic);
                cur_error = root_finder(P, estimating_functions_cubic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of cubics //
                {
                    max = cur_error;
                    worst_case_cubic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_cubic, max);
    }









    //---Quartic polynomials with 4 different real roots---//
    std::cout << "Quartic polynomials with 4 different real roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;
    vector<float> worst_case_quartic(4);

    vector<float> random_roots_quartic(4);

    vector<vector<float>(*)(fourth_degree_polynomial<float>)> estimating_functions_quartic(1);
    estimating_functions_quartic[0] = &ferrari;

    vector<string> estimating_functions_quartic_names(estimating_functions_quartic.size());
    estimating_functions_quartic_names[0] = "Ferrari's solution for quartic equations";

    for (int k = 0; k < estimating_functions_quartic.size(); k++)
    {
        cout << "Method: " << estimating_functions_quartic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            for (int j = 0; j < 3; ++j)
                random_roots_quartic[j] = range * udist(rng) / (float)1e8;
            try
            {
                fourth_degree_polynomial<float> P(random_roots_quartic);
                cur_error = root_finder(P, estimating_functions_quartic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of quartics //
                {
                    max = cur_error;
                    worst_case_quartic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_quartic, max);
    }









    //---Quartic polynomials with 1 complex conjugate---//
    std::cout << "Quartic polynomials with one complex conjugate with Re z and |Im z| within the same interval. " << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    vector<complex<float>> random_roots_quartic_c(4);
    vector<complex<float>> worst_case_quartic_c(4);

    vector<vector<complex<float>>(*)(fourth_degree_polynomial<complex<float>>)> estimating_functions_quartic_complex(1);
    estimating_functions_quartic_complex[0] = &ferrari;

    vector<string> estimating_functions_quartic_complex_names(estimating_functions_quartic_complex.size());
    estimating_functions_quartic_complex_names[0] = "Ferrari's solution for quartic equations";

    for (int k = 0; k < estimating_functions_quartic_complex.size(); k++)
    {
        cout << "Method: " << estimating_functions_quartic_complex_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            float rand_real = range * udist(rng) / (float)1e8; float rand_img = range * udist(rng) / (float)1e8;
            random_roots_quartic_c[0] = complex<float>(rand_real, rand_img);
            random_roots_quartic_c[1] = complex<float>(rand_real, -rand_img);
            random_roots_quartic_c[2] = complex<float>(range * udist(rng) / (float)1e8, 0);
            random_roots_quartic_c[3] = complex<float>(range * udist(rng) / (float)1e8, 0);
            sort(random_roots_quartic_c.begin(), random_roots_quartic_c.end(), [](complex<float> a, complex<float> b) {return (a.real() > b.real()); });
            for (int i = 0; i < 3; ++i)
            {
                if (abs(random_roots_quartic_c[i].real() - random_roots_quartic_c[i + 1].real()) < 0.000001 && random_roots_quartic_c[i].imag() < -0.000001)
                {
                    swap(random_roots_quartic_c[i], random_roots_quartic_c[i + 1]);
                    ++i;
                }
            }
            try
            {
                fourth_degree_polynomial<complex<float>> P(random_roots_quartic_c);
                cur_error = root_finder(P, estimating_functions_quartic_complex[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of quartics //
                {
                    max = cur_error;
                    worst_case_quartic_c = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;  
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_quartic_c, max);
    }









    //---Quartic polynomials with 2 complex conjugates---//
    std::cout << "Quartic polynomials with one complex conjugate with Re z and |Im z| within the same interval. " << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions_quartic_complex.size(); k++)
    {
        cout << "Method: " << estimating_functions_quartic_complex_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            for (int i = 0; i < 2; ++i)
            {
                float rand_real = range * udist(rng) / (float)1e8; float rand_img = range * udist(rng) / (float)1e8;
                random_roots_quartic_c[2 * i] = complex<float>(rand_real, rand_img);
                random_roots_quartic_c[2 * i + 1] = complex<float>(rand_real, -rand_img);
            }
            sort(random_roots_quartic_c.begin(), random_roots_quartic_c.end(), [](complex<float> a, complex<float> b) {return (a.real() > b.real()); });
            for (int i = 0; i < 3; ++i)
            {
                if (abs(random_roots_quartic_c[i].real() - random_roots_quartic_c[i + 1].real()) < 0.000001 && random_roots_quartic_c[i].imag() < -0.000001)
                {
                    swap(random_roots_quartic_c[i], random_roots_quartic_c[i + 1]);
                    ++i;
                }
            }
            try
            {
                fourth_degree_polynomial<complex<float>> P(random_roots_quartic_c);
                cur_error = root_finder(P, estimating_functions_quartic_complex[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of quartics //
                {
                    max = cur_error;
                    worst_case_quartic_c = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl; 
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_quartic_c, max);
    }










    //---Quartic polynomials with 4 same roots---//
    std::cout << "Quartic polynomials with 4 same roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions_quartic.size(); k++)
    {
        cout << "Method: " << estimating_functions_quartic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots_quartic[0] = range * udist(rng) / (float)1e8;
            random_roots_quartic[3] = random_roots_quartic[2] = random_roots_quartic[1] = random_roots_quartic[0];
            try
            {
                fourth_degree_polynomial<float> P(random_roots_quartic);
                cur_error = root_finder(P, estimating_functions_quartic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of quartics //
                {
                    max = cur_error;
                    worst_case_quartic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_quartic, max);
    }









    //---Quartic polynomials with 3 same roots---//
    std::cout << "Quartic polynomials with 3 same roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions_quartic.size(); k++)
    {
        cout << "Method: " << estimating_functions_quartic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots_quartic[0] = range * udist(rng) / (float)1e8;
            random_roots_quartic[3] = random_roots_quartic[2] = random_roots_quartic[1] = range * udist(rng) / (float)1e8;
            try
            {
                fourth_degree_polynomial<float> P(random_roots_quartic);
                cur_error = root_finder(P, estimating_functions_quartic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of quartics //
                {
                    max = cur_error;
                    worst_case_quartic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_quartic, max);
    }









    //---Quartic polynomials with 2 same roots---//
    std::cout << "Quartic polynomials with 2 same roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions_quartic.size(); k++)
    {
        cout << "Method: " << estimating_functions_quartic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            for (int i = 0; i < 3; ++i)
                random_roots_quartic[i] = range * udist(rng) / (float)1e8;
            random_roots_quartic[3] = random_roots_quartic[2];
            try
            {
                fourth_degree_polynomial<float> P(random_roots_quartic);
                cur_error = root_finder(P, estimating_functions_quartic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of quartics //
                {
                    max = cur_error;
                    worst_case_quartic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_quartic, max);
    }









    //---Quartic polynomials with 4 almost same roots---//
    std::cout << "Quartic polynomials with 4 almost same roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions_quartic.size(); k++)
    {
        cout << "Method: " << estimating_functions_quartic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots_quartic[0] = range * udist(rng) / (float)1e8;
            random_roots_quartic[1] = random_roots_quartic[0] + range * udist(rng) / (float)(1000 * 1e8);
            random_roots_quartic[2] = random_roots_quartic[0] + range * udist(rng) / (float)(1000 * 1e8);
            random_roots_quartic[3] = random_roots_quartic[0] + range * udist(rng) / (float)(1000 * 1e8);
            try
            {
                fourth_degree_polynomial<float> P(random_roots_quartic);
                cur_error = root_finder(P, estimating_functions_quartic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of quartics //
                {
                    max = cur_error;
                    worst_case_quartic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_quartic, max);
    }









    //---Quartic polynomials with 3 almost similar roots---//
    std::cout << "Quartic polynomials with 3 almost same roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions_quartic.size(); k++)
    {
        cout << "Method: " << estimating_functions_quartic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots_quartic[0] = range * udist(rng) / (float)1e8;
            random_roots_quartic[1] = random_roots_quartic[0] + range * udist(rng) / (float)(1000 * 1e8);
            random_roots_quartic[2] = random_roots_quartic[0] + range * udist(rng) / (float)(1000 * 1e8);
            random_roots_quartic[3] = range * udist(rng) / (float)1e8;
            try
            {
                fourth_degree_polynomial<float> P(random_roots_quartic);
                cur_error = root_finder(P, estimating_functions_quartic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of quartics //
                {
                    max = cur_error;
                    worst_case_quartic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_quartic, max);
    }









    //---Quartic polynomials with 2 almost similar roots---//
    std::cout << "Quartic polynomials with 2 almost same roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions_quartic.size(); k++)
    {
        cout << "Method: " << estimating_functions_quartic_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            for (int i = 0; i < 3; ++i)
                random_roots_quartic[i] = range * udist(rng) / (float)1e8;
            random_roots_quartic[3] = random_roots_quartic[2] + range * udist(rng) / (float)(1000 * 1e8);
            try
            {
                fourth_degree_polynomial<float> P(random_roots_quartic);
                cur_error = root_finder(P, estimating_functions_quartic[k], sum_avg, exceptions);
                if (cur_error > max) // Trying to find the worst case of quartics //
                {
                    max = cur_error;
                    worst_case_quartic = P.get_roots();
                }
                error_est_sum.push_back(cur_error);
            }
            catch (division_by_zero& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (sqrt_of_negative_number& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (arccos_out_of_range& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
            catch (nan_value& e)
            {
                //std::cout << "exception caught" << std::endl;
                //std::cout << e.what() << std::endl;
                ++exceptions;
            }
        }
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case_quartic, max);
    }
}