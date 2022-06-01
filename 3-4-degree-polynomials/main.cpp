#include "finding_roots.h"
#include <random>
#include <chrono>

using namespace std;

//setting the range for the RNG//
uniform_int_distribution<mt19937::result_type> udist(0,(float)1e8);
mt19937 rng;

template <typename T>
inline float root_finder(third_degree_polynomial<T>& polynomial, vector<T>(*method)(third_degree_polynomial<T>), float& sum_avg, int& exceptions)
{
        vector<T> estimated_roots(3);
        estimated_roots = method(polynomial);
        float error = polynomial.error_est_sum(estimated_roots);
        sum_avg += error;
        if (isnan(error) || isinf(error))
            throw division_by_zero();
        return error;
}

template <typename T>
inline void error_estimation_info(vector<float>& error_estimations, float& sum_avg, int& exceptions, int& tests, vector<T>& worst_case, float& max)
{
    sort(error_estimations.begin(), error_estimations.end());
    cout << "Error estimation sum| x-x'|" << endl;
    cout << "Min : " << error_estimations[0] << endl;
    cout << "Max : " << max << endl;
    sum_avg /= tests; cout << "Average : " << sum_avg << endl;
    cout << "Median : " << error_estimations[(tests - exceptions) / 2] << endl;
    cout << "Number of exceptions : " << exceptions << endl;
    cout << "Roots with the worst error estimation : ";
    for (auto v : worst_case)
        cout << v << " ";
    cout << endl << endl << endl;

    error_estimations.clear(); sum_avg = 0; exceptions = 0; max = 0;
}

int main()
{
    
    // counting the current time in ms to use it as a seed for the rng   //
    auto time = chrono::system_clock::now();
    auto time_ms = chrono::time_point_cast<chrono::milliseconds>(time);
    auto value = time_ms.time_since_epoch();
    long dur = value.count();


    /*vector<complex<float>> crr(3);
    crr[0] = complex<float>(-0.62, 1.18);
    crr[1] = complex<float>(-0.62, -1.18);
    crr[2] = complex<float>(2.24, 0);
    third_degree_polynomial<complex<float>> pol(crr);
    vector<complex<float>> roooots = tomas_co(pol);
    for (int i = 0; i < 3; ++i)
        cout << crr[i] << " " << roooots[i] << endl;*/

    // giving the seed to our RNG //
    std::mt19937::result_type const seedval = dur;
    rng.seed(seedval);


    // estimating the error for polynomial with 3 different real roots //
    
    std::cout << "Cubic polynomials with 3 different real roots within the interval [0,x]" << std::endl;
    int tests; std::cout << "Enter the number of tests: "; std::cin >> tests;
    float range; std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;
    int exceptions = 0;
    vector<float> worst_case(3); float max = 0; float cur_error;

    vector<float> random_roots(3);
    vector<float> estimated_roots(3);
    vector<float> error_est_sum; float sum_avg = 0;

    vector<vector<float> (*)(third_degree_polynomial<float>)> estimating_functions(3);
    estimating_functions[0] = &cardon; estimating_functions[1] = &tiruneh; estimating_functions[2] = &tomas_co;

    vector<string> estimating_functions_names(estimating_functions.size());
    estimating_functions_names[0] = "Cardons method to solve a cubic equation";
    estimating_functions_names[1] = "A simplified expression for the solution of cubic polynomial equations using function evaluation-Tiruneh-2020";
    estimating_functions_names[2] = "Tomas Co Real roots for cubic equation";

    for (int k = 0; k < estimating_functions.size(); k++)
    {
        cout << "Method: " << estimating_functions_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            for (int j = 0; j < 3; ++j)
                random_roots[j] = range * udist(rng) / (float)1e8;
            try
            {
                third_degree_polynomial<float> P(random_roots);
                cur_error = root_finder(P, estimating_functions[k], sum_avg, exceptions);
                if (cur_error > max)
                {
                    max = cur_error;
                    worst_case = P.get_roots();
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
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case, max);
    }
    




    


    //Cubic polynomials with a complex conjugate pair //

    std::cout << "Cubic polynomials with real root within the interval [0,x] and a complex conjugate with Re z and |Im z| within the same interval." << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;
    exceptions = 0;

    vector<complex<float>> random_roots_c(3);
    vector<complex<float>> estimated_roots_c(3);
    vector<complex<float>> worst_case_c(3);

    vector<vector<complex<float>>(*)(third_degree_polynomial<complex<float>>)> estimating_functions_complex(3);
    estimating_functions_complex[0] = &cardon; estimating_functions_complex[1] = &tiruneh; estimating_functions_complex[2] = &tomas_co;
    vector<string> estimating_functions_complex_names(estimating_functions.size());
    estimating_functions_complex_names[0] = "Cardon’s method to solve a cubic equation";
    estimating_functions_complex_names[1] = "A simplified expression for the solution of cubic polynomial equations using function evaluation-Tiruneh-2020";
    estimating_functions_complex_names[2] = "modified for complex roots Tomas Co Real roots for cubic equation";

    for (int k = 0; k < estimating_functions_complex.size(); k++)
    {
        cout << "Method: " << estimating_functions_complex_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            float rand_real = range * udist(rng) / (float)1e8; float rand_img = range * udist(rng) / (float)1e8;
            random_roots_c[0] = complex<float>(rand_real, rand_img);
            random_roots_c[1] = complex<float>(rand_real, -rand_img);
            random_roots_c[2] = complex<float>(range * udist(rng) / (float)1e8, 0);
            try
            {
                third_degree_polynomial<complex<float>> P(random_roots_c);
                cur_error = root_finder(P, estimating_functions_complex[k], sum_avg, exceptions);
                if (cur_error > max)
                {
                    max = cur_error;
                    worst_case_c = P.get_roots();
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
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case, max);
    }









    //---Cubic polynomials with 3 same real roots within the interval [0,x]---//

    std::cout << "Cubic polynomials with 3 same real roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions.size(); k++)
    {
        cout << "Method: " << estimating_functions_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots[0] = range * udist(rng) / (float)1e8;
            random_roots[2] = random_roots[1] = random_roots[0];
            try
            {
                third_degree_polynomial<float> P(random_roots); cur_error = root_finder(P, estimating_functions[k], sum_avg, exceptions);
                if (cur_error > max)
                {
                    max = cur_error;
                    worst_case = P.get_roots();
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
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case, max);
    }








    
    //---Cubic polynomials with 2 same real roots---//
    
    std::cout << "Cubic polynomials with 2 same real roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions.size(); ++k)
    {
        cout << "Method: " << estimating_functions_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots[0] = range * udist(rng) / (float)1e8;
            random_roots[2] = random_roots[1] = range * udist(rng) / (float)1e8;
            try
            {
                third_degree_polynomial<float> P(random_roots); 
                cur_error = root_finder(P, estimating_functions[k], sum_avg, exceptions);
                if (cur_error > max)
                {
                    max = cur_error;
                    worst_case = P.get_roots();
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
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case, max);
    }
    








    //---Cubic polynomials with 3 almost similar real roots---//

    std::cout << "Cubic polynomials with 3 almost similar (<= 0.001 difference) real roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions.size(); ++k)
    {
        cout << "Method: " << estimating_functions_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots[0] = range * udist(rng) / (float)1e8;
            random_roots[1] = random_roots[0] + range * udist(rng) / (float)(1000 * 1e8);
            random_roots[2] = random_roots[0] + range * udist(rng) / (float)(1000 * 1e8);
            try
            {
                third_degree_polynomial<float> P(random_roots); 
                cur_error = root_finder(P, estimating_functions[k], sum_avg, exceptions);
                if (cur_error > max)
                {
                    max = cur_error;
                    worst_case = P.get_roots();
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
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case, max);
    }
    


    





    //---Cubic polynomials with 2 almost similar real roots---//
    std::cout << "Cubic polynomials with 2 almost similar (<= 0.001 difference) real roots within the interval [0,x]" << std::endl;
    std::cout << "Enter the number of tests: "; std::cin >> tests;
    std::cout << "Enter the maximum value of a root: "; std::cin >> range; cout << endl << endl;

    for (int k = 0; k < estimating_functions.size(); ++k)
    {
        cout << "Method: " << estimating_functions_names[k] << std::endl;
        for (int i = 0; i < tests; ++i)
        {
            random_roots[0] = range * udist(rng) / (float)1e8;
            random_roots[1] = random_roots[0] + range * udist(rng) / (float)(1000 * 1e8);
            random_roots[2] = range * udist(rng) / (float)1e8;
            try
            {
                third_degree_polynomial<float> P(random_roots);
                cur_error = root_finder(P, estimating_functions[k], sum_avg, exceptions);
                if (cur_error > max)
                {
                    max = cur_error;
                    worst_case = P.get_roots();
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
        error_estimation_info(error_est_sum, sum_avg, exceptions, tests, worst_case, max);
    }
}