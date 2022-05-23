#include "finding_roots.h"
#include <random>
#include <chrono>

using namespace std;

//setting the range for the RNG//
uniform_int_distribution<mt19937::result_type> udist(0,(float)1e8);
mt19937 rng;


int main()
{
    // counting the current time in ms to use it as a seed for the rng //
    auto time = chrono::system_clock::now();
    auto time_ms = chrono::time_point_cast<chrono::milliseconds>(time);
    auto value = time_ms.time_since_epoch();
    long dur = value.count();

    // giving the seed to our RNG //
    std::mt19937::result_type const seedval = dur;
    rng.seed(seedval);


    // estimating the error for polynomial with 3 different real roots //
    
    std::cout << "Cubic polynomials with 3 different real roots within the interval [0,x]" << std::endl;
    int tests; std::cout << "Enter the number of tests: "; std::cin >> tests;
    float range; std::cout << "Enter the maximum value of a root: "; std::cin >> range;
    int exceptions = 0;

    vector<float> random_roots(3);
    vector<float> estimated_roots(3);
    vector<float> error_est_sum(tests); float sum_avg = 0;
    vector<float> error_est_max(tests); float max_avg = 0;
    for (int i = 0; i < tests; ++i)
    {
        for (int j = 0; j < 3; ++j)
            random_roots[j] = range * udist(rng) / (float)1e8;
        third_degree_polynomial<float> P(random_roots);
        //P.info();
        try
        {
        estimated_roots = cardon(P);
        sort(estimated_roots.begin(), estimated_roots.end());
        error_est_sum[i] = P.error_est_sum(estimated_roots);
        sum_avg += error_est_sum[i];
        error_est_max[i] = P.error_est_max(estimated_roots);
        max_avg += error_est_max[i];
        }
        catch (division_by_zero& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (sqrt_of_negative_number& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (arccos_out_of_range& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
    }
    sort(error_est_sum.begin(), error_est_sum.end());
    cout << "Error estimation sum| x-x'|" << std::endl;
    cout << "Min : " << error_est_sum[0] << std::endl;
    cout << "Max : " << error_est_sum.back() << std::endl;
    sum_avg /= tests; cout << "Average : " << sum_avg << std::endl;
    cout << "Median : " << error_est_sum[tests / 2] << std::endl << std::endl;

    sort(error_est_max.begin(), error_est_max.end());
    cout << "Error estimation sum( | x-x'| / max(x,x') )" << std::endl;
    cout << "Min : " << error_est_max[0] << std::endl;
    cout << "Max : " << error_est_max.back() << std::endl;
    max_avg /= tests; cout << "Average : " << max_avg << std::endl;
    cout << "Median : " << error_est_max[tests / 2] << std::endl;

    cout << "Number of exceptions : " << exceptions << std::endl;

    







    //Cubic polynomials with a complex conjugate pair //

    /*std::cout << "Cubic polynomials with real root within the interval [0,x] and a complex conjugate with Re z and |Im z| within the same interval." << std::endl;
    int tests; std::cout << "Enter the number of tests: "; std::cin >> tests;
    float range; std::cout << "Enter the maximum value of a root: "; std::cin >> range;
    int exceptions = 0;

    vector<complex<float>> random_roots(3);
    vector<complex<float>> estimated_roots(3);
    vector<float> error_est_sum(tests); float sum_avg = 0;
    vector<float> error_est_max(tests); float max_avg = 0;
    for (int i = 0; i < tests; ++i)
    {
        float rand_real =  range * udist(rng) / (float)1e8; float rand_img = range * udist(rng) / (float)1e8;
        random_roots[0] = complex<float>(rand_real, rand_img);
        random_roots[1] = complex<float>(rand_real, -rand_img);
        random_roots[2] = complex<float>(range * udist(rng) / (float)1e8,0);
        third_degree_polynomial<complex<float>> P(random_roots);
        try 
        {
            estimated_roots = tiruneh(P);
            error_est_sum[i] = P.error_est_sum(estimated_roots);
            sum_avg += error_est_sum[i];
        }
        catch (division_by_zero& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (sqrt_of_negative_number& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (arccos_out_of_range& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
    }
    sort(error_est_sum.begin(), error_est_sum.end());
    cout << "Error estimation sum| x-x'|" << std::endl;
    cout << "Min : " << error_est_sum[0] << std::endl;
    cout << "Max : " << error_est_sum.back() << std::endl;
    sum_avg /= tests; cout << "Average : " << sum_avg << std::endl;
    cout << "Median : " << error_est_sum[tests / 2] << std::endl << std::endl;

    cout << "Number of exceptions : " << exceptions << std::endl;*/









    //---Cubic polynomials with 3 same real roots within the interval [0,x]---//

    /*std::cout << "Cubic polynomials with 3 same real roots within the interval [0,x]" << std::endl;
    int tests; std::cout << "Enter the number of tests: "; std::cin >> tests;
    float range; std::cout << "Enter the maximum value of a root: "; std::cin >> range;
    int exceptions = 0; 

    vector<float> random_roots(3);
    vector<float> estimated_roots(3);
    vector<float> error_est_sum(tests); float sum_avg = 0;
    vector<float> error_est_max(tests); float max_avg = 0;
    for (int i = 0; i < tests; ++i)
    {
        random_roots[0] = range * udist(rng) / (float)1e8;
        random_roots[2] = random_roots[1] = random_roots[0];
        third_degree_polynomial<float> P(random_roots);
        //P.info();
        try
        {
        estimated_roots = tiruneh(P);
        sort(estimated_roots.begin(), estimated_roots.end());
        error_est_sum[i] = P.error_est_sum(estimated_roots);
        sum_avg += error_est_sum[i];
        error_est_max[i] = P.error_est_max(estimated_roots);
        max_avg += error_est_max[i];
        }
        catch (division_by_zero& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (sqrt_of_negative_number& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (arccos_out_of_range& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
    }
    sort(error_est_sum.begin(), error_est_sum.end());
    cout << "Error estimation sum| x-x'|" << std::endl;
    cout << "Min : " << error_est_sum[0] << std::endl;
    cout << "Max : " << error_est_sum.back() << std::endl;
    sum_avg /= tests; cout << "Average : " << sum_avg << std::endl;
    cout << "Median : " << error_est_sum[tests / 2] << std::endl << std::endl;

    sort(error_est_max.begin(), error_est_max.end());
    cout << "Error estimation sum( | x-x'| / max(x,x') )" << std::endl;
    cout << "Min : " << error_est_max[0] << std::endl;
    cout << "Max : " << error_est_max.back() << std::endl;
    max_avg /= tests; cout << "Average : " << max_avg << std::endl;
    cout << "Median : " << error_est_max[tests / 2] << std::endl;

    cout << "Number of exceptions : " << exceptions << std::endl;*/








    
    //---Cubic polynomials with 2 same real roots---//
    
    /*std::cout << "Cubic polynomials with 2 same real roots within the interval [0,x]" << std::endl;
    int tests; std::cout << "Enter the number of tests: "; std::cin >> tests;
    float range; std::cout << "Enter the maximum value of a root: "; std::cin >> range;
    int exceptions = 0; 

    vector<float> random_roots(3);
    vector<float> estimated_roots(3);
    vector<float> error_est_sum(tests); float sum_avg = 0;
    vector<float> error_est_max(tests); float max_avg = 0;
    for (int i = 0; i < tests; ++i)
    {
        random_roots[0] = range * udist(rng) / (float)1e8;
        random_roots[2] = random_roots[1] = range * udist(rng) / (float)1e8;
        third_degree_polynomial<float> P(random_roots);
        //P.info();
        try
        {
        estimated_roots = tomas_co(P);
        sort(estimated_roots.begin(), estimated_roots.end());
        error_est_sum[i] = P.error_est_sum(estimated_roots);
        sum_avg += error_est_sum[i];
        error_est_max[i] = P.error_est_max(estimated_roots);
        max_avg += error_est_max[i];
        }
        catch (division_by_zero& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (sqrt_of_negative_number& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (arccos_out_of_range& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
    }
    sort(error_est_sum.begin(), error_est_sum.end());
    cout << "Error estimation sum| x-x'|" << std::endl;
    cout << "Min : " << error_est_sum[0] << std::endl;
    cout << "Max : " << error_est_sum.back() << std::endl;
    sum_avg /= tests; cout << "Average : " << sum_avg << std::endl;
    cout << "Median : " << error_est_sum[tests / 2] << std::endl << std::endl;

    sort(error_est_max.begin(), error_est_max.end());
    cout << "Error estimation sum( | x-x'| / max(x,x') )" << std::endl;
    cout << "Min : " << error_est_max[0] << std::endl;
    cout << "Max : " << error_est_max.back() << std::endl;
    max_avg /= tests; cout << "Average : " << max_avg << std::endl;
    cout << "Median : " << error_est_max[tests / 2] << std::endl;

    cout << "Number of exceptions : " << exceptions << std::endl;*/
    








    //---Cubic polynomials with 3 almost similar real roots---//

    /*std::cout << "Cubic polynomials with 3 almost similar (<= 0.001 difference) real roots within the interval [0,x]" << std::endl;
    int tests; std::cout << "Enter the number of tests: "; std::cin >> tests;
    float range; std::cout << "Enter the maximum value of a root: "; std::cin >> range;
    int exceptions = 0; 

    vector<float> random_roots(3);
    vector<float> estimated_roots(3);
    vector<float> error_est_sum(tests); float sum_avg = 0;
    vector<float> error_est_max(tests); float max_avg = 0;
    for (int i = 0; i < tests; ++i)
    {
        random_roots[0] = range * udist(rng) / (float)1e8;
        random_roots[1] = random_roots[0] + range * udist(rng) / (float)(1000 * 1e8);
        random_roots[2] = random_roots[0] + range * udist(rng) / (float)(1000 * 1e8);
        third_degree_polynomial<float> P(random_roots);
        //P.info();
        try
        {
        estimated_roots = tiruneh(P);
        sort(estimated_roots.begin(), estimated_roots.end());
        error_est_sum[i] = P.error_est_sum(estimated_roots);
        sum_avg += error_est_sum[i];
        error_est_max[i] = P.error_est_max(estimated_roots);
        max_avg += error_est_max[i];
        }
        catch (division_by_zero& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (sqrt_of_negative_number& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (arccos_out_of_range& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
    }
    sort(error_est_sum.begin(), error_est_sum.end());
    cout << "Error estimation sum| x-x'|" << std::endl;
    cout << "Min : " << error_est_sum[0] << std::endl;
    cout << "Max : " << error_est_sum.back() << std::endl;
    sum_avg /= tests; cout << "Average : " << sum_avg << std::endl;
    cout << "Median : " << error_est_sum[tests / 2] << std::endl << std::endl;

    sort(error_est_max.begin(), error_est_max.end());
    cout << "Error estimation sum( | x-x'| / max(x,x') )" << std::endl;
    cout << "Min : " << error_est_max[0] << std::endl;
    cout << "Max : " << error_est_max.back() << std::endl;
    max_avg /= tests; cout << "Average : " << max_avg << std::endl;
    cout << "Median : " << error_est_max[tests / 2] << std::endl;

    cout << "Number of exceptions : " << exceptions << std::endl;*/
    


    





    //---Cubic polynomials with 2 almost similar real roots---//
    /*std::cout << "Cubic polynomials with 2 almost similar (<= 0.001 difference) real roots within the interval [0,x]" << std::endl;
    int tests; std::cout << "Enter the number of tests: "; std::cin >> tests;
    float range; std::cout << "Enter the maximum value of a root: "; std::cin >> range;
    int exceptions = 0; 

    vector<float> random_roots(3);
    vector<float> estimated_roots(3);
    vector<float> error_est_sum(tests); float sum_avg = 0;
    vector<float> error_est_max(tests); float max_avg = 0;
    for (int i = 0; i < tests; ++i)
    {
        random_roots[0] = range * udist(rng) / (float)1e8;
        random_roots[1] = random_roots[0] + range * udist(rng) / (float)(1000 * 1e8);
        random_roots[2] = range * udist(rng) / (float)1e8;
        third_degree_polynomial<float> P(random_roots);
        //P.info();
        try
        {
        estimated_roots = tomas_co(P);
        sort(estimated_roots.begin(), estimated_roots.end());
        error_est_sum[i] = P.error_est_sum(estimated_roots);
        sum_avg += error_est_sum[i];
        error_est_max[i] = P.error_est_max(estimated_roots);
        max_avg += error_est_max[i];
        }
        catch (division_by_zero& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (sqrt_of_negative_number& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
        catch (arccos_out_of_range& e)
        {
            std::cout << "exception caught" << std::endl;
            std::cout << e.what() << std::endl;
            ++exceptions;
        }
    }
    sort(error_est_sum.begin(), error_est_sum.end());
    cout << "Error estimation sum| x-x'|" << std::endl;
    cout << "Min : " << error_est_sum[0] << std::endl;
    cout << "Max : " << error_est_sum.back() << std::endl;
    sum_avg /= tests; cout << "Average : " << sum_avg << std::endl;
    cout << "Median : " << error_est_sum[tests / 2] << std::endl << std::endl;

    sort(error_est_max.begin(), error_est_max.end());
    cout << "Error estimation sum( | x-x'| / max(x,x') )" << std::endl;
    cout << "Min : " << error_est_max[0] << std::endl;
    cout << "Max : " << error_est_max.back() << std::endl;
    max_avg /= tests; cout << "Average : " << max_avg << std::endl;
    cout << "Median : " << error_est_max[tests / 2] << std::endl;

    cout << "Number of exceptions : " << exceptions << std::endl;*/

    /*experiment
    float x = 1e8 * 1e8, y = 1 / (float)1e8; x = x * x; y = y * y;
    if (isnan(x / y))
        std::cout << "NaN!";
    else if (isinf(x / y))
        std::cout << "inf!";
    else
        std::cout << "What an unpleasant surpise! x/y = " << (x/y);

    //received inf - you were right

    */
}