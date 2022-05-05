#include "polynomials.h"
#include <random>
#include <chrono>
#include <numeric>
#include "cardon.cpp"
#include "tiruneh.cpp"
#include "tomas_co.cpp"
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

    // estimating the error //
    vector<float> random_roots(3);
    vector<float> estimated_roots(3);
    vector<float> error_est_sum(100000); float sum_avg = 0;
    vector<float> error_est_max(100000); float max_avg = 0;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 3; ++j)
            random_roots[j] = udist(rng) / 1e8;
        third_degree_polynomial<float> P(random_roots);
        //P.info();
        estimated_roots = tomas_co(P);
        if (!isnan(estimated_roots.back())) //when we divide very small float it becomes nan and ruins the whole estimation process
        {
            sort(estimated_roots.begin(), estimated_roots.end());
            error_est_sum[i] = P.error_est_sum(estimated_roots);
            sum_avg += error_est_sum[i];
            error_est_max[i] = P.error_est_max(estimated_roots);
            max_avg += error_est_max[i];
        }
    }
    sort(error_est_sum.begin(), error_est_sum.end());
    cout << "Error estimation sum| x-x'|" << std::endl;
    cout << "Min : " << error_est_sum[0] << std::endl;
    cout << "Max : " << error_est_sum.back() << std::endl;
    sum_avg /= 10000; cout << "Average : " << sum_avg << std::endl;
    cout << "Median : " << error_est_sum[49999] << std::endl << std::endl;

    sort(error_est_max.begin(), error_est_max.end());
    cout << "Error estimation sum( | x-x'| / max(x,x') )" << std::endl;
    cout << "Min : " << error_est_max[0] << std::endl;
    cout << "Max : " << error_est_max.back() << std::endl;
    max_avg /= 10000; cout << "Average : " << max_avg << std::endl;
    cout << "Median : " << error_est_max[49999] << std::endl;
    
}