#include "polynomials.h"
#include <random>
#include <chrono>
using namespace std;

//setting the range for the RNG//
uniform_int_distribution<mt19937::result_type> udist(0,1e8);
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

    // do something //
    vector<float> random_roots(3);
    vector<float> estimated_roots(3);
    float error = 0;
    for (int i = 0; i < 10; ++i)
    {
        for (int j = 0; j < 3; ++j)
            random_roots[j] = udist(rng) / 1e8;
        third_degree_polynomial P(random_roots);
        P.info();
        estimated_roots = P.tiruneh();
        for (auto v : estimated_roots) std::cout << v << " "; std::cout << std::endl;
        // need to measure the error
    }

}