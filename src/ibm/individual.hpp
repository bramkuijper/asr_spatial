#ifndef _INDIVIDUAL_CPP_
#define _INDIVIDUAL_CPP_

#include <random>
#include "parameters.hpp"

class Individual
{
    public:
        double T[2] = {0.0,0.0};
        int time_current_state = 0;

        Individual(double const Tf, double const Tm);

        Individual(Individual const &other);

        Individual(Individual const &mom
                ,Individual const &dad
                ,Parameters const &params
                ,std::mt19937 &rng_r);

};

#endif
