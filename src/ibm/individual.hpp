#ifndef _INDIVIDUAL_CPP_
#define _INDIVIDUAL_CPP_

#include <random>
#include "parameters.hpp"

class Individual
{
    public:
        double T[2] = {0.0,0.0};
        double Tb[2] = {0.0,0.0};
        double phen = 0.0;
        int time_current_state = 0;
        
        double natal_sr = 0.0;
        bool natal_envt = false;

        Individual(double const Tf, double const Tm);

        Individual(Individual const &other);

        Individual(Individual const &mom
                ,Individual const &dad
                ,Parameters const &params
                ,double const natal_sr
                ,bool const natal_envt
                ,std::mt19937 &rng_r);

        void operator=(Individual const &other);
};

#endif
