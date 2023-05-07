#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include "parameters.hpp"
#include "patch.hpp"


class Simulation
{
    private:
        // parameter object
        Parameters params; 
        
        // a data file containing the results of the simulation
        std::ofstream data_file;

        // keep track of the time step of the simulation
        long unsigned time_step = 0;

        // random device which is used to generate
        // proper random seeds
        std::random_device rd;

        // store the random seed
        // we need to store this so that we can output the
        // random seed, so that we could 'replay' the exact
        // same sequence of random numbers for debugging purposes etc
        unsigned int seed;

        // random number generator
        std::mt19937 rng_r;

        // uniform distribution to compare against probabilities
        std::uniform_real_distribution<double> uniform;
        std::uniform_int_distribution<int> patch_sampler;

        // metapopulation consisting of patches
        std::vector<Patch> metapop;

    public:
        Simulation(Parameters const &params);

        void general_mortality();
        void mate_mortality(Patch &patch_i);
        void care_mortality(Patch &patch_i);
        void juvenile_mortality(Patch &patch_i);
        void maturation();
        void mating();
        double calculate_parental_care(Individual const &mom
                ,Individual const &dad);

        double care_survival_prob(double const Ttot, int const local_pop_size);
};

#endif
