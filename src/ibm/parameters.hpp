#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

#include <string>

enum Sex
{
    F = 0, 
    M = 1
};

class Parameters
{
    public:

        // initial allelic value of the amount of time spent 
        // caring for males and females
        double init_Tf = 20.0;
        double init_Tm = 20.0;

        // number of patches in the population
        int npatches = 1;

        // initial number of breeders in each patch
        int nf_per_patch_init = 2500;
        int nm_per_patch_init = 2500;

        // prefix of the simulation file
        std::string base_name = "sim_asr";

        long unsigned max_time_step = 500000;

        // mortalities
        double mu_mate[2] = {0.001,0.001}; // female and male mortalities while mating
        double mu_care[2] = {0.001,0.001}; // female and male mortalities while caring
        double mu_juv[2] = {0.001,0.001}; // female and male mortalities while juvenile

        // mutation rates
        double mutate_T[2] = {0.005,0.005}; // mutation rates
        double sdmu = 0.02; // mutational effect size distribution: standard deviation

        int max_time_juv[2] = {20,20};

        // synergy parameter
        double sigma = 0.0; 
        
        double D = 20.0;

        double gamma = 0.003;

        double prob_male = 0.5;

        double d[2] = {0.0,0.0};

        int output_interval = 1;
};

#endif
