#include "patch.hpp"
#include "parameters.hpp"

Patch::Patch(
        Parameters const &parameters)
{
    envt2 = false;
    // make a vector of female breeders
    std::vector <Individual> breeders_f(parameters.nf_per_patch_init,
        Individual{parameters.init_Tf, parameters.init_Tm});

    // make a vector of male breeders
    std::vector <Individual> breeders_m(parameters.nf_per_patch_init,
        Individual{parameters.init_Tf, parameters.init_Tm});

    // make 2d vector containing both previous vecotrs
    adult_mate.push_back(breeders_f);
    adult_mate.push_back(breeders_m);
    adult_care.push_back(breeders_f);
    adult_care.push_back(breeders_m);
    juvenile.push_back(breeders_f);
    juvenile.push_back(breeders_m);
}

