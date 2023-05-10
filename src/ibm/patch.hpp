#ifndef _PATCH_HPP_
#define _PATCH_HPP_

#include <vector>
#include "individual.hpp"
#include "parameters.hpp"

class Patch
{
    public:
        Patch(Parameters const &params);
        std::vector < std::vector<Individual> > adult_care;
        std::vector < std::vector<Individual> > adult_mate;
        std::vector < std::vector<Individual> > juvenile;

        bool envt2 = false;

};

#endif
