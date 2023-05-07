#include <stdexcept>
#include <random>
#include <vector>
#include <cassert>
#include <algorithm>
#include "simulation.hpp"
#include "patch.hpp"


// constructor: starts a simulation
// Sets up a simulation object that runs the simulation
//
// @param params_arg A Parameters struct containing all the parameters
//      of the current simulation
Simulation::Simulation(Parameters const &params_arg) :
    params{params_arg} // store parameters as a data member
    ,rd{} // initialize random device, see *.h file
    ,seed{rd()} // initialize seed
    ,rng_r{seed} // initialize the random number generator
    ,uniform{0.0,1.0} // initialize the uniform distribution
    ,patch_sampler{params.npatches - 1} // uniform distribution to sample patches
    ,metapop(params.npatches,Patch(params)) // set up all the patches of the metapop
    ,data_file{params.base_name.c_str()} // start the file to write data to
{
    // go over all time steps
    for (time_step = 0; 
            time_step < params.max_time_step; ++time_step)
    {
        general_mortality();

        maturation();

        mating();
    }
} // end Simulation::Simulation()

// calculates the amount of parental care (in terms of time steps)
// produced by a pair of Individuals
//
// @param mom An Individual object that is the mother
// @param dad An Individual object that is the father

double Simulation::calculate_parental_care(Individual const &mom
                ,Individual const &dad)
{
    // get Tm and Tf values from dad and mom respectively
    double Tm = dad.T[M];
    double Tf = mom.T[F];

    assert(Tm >= 0.0);
    assert(Tf >= 0.0);

    // return the amount of time spent caring
    return(Tm + Tf + params.sigma * Tf * Tm);
} // end calculate_parental_care


// calculates the survival probability of an individual that has received
// Ttot care in a population with population size local_pop_size
//
// @param Ttot total care received
// @param local_pop_size local population size of the patch (for density dependence)
double Simulation::care_survival_prob(double const Ttot, int const local_pop_size)
{
    // S0 * S(Ttot, N) - see L206, 208 in Long et al 2022 bioRxiv
    return(Ttot * Ttot / (
                (Ttot * Ttot + params.D * params.D) * 
                (1.0 + params.gamma * local_pop_size)
                )
            );
} // end care_survival_prob

// let individuals randomly choose a mate: each time step each individual in the mating
// state has one mating opportunity. Pairs are formed at random until only 
// one sex is left in the mating state
//
void Simulation::mating()
{
    // make a vector with integers then shuffle
    int mother_idx, father_idx;

    // aux variable to calculate amount of parental care
    double parental_care;

    // aux variable to store patch id of destination patch
    int destination_patch;

    for (int patch_idx = 0; patch_idx < params.npatches; ++patch_idx)
    {
        // see which one is the limiting sex
        int nmales = metapop[patch_idx].adult_mate[M].size();
        int nfemales = metapop[patch_idx].adult_mate[F].size();

        // get numbers & identifiers for common and rare sexes
        bool common_sex_is_male = nmales >= nfemales;

        // identifiers for the rare and common sex in this patch
        Sex rare_sex = common_sex_is_male ? F : M;
        Sex common_sex = common_sex_is_male ? M : F;

        int ncommon_sex = common_sex_is_male ? nmales : nfemales;
        int nrare_sex = common_sex_is_male ? nfemales : nmales;

        // total number of individuals in the local patch
        int local_density = nmales + nfemales + 
            metapop[patch_idx].adult_care[M].size() + 
            metapop[patch_idx].adult_care[F].size() + 
            metapop[patch_idx].juvenile[M].size() + 
            metapop[patch_idx].juvenile[F].size();

        //// make a vector with unique indices of the common sex
        //std::vector <int> random_partner(0, ncommon_sex);

        //// fill vector
        //for (int common_sex_idx = 0; 
        //        common_sex_idx < ncommon_sex; ++common_sex_idx)
        //{
        //    random_partner.push_back(common_sex_idx);
        //}

        // shuffle list elements containing members of the common sex
        std::shuffle(
                metapop[patch_idx].adult_mate[common_sex].begin(), 
                metapop[patch_idx].adult_mate[common_sex].end(), 
                rng_r);

        // we don't need to shuffle the rare sex (as the common sex is already
        // shuffled, guaranteeing random mating)
        for (int rare_sex_idx = 0;
                rare_sex_idx < nrare_sex
                ++rare_sex_idx
                )
        {
            // one of these arrays (the common sex is randomly shuffled
            // the other is simply ordered). This should be fine.
            parental_care = calculate_parental_care(
                    metapop[patch_idx].adult_mate[F][rare_sex_idx],
                    metapop[patch_idx].adult_mate[M][rare_sex_idx]);

            // check whether individual survives
            if (uniform(rng_r) < care_survival_prob(parental_care, local_density))
            {
                // make a kid via the birth constructor
                Individual Kid(
                    metapop[patch_idx].adult_mate[F][rare_sex_idx]
                    ,metapop[patch_idx].adult_mate[M][rare_sex_idx]
                    ,params
                    ,rng_r
                    );

                // determine sex
                Sex offspring_sex = uniform(rng_r) < params.prob_male ? M : F;

                // by default local patch is the destination
                destination_patch = patch_idx;

                // but offspring might disperse to a random patch
                if (uniform(rng_r) < params.d[offspring_sex] && params.npatches > 1)
                {
                    destination_patch = patch_sampler(rng_r);
                }

                // add offspring to stack
                metapop[destination_patch].adult_mate[offspring_sex].push_back(Kid);
            }
        } //  end for int rare sex

        // all done - now we can move individuals to caring state

        // first members of rare sex - all are caring now
        for (
                std::vector<Individual>::iterator ind_it = 
                metapop[patch_idx].adult_mate[rare_sex].begin();
                ind_it != metapop[patch_idx].adult_mate[rare_sex].end();
                // no increment as we erase
                )
        {
            // set their state to 0
            ind_it->time_current_state = 0;

            // add to caring pool, we use * to get at the actual
            // individual pointed to by the iterator
            metapop[patch_idx].adult_care[rare_sex].push_back(*ind_it);

            // remove from mating pool
            ind_it = metapop[patch_idx].adult_mate[rare_sex].erase(ind_it);
        }

        // remove the first members nrare_sex individuals 
        // of the common sex that have now mated
        // these are the members [0, nrare_sex) 
        // of the random_partner[] vector
        for (
            std::vector<Individual>::iterator ind_it = 
                    metapop[patch_idx].adult_mate[common_sex].begin();
                    ind_it != metapop[patch_idx].adult_mate[common_sex].begin() + nrare_sex;
        )
        {
            // add to care pop
            ind_it->time_current_state = 0;
            metapop[patch_idx].adult_care[common_sex].push_back(*ind_it);

            ind_it = metapop[patch_idx].adult_mate[common_sex].erase(ind_it);
        } // common sex individuals removed
    } // end for patch_idx
} // mating()

// check which juvs are ready to mature
void Simulation::maturation()
{
    for (int patch_idx = 0; patch_idx < params.npatches; ++patch_idx)
    {
        for (int sex_idx  = 0; sex_idx < 2; ++sex_idx)
        {
            for (std::vector<Individual>::iterator ind_it = 
                    metapop[patch_idx].juvenile[sex_idx].begin();
                    ind_it != metapop[patch_idx].juvenile[sex_idx].end();
                    )
            {
                // individual has been juvenile for longer than max time steps
                // put into adult_mate and erase from the juvenile stack
                if (ind_it->time_current_state >= params.max_time_juv[sex_idx])
                {
                    // set time in current state already to 0 as it is going to 
                    // be an adult for 0 days
                    ind_it->time_current_state = 0;

                    // copy over to adult mate stack
                    metapop[patch_idx].adult_mate[sex_idx].push_back(*ind_it);

                    // remove from juv stack
                    metapop[patch_idx].juvenile[sex_idx].erase(ind_it);

                    // check whether this corresponds to the last element of the stack
                    assert(metapop[patch_idx].adult_mate[
                            sex_idx][
                                    metapop[patch_idx].adult_mate[sex_idx].size() - 1
                                        ].time_current_state == 0);
                }
                else
                {
                    ++ind_it;
                }
            }
        } // end for sex idx
    } // end for patch_idx
} // end Simulation::maturation()

void Simulation::juvenile_mortality(
        Patch &patch_i)
{
    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        // juvenile mortality
        for (std::vector<Individual>::iterator ind_it = 
                patch_i.juvenile[sex_idx].begin();
                ind_it != patch_i.juvenile[sex_idx].end();
                )
        {
            if (uniform(rng_r) < params.mu_juv[sex_idx])
            {
                // remove individual from list
                patch_i.juvenile[sex_idx].erase(ind_it);
            }
            else
            {
                ind_it->time_current_state++;
                // no removal? Increase element
                ++ind_it;
            }
            
        } // for individual
    } // for sex idx
} // end juv_mortality()

void Simulation::mate_mortality(
        Patch &patch_i)
{
    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        // mate mortality
        for (std::vector<Individual>::iterator ind_it = 
                patch_i.adult_mate[sex_idx].begin();
                ind_it != patch_i.adult_mate[sex_idx].end();
                )
        {
            if (uniform(rng_r) < params.mu_mate[sex_idx])
            {
                // remove individual from list
                patch_i.adult_mate[sex_idx].erase(ind_it);
            }
            else
            {
                ind_it->time_current_state++;
                // no removal? Increase element
                ++ind_it;
            }
            
        } // for individual
    } // for sex idx
} // end mate_mortality()


void Simulation::general_mortality()
{
    for (int patch_idx = 0; 
            patch_idx < params.npatches; ++patch_idx)
    {
        mate_mortality(metapop[patch_idx]);
        care_mortality(metapop[patch_idx]);
    } // for patches
} // end general_mortality()

void Simulation::care_mortality(
        Patch &patch_i)
{
    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        // mate mortality
        for (std::vector<Individual>::iterator ind_it = 
                patch_i.adult_care[sex_idx].begin();
                ind_it != patch_i.adult_care[sex_idx].end();
                )
        {
            if (uniform(rng_r) < params.mu_care[sex_idx])
            {
                // remove individual from list
                patch_i.adult_care[sex_idx].erase(ind_it);
            }
            else
            {
                // no removal? Increase element
                ind_it->time_current_state++;
                ++ind_it;
            }
            
        } // for individual
    } // for sex idx
} // end care_mortality()
