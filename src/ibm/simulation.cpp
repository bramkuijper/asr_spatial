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
    write_data_headers();

    initialize_environmental_variation();

    // go over all time steps
    for (time_step = 0; 
            time_step < params.max_time_step; ++time_step)
    {
        // remove all individuals whose care has ended
        // and put them again in mating pool
        care_ends();

        general_mortality();

        // maturation of juveniles
        maturation();

        mating();

        environmental_change();

        if (time_step % params.output_interval == 0)
        {
            write_data();
        }
    }

    write_parameters();
} // end Simulation::Simulation()


void Simulation::environmental_change()
{
    double prob_envt2 = params.environmental_change[0] / 
        (params.environmental_change[0] + params.environmental_change[1]);

    // get current envt from first patch
    bool current_envt = metapop[0].envt2;

    bool envt2_if_coarse = 
        uniform(rng_r) < params.environmental_change[current_envt] ? !current_envt : current_envt;

    for (std::vector<Patch>::iterator patch_iter = metapop.begin();
            patch_iter != metapop.end();
            ++patch_iter)
    {
        current_envt = patch_iter->envt2;

        // if envtal grain is coarse, give it general envtal state
        if (params.environment_coarse)
        {
            patch_iter->envt2 = envt2_if_coarse;
        }
        else if (uniform(rng_r) < params.environmental_change[current_envt])
        {
            patch_iter->envt2 = !patch_iter->envt2;
        }
    }
} // end Simulation::environmental_change

// initialize the environment of each patch
void Simulation::initialize_environmental_variation()
{
    double prob_envt2 = params.init_prob_envt2;

    bool envt2_if_coarse = uniform(rng_r) < prob_envt2;

    for (std::vector<Patch>::iterator patch_iter = metapop.begin();
            patch_iter != metapop.end();
            ++patch_iter)
    {
        // if envtal grain is coarse, give it general envtal state
        patch_iter->envt2 = params.environment_coarse ? 
            envt2_if_coarse
            :
            uniform(rng_r) < prob_envt2; // otherwise determine patch state individually
    }
} // end initialize_environmental_variation


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
    return((Ttot * Ttot / (
                (Ttot * Ttot + params.D * params.D) 
                ))/
                (1.0 + params.gamma * local_pop_size)
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

        assert(rare_sex != common_sex);

        // total number of individuals in the local patch
        int local_density = nmales + nfemales + 
            metapop[patch_idx].adult_care[M].size() + 
            metapop[patch_idx].adult_care[F].size() + 
            metapop[patch_idx].juvenile[M].size() + 
            metapop[patch_idx].juvenile[F].size();

        // shuffle list elements containing members of the common sex
        std::shuffle(
                metapop[patch_idx].adult_mate[common_sex].begin(), 
                metapop[patch_idx].adult_mate[common_sex].end(), 
                rng_r);

        // we don't need to shuffle the rare sex (as the common sex is already
        // shuffled, guaranteeing random mating)
        for (int rare_sex_idx = 0;
                rare_sex_idx < nrare_sex;
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
                metapop[destination_patch].juvenile[offspring_sex].push_back(Kid);
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

        // remove the first nrare_sex individuals 
        // of the common sex that have now mated
        // (remember, the common sex vector was randomized before mating
        // began, then was mated in the order they were in and now
        // we thus need to remove the first few common_sex individuals)
        
        std::vector<Individual>::iterator ind_it = 
                    metapop[patch_idx].adult_mate[common_sex].begin();

        assert(nrare_sex <= metapop[patch_idx].adult_mate[common_sex].size());

        for (
                int common_sex_idx = 0;
                common_sex_idx < nrare_sex;
                ++common_sex_idx)
        {
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
                    ind_it = metapop[patch_idx].juvenile[sex_idx].erase(ind_it);
                }
                else
                {
                    ++ind_it;
                    // update of time in current state done in juv survival function
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
            if (uniform(rng_r) < params.mu_juv[patch_i.envt2][sex_idx])
            {
                // remove individual from list
                ind_it = patch_i.juvenile[sex_idx].erase(ind_it);
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
                ind_it = patch_i.adult_mate[sex_idx].erase(ind_it);
            }
            else
            {
                // no removal? Increase element
                ++ind_it;
            }
            
        } // for individual
    } // for sex idx
} // end mate_mortality()


// mortality of both sexes
void Simulation::general_mortality()
{
    for (int patch_idx = 0; 
            patch_idx < params.npatches; ++patch_idx)
    {

        mate_mortality(metapop[patch_idx]);
        care_mortality(metapop[patch_idx]);
        juvenile_mortality(metapop[patch_idx]);

    } // for patches
} // end general_mortality()

// mortality due to the amount of care
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
                ind_it = patch_i.adult_care[sex_idx].erase(ind_it);
            }
            else
            {
                ++ind_it;
            }
        } // for individual
    } // for sex idx
} // end care_mortality()


bool Simulation::ran_out_of_caring_time(
        int const current_time, 
        double const evolved_time)
{
    int evolved_time_i = floor(evolved_time);

    // breeding period not yet over
    if (current_time < evolved_time_i)
    {
        return(false);
    }

    // this assumes that current_time > evolved_time_i

    // if there is still a bit of evolved time left
    // then draw a random number whether it is enough to care another day
    if (uniform(rng_r) < evolved_time - evolved_time_i)
    {
        return(false);
    }

    return(true);
} // ran_out_of_caring_time

// checks whether care has ended for some pairs
// if so, removes individuals back into mating pool
void Simulation::care_ends()
{
    for (int patch_idx = 0; patch_idx < params.npatches; ++patch_idx)
    {
        for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
        {
            for (
                std::vector<Individual>::iterator ind_it = 
                        metapop[patch_idx].adult_care[sex_idx].begin();
                        ind_it != metapop[patch_idx].adult_care[sex_idx].end();
            )
            {
                if (ran_out_of_caring_time(
                            ind_it->time_current_state
                            ,ind_it->T[sex_idx]))
                {
                    ind_it->time_current_state = 0;
                    metapop[patch_idx].adult_mate[sex_idx].push_back(*ind_it);

                    ind_it = metapop[patch_idx].adult_care[sex_idx].erase(ind_it);
                }
                else
                {
                    ind_it->time_current_state++;
                    ++ind_it;
                }
            } // end for ind_it
        } // end for sex_idx
    } // end for patch_idx
} // end care_ends

void Simulation::write_data_headers()
{
    data_file << "time;";

    std::string sex_identifier[2] = {"",""};

    sex_identifier[F] = "f";
    sex_identifier[M] = "m";

    for (int sex_trait_idx = 0; sex_trait_idx < 2; ++sex_trait_idx)
    {
        data_file << "mean_T_" << sex_identifier[sex_trait_idx] << ";"
            "ss_T_" << sex_identifier[sex_trait_idx] << ";"
            "n_care_" << sex_identifier[sex_trait_idx] << ";"
            "n_mate_" << sex_identifier[sex_trait_idx] << ";"
            "n_juv_" << sex_identifier[sex_trait_idx] << ";";
    }

    data_file << std::endl;
} // end void write_data_headers()

// write out summary stats
void Simulation::write_data()
{
    double mean_T[2] = {0.0,0.0};
    double ss_T[2] = {0.0,0.0};

    double x;

    int n_care[2] = {0,0};
    int n_mate[2] = {0,0};
    int n_juv[2] = {0,0};

    for (int patch_idx = 0; patch_idx < params.npatches; ++patch_idx)
    {
        for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
        {
            //  averaging over juvenile population
            for (
                std::vector<Individual>::iterator ind_it = 
                        metapop[patch_idx].juvenile[sex_idx].begin();
                        ind_it != metapop[patch_idx].juvenile[sex_idx].end();
                        ++ind_it
            )
            {
                for (int sex_trait_idx = 0; sex_trait_idx < 2; ++sex_trait_idx)
                {
                    x = ind_it->T[sex_trait_idx];
                    mean_T[sex_trait_idx] += x;
                    ss_T[sex_trait_idx] += x * x;
                }
            }
            // now averaging over care population
            for (
                std::vector<Individual>::iterator ind_it = 
                        metapop[patch_idx].adult_care[sex_idx].begin();
                        ind_it != metapop[patch_idx].adult_care[sex_idx].end();
                        ++ind_it
            )
            {
                for (int sex_trait_idx = 0; sex_trait_idx < 2; ++sex_trait_idx)
                {
                    x = ind_it->T[sex_trait_idx];
                    mean_T[sex_trait_idx] += x;
                    ss_T[sex_trait_idx] += x * x;
                }

            }
            
            // now averaging over to-mate population
            for (
                std::vector<Individual>::iterator ind_it = 
                        metapop[patch_idx].adult_mate[sex_idx].begin();
                        ind_it != metapop[patch_idx].adult_mate[sex_idx].end();
                        ++ind_it
            )
            {
                for (int sex_trait_idx = 0; sex_trait_idx < 2; ++sex_trait_idx)
                {
                    x = ind_it->T[sex_trait_idx];
                    mean_T[sex_trait_idx] += x;
                    ss_T[sex_trait_idx] += x * x;
                }
            }
                
            n_care[sex_idx] += metapop[patch_idx].
                adult_care[sex_idx].size();

            n_mate[sex_idx] += metapop[patch_idx].
                adult_mate[sex_idx].size();
            
            n_juv[sex_idx] += metapop[patch_idx].
                juvenile[sex_idx].size();
        }
    }
    
    int n_tot[2] = {0,0};

    data_file << time_step << ";";

    for (int sex_trait_idx = 0; sex_trait_idx < 2; ++sex_trait_idx)
    {
        n_tot[sex_trait_idx] = n_care[sex_trait_idx] + 
            n_mate[sex_trait_idx] + n_juv[sex_trait_idx];

        if (n_tot[sex_trait_idx] < 1)
        {
            write_parameters();
            return;
        }
    }

    for (int sex_trait_idx = 0; sex_trait_idx < 2; ++sex_trait_idx)
    {
        mean_T[sex_trait_idx] /= n_tot[M] + n_tot[F];

        ss_T[sex_trait_idx] /= n_tot[M] + n_tot[F];
        
        data_file << mean_T[sex_trait_idx] << ";"
            << ss_T[sex_trait_idx] - mean_T[sex_trait_idx] * mean_T[sex_trait_idx] << ";"
            << n_care[sex_trait_idx] << ";"
            << n_mate[sex_trait_idx] << ";"
            << n_juv[sex_trait_idx] << ";";
    }

    data_file << std::endl;
} // end void write_data

void Simulation::write_parameters()
{
    data_file 
        << std::endl 
        << std::endl
        << "init_Tf;" << params.init_Tf << std::endl
        << "init_Tm;" << params.init_Tm << std::endl
        << "environmental_change_0;" << params.environmental_change[0] << std::endl
        << "environmental_change_1;" << params.environmental_change[1] << std::endl
        << "init_prob_envt2;" << params.init_prob_envt2 << std::endl
        << "npatches;" << params.npatches << std::endl
        << "nf_per_patch_init;" << params.nf_per_patch_init << std::endl
        << "nm_per_patch_init;" << params.nm_per_patch_init << std::endl
        << "sdmu;" << params.sdmu << std::endl
        << "sigma;" << params.sigma << std::endl
        << "D;" << params.D << std::endl
        << "gamma;" << params.gamma << std::endl
        << "prob_male;" << params.prob_male << std::endl;

    std::string sex_identifier[2] = {"",""};

    sex_identifier[F] = "f";
    sex_identifier[M] = "m";

    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        data_file << "mu_mate_" << sex_identifier[sex_idx] 
                    << ";" << params.mu_mate[sex_idx] << std::endl
                    << "mu_care_" << sex_identifier[sex_idx] 
                        << ";" << params.mu_care[sex_idx] << std::endl
                    << "mu_juv_0_" << sex_identifier[sex_idx] 
                        << ";" << params.mu_juv[0][sex_idx] << std::endl
                    << "mu_juv_1_" << sex_identifier[sex_idx] 
                        << ";" << params.mu_juv[1][sex_idx] << std::endl
                    << "mutate_T_" << sex_identifier[sex_idx] 
                        << ";" << params.mutate_T[sex_idx] << std::endl
                    << "max_time_juv_" << sex_identifier[sex_idx] 
                        << ";" << params.max_time_juv[sex_idx] << std::endl
                    << "d_" << sex_identifier[sex_idx] 
                        << ";" << params.d[sex_idx] << std::endl;
    }


} // end void write_parameters()
