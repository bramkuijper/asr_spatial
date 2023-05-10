#include "individual.hpp"
#include "parameters.hpp"

Individual::Individual(
        double const init_Tf
        ,double const init_Tm
        ) :
    T{init_Tf,init_Tm}
{}

// copy constructor
Individual::Individual(Individual const &other) :
    T{other.T[F],other.T[M]}
    ,time_current_state{other.time_current_state}
{}

// birth constructor
Individual::Individual(Individual const &mom
        ,Individual const &dad
        ,Parameters const &params
        ,std::mt19937 &rng_r)  :
    time_current_state{0}
{
    std::uniform_real_distribution <double> uniform{0.0,1.0};

    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        T[sex_idx] = uniform(rng_r) < 0.5 ? dad.T[sex_idx] : mom.T[sex_idx];

        if (uniform(rng_r) < params.mutate_T[sex_idx])
        {
//            std::normal_distribution <double> mutational_effect_size{0.0,params.sdmu};

            T[sex_idx] += uniform(rng_r) < 0.5 ? 1.0 : -1.0;

            if (T[sex_idx] <= 0.0)
            {
                T[sex_idx] = 0.0;
            }
        }
    }
} // end birth constructor


void Individual::operator=(Individual const &other)
{
    T[F] = other.T[F];
    T[M] = other.T[M];
    time_current_state = other.time_current_state;
}
