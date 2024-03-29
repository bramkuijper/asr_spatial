#include "individual.hpp"
#include "parameters.hpp"

Individual::Individual(
        double const init_Tf
        ,double const init_Tm
        ) :
    T{init_Tf,init_Tm}
    ,Tb{0,0}
    ,phen{init_Tf}
{}

// copy constructor
Individual::Individual(Individual const &other) :
    T{other.T[F],other.T[M]}
    ,Tb{other.Tb[F],other.Tb[M]}
    ,phen{other.phen}
    ,time_current_state{other.time_current_state}
    ,natal_sr{other.natal_sr}
    ,natal_envt{other.natal_envt}
{}

// birth constructor
Individual::Individual(Individual const &mom
        ,Individual const &dad
        ,Parameters const &params
        ,double const natal_sr
        ,bool const natal_envt
        ,std::mt19937 &rng_r)  :
    time_current_state{0}
    ,phen{0.0}
    ,natal_sr{natal_sr}
    ,natal_envt{natal_envt}
{
    std::uniform_real_distribution <double> uniform{0.0,1.0};
    std::normal_distribution <double> norm{0.0,params.sdmu_Tb};

    for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
    {
        T[sex_idx] = uniform(rng_r) < 0.5 ? dad.T[sex_idx] : mom.T[sex_idx];
        Tb[sex_idx] = uniform(rng_r) < 0.5 ? dad.Tb[sex_idx] : mom.Tb[sex_idx];

        if (uniform(rng_r) < params.mutate_T[sex_idx])
        {
//            std::normal_distribution <double> mutational_effect_size{0.0,params.sdmu};

            T[sex_idx] += uniform(rng_r) < 0.5 ? 1.0 : -1.0;

            if (T[sex_idx] <= 0.0)
            {
                T[sex_idx] = 0.0;
            }
        }

        if (uniform(rng_r) < params.mutate_Tb)
        {
            Tb[sex_idx] += norm(rng_r);
        }
    }
} // end birth constructor


void Individual::operator=(Individual const &other)
{
    T[F] = other.T[F];
    T[M] = other.T[M];
    Tb[F] = other.Tb[F];
    Tb[M] = other.Tb[M];
    phen = other.phen;
    time_current_state = other.time_current_state;
}
