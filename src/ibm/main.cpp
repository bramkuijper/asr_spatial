#include <string>
#include "simulation.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{
    Parameters params;

    params.mu_juv[0][F] = std::stod(argv[1]);
    params.mu_juv[1][F] = std::stod(argv[2]);
    params.mu_juv[0][M] = std::stod(argv[3]);
    params.mu_juv[1][M] = std::stod(argv[4]);
    params.sigma = std::stod(argv[5]);
    params.environmental_change[0] = std::stod(argv[6]);
    params.environmental_change[1] = std::stod(argv[7]);
    params.init_prob_envt2 = std::stod(argv[8]);
    params.base_name = argv[9];
    
    Simulation sim(params);
}
