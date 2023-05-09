#include <string>
#include "simulation.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{
    Parameters params;

    params.mu_juv[F] = std::stod(argv[1]);
    params.mu_juv[M] = std::stod(argv[2]);
    params.base_name = argv[3];
    
    Simulation sim(params);
}