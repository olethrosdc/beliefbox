#include "POMDPWrapper.h"
#include "POMDPGridworld.h"
#include "MersenneTwister.h"
#include "RandomSourceRNG.h"

/// Compilation test.
int main(void)
{
    
    RandomSourceRNG rng(false);
    POMDPGridworld gridworld(&rng, "/home/olethros/projects/beliefbox/data/meze1", 2, 0.0, -1.0, 1.0, 0.0);
    POMDPWrapper<int, int, POMDPGridworld> wrapper(gridworld);

    return 0;
}
