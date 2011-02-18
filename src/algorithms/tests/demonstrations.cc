#include <cstdio>

#include "Demonstrations.h"

int main()
{
    Demonstrations<int, int> demo;

    for (int episode = 0; episode<5; episode++) {
        if (episode) {
            demo.NewEpisode();
        }
        for (int i=0; i<10; ++i) {
            demo.Observe(i, episode);
        }
    }
}

