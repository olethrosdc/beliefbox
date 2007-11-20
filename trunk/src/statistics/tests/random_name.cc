#ifdef MAKE_MAIN

#include "BayesianMarkovChain.h"
#include "Random.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "debug.h"
#include "MarkovChain.h"

typedef struct State_ {
    char c;
} State;

#define N_STATES 27
int main (int argc, char** argv) {
    State* state;
    char* fname;
    FILE* fstream;
    int n_states = N_STATES;
    BayesianMarkovChain* chain = NULL;
    int mem_size = 1;
    int i;

    printf ("\n\
Random Name Generator v0.3: Using Bayesian Markov Chains to create random names from examples\n\n2007, Christos Dimitrakakis\n\n");

    srand48(time(NULL));

    if (argc!=4) {
        fprintf (stderr, "Usage : randomname memsize namelist prior\n\n\
Where memsize is a positive integer and namelist is the file that\n\
should be used for taking examples from.\n\n");
        return -1;
    }

    mem_size = atoi (argv[1]);
    printf ("Using a memory size of %d\n", mem_size);
    printf ("number of states: %d\n", n_states);
    fname = argv[2];
    if (!fname) {
        fprintf (stderr, "You must supply a filename\n");
    }

    float prior = atof(argv[3]);

    printf ("Creating Markov Chain\n");

    chain = new BayesianMarkovChain(n_states, mem_size, prior);

    state = (State*) malloc (n_states * sizeof(State));
    if (state==NULL) {
        delete chain;
        Serror ("malloc failed");
        return -1;
    }

  
    for (i=1; i<n_states; i++) {
        state[i].c = i-1+'a';
    }
    state[0].c = ' ';

    // clear state transition tables.
    printf ("Clearing state transitions\n");

    printf ("Opening file\n");
    fstream = fopen (fname, "r");
    if (fstream) {
        //int prev = 0;
        int curr = 0;
        int ret = 0;
        int cnt = 0;
        printf ("Processing input file...\n");
        do {
            ret = fgetc(fstream);
            cnt++;
            if (ret!=EOF) {
                //prev = curr;
                if ((ret>='A')&&(ret<='Z')) {
                    ret = ret - 'A' + 'a';
                }
                if ((ret<'a')||(ret>'z')) {
                    curr = 0;

                } else {
                    curr = 1+ret-'a';
                }

                chain->ObserveNextState(curr);
                if (curr==0) {
                    chain->Reset();
                }
            }
        } while (ret!=EOF);

        printf ("Closing file after reading %d characters\n",cnt);
        fclose (fstream);
    } else {
        fprintf (stderr, "Error, could not open file\n");
    }

  
#ifdef DISPLAY_TRANSITION_TABLES
    for (i=0; i<n_states; i++) {
        int j;
        for (j=0; j<n_states; j++) {
            printf ("%c->%c=%f\n",state[i].c,state[j].c,trans[i+j*n_states]);
        }
    }
#endif
    printf ("Probabilities of k-th order models\n");
    for (int j=0; j<mem_size; ++j) {
        printf (" %f", chain->Pr[j]);
        //chain->Pr[j]=0.0;
    }
    //chain->Pr[mem_size-1]=1.0;
    printf("\n");


    printf ("Generating names from the chain\n");
    printf ("\n-----------------------------\n");
    {
        int next = n_states-1;
        int curr = n_states-1;
        int cnt=1000;

        /* MarkovChainReset() is called at the beginning of the generation
           and also at points where a whitespace is generated. It is also
           called when nothing is generated (although this should never
           happen, as transition tables exist for even non-encountered
           states). MarkovChainReset() ensures that we only take into account
           letter order when generating words. Thus, after whitespace or end
           of generation a new word is created from scratch. */
        chain->Reset();

        for (i=0; i<50;) {
            next = chain->generate();
            //next = chain->mc[mem_size-1]->generate();
            if (next<=0) {
                chain->Reset();
            }
            if (next<0) {
                next = 0;
            }
            if ((curr!=next)||(next)) {
                printf ("%c",state[next].c);
                i++;
                cnt = 1000;
            }
            cnt--;
            curr = next;
            if (!cnt) {
                fprintf (stderr, "Waiting too long for something to be generated. Aborting\n");
                chain->ShowTransitions();
                break;
            }
        }
    }
    
    //chain->ShowTransitions();
    printf ("\n-----------------------------\n");
    printf ("Freeing memory\n");

    free (state);

    delete chain;

    printf("Exiting\n");
    return 0;
  
}

#endif
