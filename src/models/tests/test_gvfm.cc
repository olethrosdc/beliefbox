#include "GaussianValueFunctionModel.h"
#include "Pendulum.h"
#include "Grid.h"

int main (void)
{
	Pendulum pendulum;
	Environment<Vector, int>& environment = pendulum;
	GaussianValueFunctionModel gvfm(environment.getNStates(),
									environment.getNActions());
		
	int T = 100;
	int evaluation_grid_size = 32;

	int n_iter = 10000;
	real gamma = 0.95;
	
	for (int k=0; k<n_iter; k++) {
		real U = 0;
		real discount = 1 ;
		environment.Reset();
		bool action_ok = true;
		Vector state = environment.getState();
		int action = urandom(0, environment.getNActions());
		action_ok = environment.Act(action);
		// get the return of a uniform policy
		for (int t=0; t<T && action_ok; t++) {
			real reward = environment.getReward();
			U += discount * reward;
			discount *= gamma;
			action_ok = environment.Act(urandom(0, environment.getNActions()));
			//printf ("%d %f\n",t, reward);
		}
		real reward = environment.getReward();
		U += discount * reward;
		//printf ("END %f | %f, %f\n", reward, discount, U);
		printf ("%f -- ", U);
		state.print(stdout);
		gvfm.AddReturnSample(state, action, -U);
	}
	gvfm.CalculateValues();
    FILE* outfile = fopen("Pendulum-gvfm.values", "w");
    if (outfile) {
        EvenGrid evaluation_grid(environment.StateLowerBound(),
                                 environment.StateUpperBound(),
                                 evaluation_grid_size);
        for (int i=0; i<evaluation_grid.getNIntervals(); ++i) {
            Vector state = evaluation_grid.getCenter(i);
            fprintf(outfile, "%f ", gvfm.getValue(state));
            state.print(outfile);
        }
        fclose(outfile);
    } else {
        Serror("Failed to write to file\n");
    }
	
	

	return 0;
}
