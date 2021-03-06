#include "GaussianValueFunctionModel.h"
#include "KNNValueFunctionModel.h"
#include "Pendulum.h"
#include "Grid.h"
#include <sstream>

int Evaluate(Environment<Vector, int>& environment,
			 ValueFunctionModel<Vector, int>& vfm,
			 int grid_size,
			 int evaluation_grid_size,
			 bool use_kernel,
			 RBFBasisSet& kernel,
			 int n_iter,
			 real gamma,
			 int T)
{
	for (int k=0; k<n_iter; k++) {
		real U = 0;
		real discount = 1 ;
		environment.Reset();
		environment.setState(urandom(environment.StateLowerBound(),
									 environment.StateUpperBound()));
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
		//printf ("%f -- ", U);
		//state.print(stdout);
		if (use_kernel) { 
			kernel.Evaluate(state);
			vfm.AddReturnSample(kernel.F(), action, U);
		} else {
			vfm.AddReturnSample(state, action, U);
		}
	}
	vfm.CalculateValues();
	std::stringstream fname;
	fname  << environment.Name() << vfm.Name();
	if (use_kernel) {
		fname << grid_size << "g";
	} else {
		fname << "raw";
	}
	fname << ".values";
	
	std::string fname_string = fname.str();
    FILE* outfile = fopen(fname_string.c_str(), "w");
    if (outfile) {
        EvenGrid evaluation_grid(environment.StateLowerBound(),
                                 environment.StateUpperBound(),
                                 evaluation_grid_size);
        for (int i=0; i<evaluation_grid.getNIntervals(); ++i) {
            Vector state = evaluation_grid.getCenter(i);
			state.printf(outfile);
			if (use_kernel) {
				kernel.Evaluate(state);
				fprintf(outfile, " %f", vfm.getValue(kernel.F()));
			} else {
				fprintf(outfile, " %f", vfm.getValue(state));
			}
			fprintf(outfile, "\n");
        }
        fclose(outfile);
    } else {
        Serror("Failed to write to file\n");
    }
	return 0;
}

int main (int argc, char* argv[])
{
	Pendulum pendulum;
	Environment<Vector, int>& environment = pendulum;
		
	int T = 100;
	int grid_size = 8;
	int evaluation_grid_size = 32;

	int n_iter = 10000;
	real gamma = 0.95;
	bool use_kernel = false;


	logmsg("Usage: test_gvfm use_kernel [{0,1}] grid_size\n");
	
	use_kernel = (atoi(argv[1])!=0);
	grid_size = atoi(argv[2]);
	
    logmsg("Generating grid\n");
    EvenGrid grid(environment.StateLowerBound(),
                  environment.StateUpperBound(),
                  grid_size);

    logmsg("Creating kernel\n");
    // use an RBF basis for the kernel fromthe grid
    RBFBasisSet kernel(grid, 1.0);
	int n_features = kernel.size();

	if (!use_kernel) {
		n_features = environment.getNStates();
	}
	GaussianValueFunctionModel gvfm(n_features,
									environment.getNActions());

	Evaluate(environment,
			 gvfm,
			 grid_size,
			 evaluation_grid_size,
			 use_kernel,
			 kernel,
			 n_iter,
			 gamma,
			 T);
			 
	KNNValueFunctionModel kvfm(n_features,
							   environment.getNActions(),
							   3);


	Evaluate(environment,
			 kvfm,
			 grid_size,
			 evaluation_grid_size,
			 use_kernel,
			 kernel,
			 n_iter,
			 gamma,
			 T);

	return 0;
	//ValueFunctionModel<Vector, int>* vfm = &gvfm;
}

