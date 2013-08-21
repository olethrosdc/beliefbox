#include <iostream>
#include <bitset>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <map>
#include <sstream>
#include <fstream>
#include <cmath>
    
    
    
    
#include "PolicyEvaluation.h"
#include "PolicyIteration.h"
#include "ValueIteration.h"
#include "RandomMDP.h"
#include "Gridworld.h"
#include "DiscreteChain.h"
#include "OptimisticTask.h"
#include "OneDMaze.h"
#include "InventoryManagement.h"
#include "DiscretePolicy.h"
#include "Environment.h"
#include "ExplorationPolicy.h"
#include "Sarsa.h"
#include "QLearning.h"
#include "QLearningDirichlet.h"
#include "ModelBasedRL.h"
#include "SampleBasedRL.h"
#include "ModelCollectionRL.h"
#include "ContextBanditGaussian.h"
#include "ContextBandit.h"
#include "DiscreteMDPCollection.h"
#include "ContextBanditCollection.h"
#include "RandomNumberFile.h"
#include "MersenneTwister.h"
#include "MountainCar.h"
#include "DiscretisedEnvironment.h"
#include "HQLearning.h"
#include "MWAL.h"
#include "PolicyBelief.h"
#include "RewardPolicyBelief.h"
#include "PolicyRewardBelief.h"
#include "EasyClock.h"
#include "Blackjack.h"
    
    
//using namespace std;

#define N_COMP 4
    
    
struct Statistics
{
	
	int n_comp;
	std::vector<double> reward;
	std::vector<int> n_games;
	std::vector<std::string> name;
};
    
    
enum SamplerType {
    INVALID = 0x0,
    METROPOLIS,
    MONTE_CARLO
};
    
struct Sampler
{
    SamplerType type;
    int n_chains;
    int n_chain_samples;
};
    
Statistics EvaluateAlgorithm (int episode_steps,
                              int n_episodes,
                              uint n_steps,
                              OnlineAlgorithm<int,int>* algorithm,
                              real epsilon,
                              DiscreteEnvironment* environment,
                              real gamma,
                              int iterations,
                              Sampler sampler,
                              real accuracy);
    
    
const uint n_steps_optimal=100000;
const int n_episodes_optimal=1000;
const int n_episodes_play=1000;
int main()
{
    
	int n_actions = 2;
	int n_states = 201;
	real gamma = 0.9;
	real lambda = 0.9;
	int iterations = 1000;
	real alpha = 0.1;
	real epsilon = 0.01;
	uint n_runs = 1;
	uint n_episodes = 100;
	uint n_steps = 100000;
	uint episode_steps = 100;
	real randomness = 0.001;
    
    Sampler sampler = {MONTE_CARLO, 2, 10};
    real accuracy = 1e-2;
    
    assert (n_states > 0);
	assert (n_actions > 0);
	assert (gamma >= 0 && gamma <= 1);
	assert (lambda >= 0 && lambda <= 1);
	assert (randomness >= 0 && randomness <= 1);
	assert (n_runs > 0);
	assert (n_episodes > 0);
	assert (n_steps > 0);
    //  assert (grid_size > 0);
    if (sampler.type == INVALID) {
        Serror ("invalid sampler specified\n");
        exit(-1);
    }
	srand48(34987235);
	srand(34987235);
	//setRandomSeed(34987235);
    
    
	std::cout << "Starting test program" << std::endl;
    
	std::cout << "Starting evaluation" << std::endl;
	
	Statistics stat;
	stat.reward.resize(N_COMP,0);
	stat.n_games.resize(N_COMP,0);
	stat.n_comp=N_COMP;
	stat.name.resize(N_COMP);
	stat.name[0]="Optimal";
	stat.name[2]=" PRB ";
	stat.name[1]=" RPB";
	stat.name[3]=" Imit";
	
    
    
	for (uint run=0; run<n_runs; ++run) {
    
        DiscreteEnvironment* environment = NULL;
	      
        environment=new Blackjack();
	      
    
        n_states = environment->getNStates();
        n_actions = environment->getNActions();
	      
    
        VFExplorationPolicy* exploration_policy = NULL;
        exploration_policy = new EpsilonGreedy(n_actions, epsilon);
    
        OnlineAlgorithm<int, int>* algorithm = NULL;
        algorithm = new QLearning(n_states,
                                  n_actions,
                                  gamma,
                                  lambda,
                                  alpha,
                                  exploration_policy,
                                  1.0);
					  
    
		Statistics run_statistics = EvaluateAlgorithm(episode_steps,
                                                      n_episodes,
                                                      n_steps,
                                                      algorithm,
                                                      epsilon,
                                                      environment,
                                                      gamma,
                                                      iterations,
                                                      sampler,
                                                      accuracy);
													     
													     
	    for(int i=0;i<stat.n_comp;i++){
	    
            stat.reward[i]+=run_statistics.reward[i];
            stat.n_games[i]+=run_statistics.n_games[i];
	    
	    }
													      
													      
													      
    }
													      
    std::cout<<"Algorithm\t"<<" Number of points\t"<<" Number of games\t"<<std::endl;												  
    for(int i=0;i<stat.n_comp;i++){
	    
	    stat.reward[i]/=n_runs;
	    stat.n_games[i]/=n_runs;
	    std::cout<<stat.name[i]<<"\t\t"<<stat.reward[i]<<"\t\t"<<stat.n_games[i]<<std::endl;
	    
    }
	    
	    
    
	std::cout << "# Done" << std::endl;
    
    
    
	return 0;
}
    
    
    
/*** Evaluate an algorithm
    
     episode_steps: maximum number of steps per episode. If negative, then ignore
     n_steps: maximun number of total steps. If negative, then ignore.
     n_episodes: maximum number of episodes. Cannot be negative.
*/
    
Statistics EvaluateAlgorithm (int episode_steps,
                              int n_episodes,
                              uint n_steps,
                              OnlineAlgorithm<int, int>* algorithm,
                              real epsilon,
                              DiscreteEnvironment* environment,
                              real gamma,
                              int iterations,
                              Sampler sampler,
                              real accuracy)
{
    assert(accuracy >= 0);
	std:: cout << "evaluating..." << environment->Name() << std::endl;
    
	std:: cout << "(value iteration)" << std::endl;
    
    
	Statistics statistics;
	
	statistics.reward.resize(N_COMP,0);
	statistics.n_games.resize(N_COMP,0);
	statistics.n_comp=N_COMP;
	statistics.name.resize(N_COMP);
	statistics.name[0]="Optimal";
	statistics.name[2]=" PRB ";
	statistics.name[3]=" Imit ";
	statistics.name[1]=" RPB";
       
    
    
	int current_time = 0;
	environment->Reset();
    
	std:: cout << "(running)" << std::endl;
    int episode = -1;
    bool action_ok = false;
	    
	    
	Demonstrations<int, int> demonstrations;
	
	
	/// Observing the agent demonstrations
    
	for (uint step = 0; step < n_steps; ++step) {
        if (!action_ok) {
            if (episode > 0) {
                demonstrations.NewEpisode();
		    
                int state1 = environment->getState();
		    
                real reward1 = environment->getReward();
                //std::cout<<" Rew: "<<reward1<<std::endl;
                algorithm->Act(reward1, state1);
            }
            episode++;
            if (n_episodes >= 0 && episode >= n_episodes) {
                //fprintf (stderr, "Breaking after %d episodes,  %d steps\n",
                //episode, step);
                break;
            } else {
                environment->Reset();
                if (algorithm) {
                    algorithm->Reset();
                }
                action_ok = true;
                current_time = 0;
            }
        }
    
        int state = environment->getState();
        real reward = environment->getReward();
        //std::cout<<" Rew: "<<reward<<std::endl;
    
        int action=0;
	    if (algorithm) {
            action = algorithm->Act(reward, state);
	    
		    
		    action_ok = environment->Act(action);
		    current_time++;
            if (step > n_steps / 2 || episode > n_episodes / 2) {
                demonstrations.Observe(state, action);
            }
    
        }
	}
    
    
    
	
    double start_time;
    double end_time;
    //fprintf (stderr, "Trying to guess policy!\n");
    DiscreteMDP* mdp = environment->getMDP();
    if (!mdp->Check()) {
		Serror("MDP model is nonsense\n");
		mdp->ShowModel();
    }
    
    int n_states = mdp->getNStates();
    int n_actions = mdp->getNActions();
	    
	    
	    
    // ----- Approximate Optimal Policy  by QLearning --------//
	    
    algorithm->Reset();
    action_ok=false;
    episode=-1;
    double point1=0;
    double num1=0;
	    
    for (uint step = 0; step < n_steps_optimal; ++step) {
        if (!action_ok) {
            if (episode > 0) {
		    
		    
                int state1 = environment->getState();
		    
                real reward1 = environment->getReward();
                //std::cout<<" Rew: "<<reward1<<std::endl;
                algorithm->Act(reward1, state1);
                if(episode>n_episodes_optimal-n_episodes_play){
		      
                    num1++;
                    point1+=reward1;
		      
                }
            }
            episode++;
			    
            if (n_episodes_optimal >= 0 && episode >= n_episodes_optimal) {
				    
                break;
            } else {
				    
                environment->Reset();
		   
                action_ok = true;
                current_time = 0;
            }
        }
    
        int state = environment->getState();
        real reward = environment->getReward();
        //std::cout<<" Rew: "<<reward<<std::endl;
    
    
        int action=0;
	    if (algorithm) {
            action = algorithm->Act(reward, state);
	    
		    action_ok = environment->Act(action);
		    current_time++;
    
        }
	}
	
	statistics.reward[0]+=point1;
	statistics.n_games[0]+=num1;
	
	std::cout<<point1<<" Points after "<<num1<<" Plays"<<std::endl;
	
	    
    
    // -------- RPB : Reward | Policy Belief -------- //
    real expected_optimality = 1.0;
    DirichletDistribution dirichlet(n_states * n_actions);
    start_time = GetCPU();
    RewardPolicyBelief reward_policy_belief (expected_optimality,
                                             gamma,
                                             *mdp,
                                             dirichlet,
                                             sampler.n_chain_samples);
    
    reward_policy_belief.setAccuracy(accuracy);
    
    DiscretePolicy* rpb_policy = reward_policy_belief.CalculatePosterior(demonstrations);
    end_time = GetCPU();
    printf("%f # T_RPB\n", end_time - start_time);
    printf ("Posterior policy:\n---------------\n");
    rpb_policy->Show();
	    
	    
    std::cout<<"Statistics  "<<std::endl;
	    
    double point=0;
    int num=0;
	    
    episode = -1;
    action_ok = false;
	    
	    
	    
    for (uint step = 0; step < n_steps; ++step) {
        if (!action_ok) {
            if (episode > 0) {
		  
		  
		    
                real reward1 = environment->getReward();
                point+=reward1;
                num++;
		    
            }
            episode++;
            if (n_episodes_play >= 0 && episode >= n_episodes_play) {
				    
                break;
            } else {
				    
                environment->Reset();
		  
                action_ok = true;
                current_time = 0;
            }
        }
    
        int state = environment->getState();
        real reward = environment->getReward();
        //std::cout<<" Rew: "<<reward<<std::endl;
    
		    
        point+=reward;
		    
    
        int action;
		    
        rpb_policy->Observe(reward, state);
        action = rpb_policy->SelectAction();
		    
	   
        action_ok = environment->Act(action);
        current_time++;
    
	}
	
    statistics.reward[1]+=point;
	statistics.n_games[1]+=num;
	
	std::cout<<point<<" Points after "<<num<<" Plays"<<std::endl;
	    
	    
	    
	    
	    
	    
    /* Calculating the Value function
       PolicyEvaluation rpb_evaluator(rpb_policy, mdp, gamma);
       rpb_evaluator.ComputeStateValues(accuracy);
       for (int i=0; i<n_states; ++i) {
       printf ("%f ", rpb_evaluator.getValue(i));
       }*/
    printf ("# V_RPB\n");
    delete rpb_policy;
    

    // ----- PRB: Policy | Reward Belief ------ //
    start_time = GetCPU();
    PolicyRewardBelief prb(1.0, gamma, *mdp);
    
    if (sampler.type == METROPOLIS) {
        logmsg("Metropolis sampler: %d %d\n", sampler.n_chain_samples, sampler.n_chains);
        prb.MHSampler(demonstrations,
                      sampler.n_chain_samples,
                      sampler.n_chains);
    } else if (sampler.type == MONTE_CARLO) {
        logmsg("Monte Carlo sampler: %d %d\n", sampler.n_chain_samples, sampler.n_chains);
        prb.MonteCarloSampler(demonstrations, sampler.n_chain_samples);
    }
    
    DiscretePolicy* prb_policy = prb.getPolicy();
    end_time = GetCPU();
    printf("%f # T_PRB\n", end_time - start_time);
    printf ("MH sampler policy:\n---------------\n");
    prb_policy->Show();
    
    /*PolicyEvaluation prb_evaluator(prb_policy, mdp, gamma);
      prb_evaluator.ComputeStateValues(accuracy);
      for (int i=0; i<n_states; ++i) {
      printf ("%f ", prb_evaluator.getValue(i));
      }*/
    num=0;
    point=0;
	    
    episode = -1;
    action_ok = false;
	    
	    
    for (uint step = 0; step < n_steps; ++step) {
        if (!action_ok) {
            if (episode > 0) {
                real reward1 = environment->getReward();
                point+=reward1;
                num++;
            }
            episode++;
            if (n_episodes_play >= 0 && episode >= n_episodes_play) {
                break;
            } else {
                environment->Reset();
                action_ok = true;
                current_time = 0;
            }
        }
    
        int state = environment->getState();
        real reward = environment->getReward();
        point+=reward;
        int action;
        prb_policy->Observe(reward, state);
        action = prb_policy->SelectAction();
        action_ok = environment->Act(action);
        current_time++;
	}
	
    statistics.reward[2]+=point;
	statistics.n_games[2]+=num;
	std::cout<<point<<" Points after "<<num<<" Plays "<<std::endl;
    printf ("# V_PRB\n");
    delete prb_policy;

    // ----- PRB: Policy | Reward Belief ------ //
    start_time = GetCPU();
    FixedDiscretePolicy imitating_policy(n_states, n_actions, demonstrations);
    end_time = GetCPU();
    printf("%f # T_IMIT\n", end_time - start_time);
    printf ("IMIT sampler policy:\n---------------\n");
    imitating_policy.Show();
    
    num=0;
    point=0;
    episode = -1;
    action_ok = false;
	    
    for (uint step = 0; step < n_steps; ++step) {
        if (!action_ok) {
            if (episode > 0) {
                real reward1 = environment->getReward();
                point+=reward1;
                num++;
            }
            episode++;
            if (n_episodes_play >= 0 && episode >= n_episodes_play) {
                break;
            } else {
                environment->Reset();
                action_ok = true;
                current_time = 0;
            }
        }
    
        int state = environment->getState();
        real reward = environment->getReward();
        point+=reward;
        int action;
        imitating_policy.Observe(reward, state);
        action = imitating_policy.SelectAction();
        action_ok = environment->Act(action);
        current_time++;
	}
	
    statistics.reward[3]+=point;
	statistics.n_games[3]+=num;
	std::cout<<point<<" Points after "<<num<<" Plays "<<std::endl;
    printf ("# V_IMIT\n");
    
    delete mdp;
    
	
    
	return statistics;
    
}
