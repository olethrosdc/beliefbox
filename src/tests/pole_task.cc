/* -*- Mode: C++; -*- */
// $Id: test_pole_task.c,v 1.1 2006/10/23 08:33:32 olethros Exp cdimitrakakis $

/*! \file test_pole_task.cc
  \brief Test the pole balancing task.

  Use this testbench to evaluate the performance of discrete-space
  reinforcement learning algorithms in various tasks. Currently a
  maze task and a bandit task are implemented.
*/

#ifdef MAKE_MAIN

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <getopt>
#include "Distribution.h"
#include "MathFunctions.h"
#include "DiscretePolicy.h"


/// Pole-cart system, where the pole is modelled as a stiff spring.
/// Simulations show that energy remains stable in the system.
class PoleCart {
private:
    real prevL; // previous spring length
    real prev_theta;
    static const real Border;
public:
    real x; ///< position
    real u; ///< speed
    real theta; ///< angle
    real w; ///< angular velocity
    real m_cart; ///< cart mass
    real m_ball; ///< ball mass
    real n; ///< friction cart/ground
    real g;///< gravitational acceleration
    real L;///<spring length
    real L0;///<spring natural length
    real K;///< spring stiffness
    real Kd;///< spring damping
    real x_ball;///< ball x coordinate
    real y_ball;///< ball y coordinate
    real ux_ball;///< ball ux speed
    real uy_ball;///< ball uy speed
    real ground_n;///< ground friction non-linearity scaling
    PoleCart()
    {

    }
    ~PoleCart() {}
    /// Reset the cart simulation
    void Reset() {
        x = 0;
        u = 0;
        w = 0.0;
        m_cart = 1;
        m_ball = .1;
        n = 0.05;
        g = 9.82;
        L0 = 1.0;
        L = 1.0;
        prevL = L;
        K = 200.0;
        Kd = 1.0;
        theta  =  (urandom()-.5)*.01;
        x_ball  = sin(theta)*L;
        y_ball = cos(theta)*L;
        ux_ball = 0.0;
        uy_ball = 0;
        ground_n = 10.0f;
    }
    /// Simulate the cart, given external force and dt
    bool Simulate(real F, real dt)
    {
        bool failure = false;
        prev_theta = theta;
        theta=atan2(x-x_ball,y_ball);
        w = (theta-prev_theta)/dt;
        prevL=L;
        L = sqrt((x-x_ball)*(x-x_ball)+y_ball*y_ball);
        real dL=(L-prevL)/dt;
        //real spring_heat = Kd*dL*(L-prevL);
        real F_spring = K*(L0-L)*fabs(L0-L) - Kd*dL;
        real F_spring_y = -F_spring * cos(theta);
        real F_spring_x = -F_spring * sin(theta);
        real N = m_cart*g - F_spring_y;
        real friction = tanh(ground_n*u)*n*N;
        F -= friction + F_spring_x;
        real F_ballx = F_spring_x;
        real F_bally = -F_spring_y - m_ball*g;
  
        u += dt*F/m_cart;
        x += dt*u;
        if ((x<-Border)||(x>Border)) {
            x = Border * sign(x);
            u = -fabs(u)*sign(x);
            failure = true;
        }

        ux_ball += dt*F_ballx/m_ball;
        uy_ball += dt*F_bally/m_ball;
        x_ball += dt*ux_ball;
        y_ball += dt*uy_ball;
        if (y_ball<0) {
            real dE = y_ball*m_ball*g;
            y_ball = 0;
            real Kin=uy_ball*uy_ball*m_ball*.5+dE;
            if (Kin>0) {
                uy_ball=sqrt(2.0*Kin/m_ball);
            } else {
                uy_ball=0;
            }
            failure = true;
        }
        return failure;
    }
};

const real PoleCart::Border = 2.0f;

/// Discretise state
int Discretise(real x, int n_states, real dx)
{
    real ax = fabs(x);
    int state = 0;
    while (state+1 < n_states-2) {
        if (ax<dx) {
            if (x>0)
                return state+1;
            else
                return state;
        }
        dx *= 2.0;
        state+=2;
    }
    if (x>0)
        return n_states-1;
    else
        return n_states-2;
}

/// number of split states in W
#define W_STATES 4
/// number of split states in U
#define U_STATES 4
/// number of split states in X
#define X_STATES 4
/// number of split states in THETA
#define THETA_STATES 6

/// Make a new state ID from real numbers
int MakeState (real w, real u, real x, real theta)
{
    int state = 0;
    int Sw = Discretise(w, W_STATES, 0.5);
    int Su = Discretise(u, U_STATES, 0.5);
    int Sx = Discretise(x, X_STATES, 0.5);
    int Stheta = Discretise(theta, THETA_STATES, 0.05);
    state = Sw*U_STATES*X_STATES*THETA_STATES;
    state += Su*X_STATES*THETA_STATES;
    state += Sx*THETA_STATES;
    state += Stheta;
    return state;
}


static const char* const copyright_text =
    "test_pole_task $Revision: 1.1 $ \n\n\
Copyright (c) 2004-2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>\n\n\
This program is free software; you can redistribute it and/or modify\n\
it under the terms of the GNU General Public License as published by\n\
the Free Software Foundation; either version 2 of the License, or\n\
(at your option) any later version.\n";

static    const char* const help_text="\
$Id: test_pole_task.c,v 1.1 2006/10/23 08:33:32 olethros Exp cdimitrakakis $\n\
\n\
A reinforcement learning maze testbench.\n\
\n\
Usage: rl_testbench [options]\n\
Options:\n\
  --method string          type of action selection, with name:\n\
      egreedy              e-greedy action selection\n\
      softmax              softmax action selection\n\
      pursuit              pursuit action selection\n\
      confidence           confidence action selection\n\
  --sarsa                  Use Sarsa estimates\n\
  --q-learning             Use Q-learning estimates\n\
  --e-learning             Use E-learnign estimates\n\
  --n_episodes N           number of episodes\n\
  --max_episode_time x     maximum length of an episode [1000 s]\n\
  --stats_file string      file name for statistics\n\
  --max_stats N            number of runs for statistics\n\
  --zeta x                 Randomness in action selection [0.1].\n\
  --accumulating_traces    Set for accumulating instead of replacing\n\
  --dt x                   Simulation delta time [0.1]\n\
  --decision_time N        Simulation steps per decision [10]\n\
  --init_eval x            Initial evaluation [0]\n\
  --theta_factor x         Set reward to -theta*abs(angle) [0]\n\
  --gamma x                Discounting [0.9]\n\
  --lambda x               Trace decay [0.7]\n\
  --lr x                   Learning rate [0.01]\n\
  --show                   Show a run\n\
  --help                   This help\n\
  --copyright              Copyright notice\n\
\n\
What zeta means:\n\
  For egreedy, the amount of times a non-greedy action is taken.\n\
  For softmax, the smoothness of the distribution (inverse of variance).\n\
  Pursuit updates its parameters faster when this value is high.\n\
  The confidence method maintains an exponentially decaying weighted sum of\n\
  sample variances. The smaller this value, the smoother the estimate is.\n\
  Because sample variance estimators are always biased, this means that\n\
  less smooth estimates will always overestimate the variance and thus\n\
  the action selection will explore more.\n\
\n";
int main (int argc, char** argv)
{

    int max_episodes=100;
    FILE* stat_file = NULL;
    char* stat_filename = NULL;
    int n_actions = 2;
    int n_states = W_STATES*U_STATES*X_STATES*THETA_STATES;
    real alpha = 0.01;
    real gamma = 0.95f;//0.95; 0.9
    real lambda = 0.6f;//0.6; 0.7
    real zeta = 0.01;
    real max_episode_time = 1000.0f;
    char* method = NULL;
    int max_stats = 100;
    bool replacing_traces = true;
    real theta_factor = 0.0f;
    int decision_time = 10;
    real dt = 0.01f;
    real init_eval =0.0f;
    DiscretePolicy* policy;
    bool show_final_policy = false;
    LearningMethod estimate = QLearning;

    setRandomSeed(time(NULL));
    int c;
    int digit_optind = 0;
    while (1) {
        int this_option_optind = optind ? optind : 1;
        int option_index = 0;
        static struct option long_options[] = {
            {"method", 1, 0, 0}, //0
            {"n_episodes", 1, 0, 0}, //1
            {"stats_file", 1, 0, 0}, //2
            {"max_stats", 1, 0, 0}, //3
            {"accumulating_traces", 0, 0, 0}, //4
            {"dt", 1, 0, 0}, //5
            {"decision_time", 1, 0, 0}, //6
            {"init_eval", 1, 0, 0}, //7
            {"theta_factor", 1, 0, 0}, //8
            {"gamma", 1, 0, 0}, //9
            {"lambda", 1, 0, 0}, //10
            {"zeta", 1, 0, 0}, //11
            {"lr", 1, 0, 0}, //12
            {"help", 0, 0, 0}, //13
            {"copyright", 0, 0, 0}, //14
            {"show", 0, 0, 0}, //15
            {"max_episode_time", 1, 0, 0}, //16
            {"sarsa", 0, 0, 0}, //17
            {"q-learning", 0, 0, 0}, //18
            {"e-learning", 0, 0, 0}, //19
            {0, 0, 0, 0}
        };

        c = getopt_long (argc, argv, "",
                         long_options, &option_index);
        if (c == -1)
            break;

        switch (c) {
        case 0:
#if 0
            printf ("option %s (%d)", long_options[option_index].name, option_index);
            if (optarg)
                printf (" with arg %s", optarg);
            printf ("\n");
#endif
            switch (option_index) {
            case 0: method = optarg; break;
            case 1: max_episodes = atoi(optarg); break;
            case 2: stat_filename = optarg; break;
            case 3: max_stats = atoi(optarg); break;
            case 4: replacing_traces = false; break;
            case 5: dt = atof(optarg); break;
            case 6: decision_time = atoi(optarg); break;
            case 7: init_eval = atof(optarg); break;
            case 8: theta_factor = atof(optarg); break;
            case 9: gamma = atof(optarg); break;
            case 10: lambda = atof(optarg); break;
            case 11: zeta = atof(optarg); break;
            case 12: alpha = atof(optarg); break;
            case 14: 
                printf ("%s", copyright_text);
                exit(0);
                break;
            case 15: show_final_policy = true; break;
            case 16: max_episode_time = atof(optarg); break;
            case 17: estimate = Sarsa; break;
            case 18: estimate = QLearning; break;
            case 19: estimate = ELearning; break;
            case 13:
            default:
                fprintf (stderr, "%s", help_text);
                exit(0);
                break;
            }
            break;
        case '0':
        case '1':
        case '2':
            if (digit_optind != 0 && digit_optind != this_option_optind)
                printf ("digits occur in two different argv-elements.\n");
            digit_optind = this_option_optind;
            printf ("option %c\n", c);
            break;
        default:
            printf ("?? getopt returned character code 0%o ??\n", c);
            exit (-1);
        }
    }
	
    if (optind < argc) {
        printf ("non-option ARGV-elements: ");
        while (optind < argc)
            printf ("%s ", argv[optind++]);
        printf ("\n");
    }


    std::cout << "# Max episodes " << max_episodes << std::endl;

    if (stat_filename) {
        if ((stat_file = fopen(stat_filename,"w"))) {
            std::cout << "# Writing stats to: " << stat_filename << std::endl;
        } else {
            std::cerr << "# Could not open file: " << stat_filename << std::endl;
            exit(-1);
        }
    } else {
        std::cerr << "Where should I save the stats?\n";
        exit(-1);
    }

    // Pole cart system;
    PoleCart pole_cart;
    pole_cart.Reset();
	
    // paramters for policy
    //DISABLED_ASSERT(zeta>=0.0 && zeta<=1.0)(zeta);
    std::cerr << "zeta : " << zeta << std::endl;
        
    real* average_score = new real [max_episodes];
    for (int i=0; i<max_episodes; i++) {
        average_score[i] = 0.0;
    }

    
    std::cout << "# Max stats " << max_stats << std::endl;

    fprintf (stat_file, "# %s %s g:%f z:%f l:%f ep:%d st:%d dt:%f(x%d) init:%f\n",
             argv[1], argv[6], gamma, zeta, lambda, max_episodes,
             max_stats, dt, decision_time, init_eval);

    if (replacing_traces) {
        fprintf (stat_file, "# replacing traces");
    } else {
        fprintf (stat_file, "# accumulating traces");
    }


    int init_pos = 0;
    int pos = init_pos;
    
    policy = new DiscretePolicy (n_states, n_actions, alpha, gamma, lambda, false, 0.1, init_eval);
    //policy
    //std::cerr << max_stats - stats << "to go.\n";

    policy->setVarianceUpdate(Counting);

    if (replacing_traces) {
        policy->setReplacingTraces(true);
    } else {
        policy->setReplacingTraces(false);
    }
    // it seems that greedy works best with q-learning. Why?
    switch (estimate) {
    case Sarsa:
        policy->setSarsa();
        break;
    case QLearning:
        policy->setQLearning();
        break;
    case ELearning:
        policy->setELearning();
        break;
    default:
        fprintf (stderr, "Unsupported learning type\n");
        exit(-1);
        break;
    }
}

void Evaluate(char* stat_file_name, Policy& policy, int n_actions)
{
    real total_r = 0.0;
    real total_time = 0.0;
    pos = 0;
    real r = 0.0;
    FILE* small_stat_file = NULL;
    
    if (stat_filename) {
        const int FLEN = 1000;
        char small_stat_filename[FLEN];
        sprintf (small_stat_filename, "%s_sm%d", stat_filename, stats);
        small_stat_file = fopen(small_stat_filename, "w");
        if (!small_stat_file) {
            std::cerr << "# Warning: Could not open file: " << small_stat_filename << std::endl;
        }
    }
    real max_t = 0.0f;
    real local_max_t = 0.0f;
    real local_average_t = 0.0f;

    for (int episode=0; episode<max_episodes; episode++) {
        int a;
        pole_cart.Reset();
        policy->Reset();
        bool termination = false;
        real t = 0.0;
        //policy->setRandomness(zeta / ((real) (1 + episode)));
        do {
            pos = MakeState (pole_cart.w, pole_cart.u, pole_cart.x, pole_cart.theta);
            assert(pos>=0 && pos<n_states);

            a = policy->SelectAction (pos, r);
            real F = 0.0;
            switch (a) {
            case 0:
                F = 10.0; break;
            case 1:
                F = -10.0; break;
            case 2:
                F = 0.0; break;
            case 3:
                F = 20.0; break;
            case 4:
                F = -20.0; break;
            case 5:
                F = 40.0; break;
            case 6:
                F = -40.0; break;
            default:
                fprintf (stderr, "Error action\n");
            }
            r = 0.0;
            r = - theta_factor * fabs(pole_cart.theta);
            for (int i=0; i<decision_time; i++) {
                t += dt;
                if ((pole_cart.Simulate(F, dt))==true) {
                    r = -1.0;
                    termination = true;	
                }
            }
            if (t>max_episode_time) {
                termination=true;
                r = 0.0;
            }
            if (show_final_policy && episode>max_episodes - 10) {
                printf("%f %f %f %f %f %f\n",
                       pole_cart.x, pole_cart.theta, pole_cart.x_ball,
                       pole_cart.y_ball, pole_cart.u, pole_cart.w);
            }
            if (termination) {
                pos = MakeState (pole_cart.w, pole_cart.u, pole_cart.x, pole_cart.theta);
                policy->SelectAction(pos,r);
            }
				
        } while (!termination);
        local_average_t += t;
        if (local_max_t < t) {
            local_max_t = t;
            if (max_t < t) {
                max_t = t;
            }
        }
        if (episode % 100 == 0) {
            std::cerr << "# Ep: " << episode 
                      << " T_av: " << local_average_t/100.0f
                      << " T_max: " << local_max_t << std::endl;
            local_max_t = 0.0f;
            local_average_t = 0.0f;
        }
        total_r = t;
        total_time += t;
        average_score[episode] += t;
        if (small_stat_file) {
            fprintf (small_stat_file, "%f\n", t);
        }
    }
    if (small_stat_file) {
        fclose (small_stat_file);
    }
    fprintf (stderr, "%f ", total_r);

    fprintf (stderr, "%f\n", total_time/((real) max_episodes));
    delete policy;
} // for (stats)
{
    real invm=1.0/((real) max_stats);
    for (int i=0; i<max_episodes; i++) {
        average_score[i] *= invm;
        if (stat_file) {
            fprintf (stat_file, "%f\n", average_score[i]);
        }
    }
}

if (stat_file) {
    fclose(stat_file);
}
	
delete [] average_score;
std::cout <<"# Done\n";

}

	
#endif
