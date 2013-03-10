/**
* Definition of the Blackjack environment
*
*/



///           Adding assert statements


#include "Blackjack.h"
#include "Distribution.h"
#include "SingularDistribution.h"
#include "NormalDistribution.h"

Blackjack::Blackjack(real win_value_, real draw_value_, real loss_value_):
                                    win_value(win_value_),draw_value(draw_value_),
                                    loss_value(loss_value_)
{
      //ctor

      /// Number of actions
      n_actions=2;
      n_states=(21-12+1)*(10-1+1)*2+1;
      terminal_state=n_states-1;
     
      mdp=getMDP();
       
      RandomNumberGenerator *rng;
      MersenneTwisterRNG mersenne_twister;
      rng = (RandomNumberGenerator*) &mersenne_twister;
      //rng->seed();
      rng->manualSeed(380731572);
      
      srand48(34987235);
    srand(34987235);
      
      Reset();
      
}


void Blackjack::getTriplet(int& pc_,int& dc_,int& isUsable_)
{
      isUsable_=0;
      int state_=state;
      if(state_>100){
            state_-=100;
            isUsable_=1;
      }
      pc_=state_%10;
      dc_=(state_-pc_)/10;
      pc_+=12;
      dc_+=1;


}

real Blackjack::calculReward() const
{
      if(pc>21)
      return loss_value;
      if(isNatural && pc==21)
      {
            if(dc==21)
            return draw_value;
            else
            return win_value;

      }

      /// The dealer is playing

      int dhv=dc;
      while(dhv<17)
      {

            /*RandomNumberGenerator *rng;
            MersenneTwisterRNG mersenne_twister;
            rng = (RandomNumberGenerator*) &mersenne_twister;
            rng->seed();*/

            //int j=std::min((int) rng->discrete_uniform(13)+1,10);
            int j=std::min(rand()%13+1,10);
            if(j==1 && dhv<12)
                  dhv+=10;
            dhv+=j;
            //deck_position++;
      }
      //dc=dhv;

      if(dhv>21)
      return win_value;
      if(dhv>pc)
      return loss_value;
      else if(dhv==pc)
      return draw_value;
      return win_value;
}

DiscreteMDP* Blackjack::getMDP() const
{
      /// Initialization of the MDP
      //DiscreteMDP* mdp;
      /*uint a=0,b=0;
      real c;
      NormalDistribution* reward_dist
                = new NormalDistribution(urandom(), 1.0);
            mdp->setRewardDistribution(a, b, reward_dist);*/
      //mdp->setFixedReward(a,b,c);
     DiscreteMDP* mdp = new DiscreteMDP (n_states, n_actions, NULL);

      /// Initialization of the reward

      for(uint i=0;i<n_states;i++)
      {
            for(int k=0;k<2;k++)
            mdp->setFixedReward(i,k,0);
      }

      /// Setting the reward and the MDP  state transition probability kernel

      real pr=1.0/10.0;

      for(uint pci=12;pci<22;pci++){
            for(uint dci=1;dci<11;dci++){
                  for(uint isUsablei=0;isUsablei<2;isUsablei++){
                        int s=getState(pci,dci,isUsablei);
                        if(s==terminal_state){
                            mdp->setTransitionProbability(s,0,terminal_state,1.0);
                            mdp->setTransitionProbability(s,1,terminal_state,1.0);
                            continue;
                        }
                        for(uint cardi=1;cardi<11;cardi++){

                              int s_next=getState(pci+cardi,dci,isUsablei);
                              mdp->setTransitionProbability(s,0,s_next,pr);
                              mdp->setTransitionProbability(s,1,terminal_state,1.0);
                              if(s_next==terminal_state){
                                    mdp->setFixedReward(s,0,1);

                              }
                              mdp->setFixedReward(s,1,1);


                        }


                  }

            }
      }





      return mdp;
}

/**
* Performing  action from state to a new state
*/

bool Blackjack::Act(const int& action)
{

      state=mdp->generateState(state,action);
      if(action==0)
      getTriplet(pc,dc,isUsable);

      if(state==terminal_state)
      {
            dc+=second_value;
            reward=calculReward();

            return false;
      }
      else
      {
            reward=0.0;
            isNatural=false;
            return true;
      }
}

/**
* Getting the state of the environment
*/


int Blackjack::getState(int pc_,int dc_,int isUsable_) const
{
      int s=(pc_-12)+(dc_-1)*10;
      if(isUsable_==1)
      {
      	s+=100;
      	s=std::min(s,terminal_state);
      }	
      
      else
      {
       	if(s>=100)
       	s=terminal_state;
      }
      
      return s;
} 

void Blackjack::Reset()
{
      isNatural=true;
      
      //dc=std::min((int) rng->discrete_uniform(13)+1,10);
      dc=std::min(rand()%13+1,10);
      
      
      //second_value=std::min((int) rng->discrete_uniform(13)+1,10);
      second_value=std::min(rand()%13+1,10);
      //state=(int) rng->discrete_uniform(n_states-1);
      state=rand()%(n_states-1);
      
      mdp->setState(state);
      getTriplet(pc,dc,isUsable);
      //std::cout<<" Pc= "<<pc<<" dc= "<<dc<<" isUsable "<<isUsable<<std::endl; 
      reward=0.0;
}

Blackjack::~Blackjack()
{
      //dtor
      delete mdp;
}
