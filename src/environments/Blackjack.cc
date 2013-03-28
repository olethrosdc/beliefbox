/**
* Definition of the Blackjack environment
*
*/





#include "Blackjack.h"
#include "Distribution.h"
#include "SingularDistribution.h"
#include "NormalDistribution.h"
#include <time.h>
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
      srand48( time(NULL) );
      srand ( time(NULL) );

      setRandomSeed( time(NULL) );

      MersenneTwisterRNG *mersenne_twister=new MersenneTwisterRNG();
      rng = (RandomNumberGenerator*) mersenne_twister;
      rng->manualSeed(time(NULL));
      Reset();

}


void Blackjack::getTriplet(int& pc_,int& dc_,int& isUsable_)
{
      isUsable_=0;
      int state_=state;
      if(state_>=100){
            state_-=100;
            isUsable_=1;
      }
      pc_=state_/10;
     dc_=state_%10;
      pc_+=12;
      dc_+=1;


}

real Blackjack::calculReward() const
{
	/// Player bust
	if(pc>21)
	return loss_value;
      
     /// Player got a natural
      if(isNatural && pc==21)
      {
            if(dc==21)
            return draw_value;
            else
            return win_value+0.5;
            

      }
      
      /// Including natural win for the dealer ???
      /*
      if(dc==21)
      return loss_value;
      
      */
      
      /// The dealer is playing

      int dhv=dc;
      while(dhv<17)
      {
            int j=std::min((int) rng->discrete_uniform(13)+1,10);
           // int j=std::min(rand()%13+1,10);
            dhv+=j;
            if(j==1 && dhv<12)
                  dhv+=10;
            
      }
     
	/// Dealer bust
      if(dhv>21)
      	return win_value;
      	
      	/// Tie
      
      if(dhv==pc)
     	 return draw_value;
      
      /// Player win
       if(dhv<pc)
      	return win_value;
      
      /// Player loose
      if(dhv>pc)
      	return loss_value;
      std::cout<<" Error should not come here ";
      exit(-1);
     return 0;
}

DiscreteMDP* Blackjack::getMDP() const
{
      /// Initialization of the MDP

     DiscreteMDP* mdp = new DiscreteMDP (n_states, n_actions, NULL);

      /// Initialization of the reward

      for(uint i=0;i<n_states;i++)
      {
            for(int k=0;k<2;k++)
            mdp->setFixedReward(i,k,0);
      }

      /// Setting the reward and the MDP  state transition probability kernel

      real pr=1.0/13.0;

      mdp->setTransitionProbability(terminal_state,0,terminal_state,1.0);
      mdp->setTransitionProbability(terminal_state,1,terminal_state,1.0);


      for(uint pci=12;pci<22;pci++){
            for(uint dci=1;dci<11;dci++){
                  for(uint isUsablei=0;isUsablei<2;isUsablei++){
                        int s=getState(pci,dci,isUsablei);

                        for(uint cardi=1;cardi<11;cardi++){
                        if(cardi==10)
                        pr=4.0/13.0;
                        else
                        pr=1.0/13.0;

                              int s_next=getState(pci+cardi,dci,isUsablei);
                              double tmp1=mdp->getTransitionProbability(s,0,s_next);
                              mdp->setTransitionProbability(s,0,s_next,pr+tmp1);
                              mdp->setTransitionProbability(s,1,terminal_state,1.0);

                        }
                        mdp->setFixedReward(s,1,0);
                        mdp->setFixedReward(s,0,0);


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
	
	int state1;
	state1=state;
	
       /// Generating the resulting state corresponding to the effect of action on the current state

     //state=mdp->generateState(state1,action);
     
     
     
    /// Hitting
     if(action==0)
     {
     	int tmp1=std::min((int) rng->discrete_uniform(13)+1,10);
     	//int tmp1=std::min(rand()%13+1,10);
     	state=getState(pc+tmp1,dc,isUsable);
     	pc+=tmp1;
     }
     
     /// Standing
     if(action==1)
     state=terminal_state;
      

      /// If action is to stand or state is terminal_state computing the reward for the player
      if(state==terminal_state || action==1)
      {
            //dc+=second_value;
            
            /// Updating dealer total point
            dc=first_value+second_value;
            reward=calculReward();

            return false;
      }
      else
      {
            reward=0.0;
            return true;
      }
}

/**
* Getting the state of the environment
*/


int Blackjack::getState(int pc_,int dc_,int isUsable_) const
{
    assert(pc_>=12 && dc_>=1);
      int s=(pc_-12)*10+(dc_-1);
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
     /// The player has not a Blackjack or Natural
     isNatural=false;
     
     /// Dealing the first card of the Dealer
     
      first_value=std::min((int) rng->discrete_uniform(13)+1,10);
      
      /// Dealing the second card of the dealer
      second_value=std::min((int) rng->discrete_uniform(13)+1,10);
      //second_value=std::min(rand()%13+1,10);
      
      /// Dealing cards to the player until he has 12 or more
      int state_ini=0;
      int isU=0;
      int num=0;
      while(state_ini<12)
      {
      	int tmp=std::min((int) rng->discrete_uniform(13)+1,10);
      	//int tmp=std::min(rand()%13+1,10);
      	
      	
      	state_ini+=tmp;
      	num++;
      	
      	if(tmp==1 && state_ini<12)
      	{
      		state_ini+=10;
      		isU=1;
      	}
      	
      }
      
      /// If the player has only a Face card and an Ace then he has a natural
      
      if(isU==1 && state_ini==21 && num==2)
      isNatural=true;
     
      
      
      /// If the dealer first two cards contains an Ace then Add 10 to his second_value
      
      if( (first_value==1)  ||( second_value==1) )
      second_value+=10;
      
      /// Initialisation of the initial state
      state=(state_ini-12)*10+first_value-1+100*(isU);
      
      /// Getting the triplet pc dc and isUsable for the initial state
      getTriplet(pc,dc,isUsable);
     
      assert(pc==state_ini && dc==first_value && isUsable==isU);
      
      /// Setting the reward to zero
      reward=0.0;
}

Blackjack::~Blackjack()
{
      //dtor
      delete mdp;
}
