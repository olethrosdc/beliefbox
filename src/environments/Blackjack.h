/**
* Blackjack Environment for the popular casino game Blackjack
* described in http://webdocs.cs.ualberta.ca/~sutton/book/ebook/node51.html
*
* STATES:
* The state of the environment consists of 3 variables
* pc: Current player sum ( between 12 and 21 )
* dc: Current dealer showing card ( between 1 and 10, 1 for ace )
* isUsable: 1 if a value of 11 has been used in the play for the ace, 0 otherwise
* an additional state is added for the terminal state
*
*
*ACTIONS:
* There is two actions 0 or 1
* 0 if the player chooses to hit
* 1 if the player chooses to stick
*/


#ifndef BLACKJACK_H
#define BLACKJACK_H

#include "DiscreteMDP.h"
//#include "MDP.h"
#include "Environment.h"
#include "MersenneTwister.h"
#include <string>
#include <vector>


class Blackjack : public Environment<int,int>
{
      public:


            DiscreteMDP* mdp;


            Blackjack(real win_value_=1.0, real draw_value_=0.0, real loss_value_=-1.0);

            virtual DiscreteMDP* getMDP() const;

            virtual void Reset();
            virtual bool Act(const int& action);
            //void Show();
            
            int getState1() const
            {
            		return state;
            }
            

            /**
            * Computing the reward of the player given
            * the state triplet and the value of isNatural
            * -1 for a loss, 1 for a win, 0.0 for a draw
            *
            */
            real calculReward() const;

            /**
            * Given the state triplet pc, dc and isUsable
            * state returns a unique identification number
           */
            int getState(int pc_,int dc_,int isUsable_) const;

            /**
            *  Compute the corresponding triplet pc,dc,isUsable
            *  from the current state
            */

            void getTriplet(int& pc_,int& dc_,int& isUsable_);

            /**
            * Return the name of the environment
            */

            virtual const char* Name()
            {
                  return "Blackjack";
            }
            virtual ~Blackjack();
      protected:
            real win_value;
            real draw_value;
            real loss_value;
            int terminal_state;
            int pc,dc,isUsable,second_value,first_value;
            /// True is the first two cards of the player is a Face and an Ace
            bool isNatural;
            RandomNumberGenerator *rng;

      private:
};

#endif // BLACKJACK_H
