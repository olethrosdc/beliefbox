#include "Ep.h"
#include "../Constants.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "../ComputerPlayer.h"
	
/*top level class of the core algorithm, it provides interface to Player
  class, they are well modeled, the granity of this class is object Node
*/


Distribution Ep::descent(Node &n){
	if(n.isVisited()){
		Node* c = n.getChild();
		if(c!=NULL){
			/*update message to c*/
			Ep::updateMessageFromParent(*c,n);
			/*continue to descent*/
			Distribution newMessage = descent(*c);
			/*update marginal of n*/
			Ep::updateMessageToParent(*c,n,newMessage);
		}
	}

	else{
		/*update message, include dividing roll out out of the distribution*/
		Node * parent  = n.getParent();
		Ep::updateMessageExceptRollOut(n,*parent);
		
		/*dealing with roll out*/
		/*assume this is returned by the roll out*/
	        
		int length = 0;
		int result = 0;
		char lastPlayer;
		std::vector<int> stored(0);
		Ep::doRollout(n,length,result,lastPlayer,stored);
		/*call function to calculated the approximate message*/
		Distribution rollOut = Ep::getMessageFromRollOut(n,length,result);	
		n.setRollOut(rollOut);
		/*check whether the rollout-ed node is a leaf node*/
		if(length==0)
			n.setGDis(rollOut);
		else
			n.setGDis(n.getGDis()*rollOut);

		/*update V distribution =  delta + V*/
		Distribution v = Ep::calculateDelta(stored,lastPlayer);
		n.setDelta(v);
		n.setVDis(v+n.getGDis());
		n.setVisited();
	}
	Distribution messageExceptP = n.getGDis()/n.getMessageFromParent();
	Ep::updateParentVDistribution(n);
	return Distribution(messageExceptP.getMean(),messageExceptP.getVar()+1);

}

/*change as different board and computer player implementation*/
void Ep::doRollout(Node &n,int &l,int &r,char &lastPlayer,std::vector<int> &stored){
	Board b = n.getBoard();
	char color = n.getColor();
	ComputerPlayer p(color);
	int num = 0;
	char result;
	if(b.checkWin(p.getColor()==WHITE ? BLACK : WHITE)){
		result = p.getColor()==WHITE ? BLACK : WHITE;
	}
	else{
	while(1){
		 num = b.numberOfNextBoardStates(p.getColor());
		if(num==0){
			p.setColor(p.getColor()==WHITE ? BLACK : WHITE);
			result=BLACK;
			break;
		}
		stored.push_back(num);
		result = p.randomPlay(b);
		if(result != CONTINUE)
			break;
		p.setColor(p.getColor()==WHITE ? BLACK : WHITE);
	}
	}
	l = p.step; 
	lastPlayer =  p.getColor();
	r = result == 'X' ? 1 : -1;
}

/*Assume branch factor and depth are constant,which obviously are not*/
Distribution Ep::calculateDelta(std::vector<int> &branch,char lastPlayer){
	if(branch.size()==0)
		return Distribution(0,0);

	Distribution delta(0,0);
	for(std::vector<int>::iterator it = branch.end()-1;it>=branch.begin();it--){
		int branchFactor = (*it);
		std::vector<Distribution> variables(0);
		for(int i=0;i<branchFactor;i++){
			variables.push_back(delta+Distribution(0,1));

		}
		if(lastPlayer == MAX){
			delta =  Distribution::getMaxOfIndependentSet(variables);
			lastPlayer = MIN;
		}
		else{
			delta =  Distribution::getMinOfIndependentSet(variables);
			lastPlayer = MAX;
		}
	}
	return delta;
}

/*Update the V distribution of parent*/
void Ep::updateParentVDistribution(Node &n){
	Node *parent = n.getParent();
	if(parent==NULL)
		return;
	
	std::vector<Node*> children = parent->getChildren();
        //when delta are unknown,temporarily use average delta to approximate the real delta
	for(std::vector<Node*>::iterator it=children.begin();it<children.end();it++){
		if((*it)->getDelta()==Distribution(0,0)){
			int length = 0;
			int result = 0;
			char lastPlayer;
			std::vector<int> stored(0);
			Ep::doRollout((**it),length,result,lastPlayer,stored);
			Distribution v = Ep::calculateDelta(stored,lastPlayer);
			(*it)->setDelta(v);
		}
	}

	Distribution newV;
	if(parent->getColor()==MAX)
		newV = Distribution::getMaxOfCorrelatedSet(children);
	else
		newV = Distribution::getMinOfCorrelatedSet(children);
	parent->setVDis(newV);
}

/*Get message from RollOut result*/
Distribution Ep::getMessageFromRollOut(Node &boundary,const int length,const int result){
	Distribution prior = Distribution(boundary.getGDis().getMean(),boundary.getGDis().getVar()+length);
	double priorMean = prior.getMean();
	double priorVar = prior.getVar();
	double priorDev = sqrt(priorVar);

	double firstMoment,secondMoment,p,k,n,d,v;
	double meanSqaure = pow(priorMean,2);
	k = priorMean/priorDev;
	if(result>0){
		/*first moment*/
		n = Distribution::normalPdf(k);
		d = Distribution::normalPhi(k);
		/*seconde moment*/
		p = (1+k*n/d);
		firstMoment = priorMean + priorDev*(n/d);
		
	}
	else{
		/*first moment*/
		n =  Distribution::normalPdf(k);
		d =  Distribution::normalPhi(-k);
		/*second moment*/
		p  = (1-k*n/d);
		firstMoment = priorMean - priorDev*(n/d);
		}
	secondMoment =  meanSqaure+ priorVar*p;
	/*variance*/
	v = secondMoment - firstMoment*firstMoment;
	return Distribution(firstMoment,v+length);
}

/*Update the message from its parent*/
void Ep::updateMessageFromParent(Node &child,Node &parent){
	child.setGDis(child.getGDis()/child.getMessageFromParent());
	Distribution messageExceptC = parent.getGDis()/child.getMessageToParent();
	Distribution messageIntegral = Distribution(messageExceptC.getMean(),messageExceptC.getVar()+1);
	child.setMessageFromParent(messageIntegral);
	child.setGDis(child.getGDis()*messageIntegral);
}

/*Update message from rollout*/
void Ep::updateMessageExceptRollOut(Node &child,Node &parent){
	Distribution newParent = parent.getGDis()/parent.getRollOut();
	child.setGDis(child.getGDis()/child.getMessageFromParent());
	Distribution messageExceptC = newParent/child.getMessageToParent();
	Distribution messageIntegral = Distribution(messageExceptC.getMean(),messageExceptC.getVar()+1);
	child.setGDis(child.getGDis()*messageIntegral);
	child.setMessageFromParent(messageIntegral);
}

/*update the messsage to parent*/
void Ep::updateMessageToParent(Node &child,Node &parent,Distribution &message){
	
	parent.setGDis(parent.getGDis()/child.getMessageToParent()*message);
	child.setMessageToParent(message);

}

