#pragma once
#include "Distribution.h"
#include "Node.h"

class Ep{
	public:
		static void updateMessageFromParent(Node &child,Node &parent);
		static void updateMessageToParent(Node &child,Node &parent,Distribution &message);
		static void updateMessageExceptRollOut(Node &child,Node &parent);
		static Distribution descent(Node &n);
		
		static Distribution getMessageFromRollOut(Node &b,const int len,const int result);
		static void doRollout(Node &b,int &length,int &result,char &lastPlayer,std::vector<int> &stored);
		static Distribution calculateDelta(std::vector<int> &branch,char lastPlayer);
		static void updateParentVDistribution(Node &n);

};
