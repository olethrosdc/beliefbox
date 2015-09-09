#pragma once
#include "../Board.h"
#include "Distribution.h"

class Node{
	protected:
		Board board;
		Distribution gDis;
		double gValue;
		double vValue;
		Distribution delta;
                Distribution vDis;
		bool visited;
		char color;
		std::vector<Node*> children;
		Node * parent;	
		Distribution messageToParent;
		Distribution messageFromParent;
		Distribution rollOut;
	public:
		Node(Board &b,Distribution &d,const char color ,const double gValue=0.0,const double vValue=0.0);
		~Node();
		double getG() const;
		double getV() const;
		char getColor();
		Distribution &getGDis() ;
		Distribution &getVDis() ;
		Distribution &getDelta();

		Distribution &getMessageToParent() ;
		Distribution &getMessageFromParent() ;
		Distribution &getRollOut() ;
	        Board &getBoard();

		Node * getParent() ;
		Node * getChild() ;	
		std::vector<Node*> getChildren();

		void setParent(Node * n);
		void setMessageToParent(const Distribution &dis);
		void setMessageFromParent(const Distribution &dis);
		void setVisited();
		void setDelta(const Distribution &dis);
		void setGDis(const Distribution &dis);
		void setVDis(const Distribution &dis);
		void setRollOut(const Distribution &dis);
		bool isVisited() const;
		void sampleG();
	;	

};
