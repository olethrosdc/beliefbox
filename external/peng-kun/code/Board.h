#pragma once
#include <vector>
class Board{
	public:
		int size;
		char board[10][10];

		Board(const int size);
		Board(Board &b);
		Board();
	
		void initialize();
		void display();
		std::vector<Board*> nextBoardStates(const char color);
		void makeMove(const int i,const int j,const char color);
		bool checkWin(const char color);
		int numberOfNextBoardStates(const char color);
};
