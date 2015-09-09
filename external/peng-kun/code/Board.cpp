#include "Board.h"
#include "ComputerPlayer.h"
#include <iostream>
#include "Constants.h"
#include <vector>
using namespace std;

/*size is the length of a side, like 3,4 or 5, but it should be less than 
 10.And board is the char array to keep pieces on the board*/ 
Board::Board(int size){
	this->size = size;
	this->initialize();
}

Board::Board(){
	;
}

/*this copy function is mandatory if algorithm strategy is used*/
Board::Board(Board &b){
	size = b.size;
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			board[i][j] = b.board[i][j];
}

/*show the board*/
void Board::display(){
	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++)
			cout << board[i][j] << ' ';
		cout << endl;
	}

	cout << endl << endl;
}

/*initialize the board*/
void Board::initialize(){
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
			board[i][j] = '_';
	
}

/*return all possible next board states.Argument color specify who make moves*/
std::vector<Board*> Board::nextBoardStates(const char color){
	std::vector<Board*> states(0);
	if(checkWin(color==BLACK?WHITE:BLACK))
		return states;
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++){
			if(board[i][j]==EMPTY){
				Board * s = new Board(size);
				*s = *this;
				s->makeMove(i,j,color);
				states.push_back(s);
			}
		}

	return states;
}

int Board::numberOfNextBoardStates(const char color){
	int counter=0;
	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++){
			if(board[i][j] == EMPTY){
				counter++;	
			}
		}
	}
	return counter;
}
/*check whether the player with one specific color wins*/
bool Board::checkWin(const char color){
	char c=color;
	bool win;
	int i,j;
	/*vertical*/
	for(j=0;j<size;j++){
		win = TRUE;
		for(i=0;i<size;i++){
			if(board[i][j] != c){	
				win = FALSE;
				break;
			}		
		}
		if(win)
			return win;
	}
	/*horizontal*/
	for(i=0;i<size;i++){
		win = TRUE;
		for(j=0;j<size;j++){
			if(board[i][j] != c){
				win =FALSE;
				break;
			}
			
		}
		if(win)
			return win;
	}


	/*diagonal*/
	win = TRUE;
	for(i=0,j=0;i<size;i++,j++){
		if(board[i][j] != c){
			win = FALSE;
			break;
		}
	}
	if(win)
		return win;

	win = TRUE;
	for(i = 0,j=size-1;i<size;i++,j--){
		if(board[i][j]!=c){
			win = FALSE;
			break;
		}
		
	}
	return win;
}

/*add one piece onto the board*/
void Board::makeMove(const int i,const int j, const char color){
	board[i][j] = color;

}



int main(){
	int rounds=100;
	int winTimes =0;
	while(rounds-->0){
		Board start(4);
		start.display();

		ComputerPlayer player1(WHITE);
		ComputerPlayer player2(BLACK);

		char result;
	
		while(1){
			result = player1.algorithmPlay(start);
			start.display();
			if(result != CONTINUE){
				std::cout << result << " win" <<std::endl;
				break;
			}
			result = player2.randomPlay(start);
			start.display();
			if(result != CONTINUE){
				std::cout << result << " win" << std::endl;
				break;
			}

		}
		if(result == WHITE)
			winTimes ++;
	}
	std::cout << winTimes << std::endl;
}
