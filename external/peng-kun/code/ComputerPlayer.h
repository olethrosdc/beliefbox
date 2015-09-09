#pragma once
#include <random>
class Board;
class ComputerPlayer{
	protected:
		char color;
		std::mt19937 seed;
	public:
		int step;
		char randomPlay(Board &b);
		ComputerPlayer(const char color);
		char algorithmPlay(Board &b);
		void setColor(const char color);
		char getColor() const;
};
