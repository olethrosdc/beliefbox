// copyright (c) 2015 by Hannes Eriksson <hannese18@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
 #include "Vector.h"
 #include "GaussianProcessBandit.h"
 #include "LinearBanditTS.h"
 
 // test of GP-UCB and LinearBandit-TS
 
 int main(int argc, char* argv[]) {
	 Matrix X = Matrix(25,3);
	 Vector Y = Vector(25);
	 X(0,0) = 0; X(0,1) = 0; X(0,2) = 0; Y(0) = 8;
	 X(1,0) = 1; X(1,1) = 0; X(1,2) = 0; Y(1) = 7;	 
	 X(2,0) = 0; X(2,1) = 1; X(2,2) = 0; Y(2) = 7;	
	 X(3,0) = 0; X(3,1) = 0; X(3,2) = 1; Y(3) = 7;	 
	 X(4,0) = 1; X(4,1) = 1; X(4,2) = 0; Y(4) = 6.5;	 
	 X(5,0) = 1; X(5,1) = 0; X(5,2) = 1; Y(5) = 6.5;
	 X(6,0) = 0; X(6,1) = 1; X(6,2) = 1; Y(6) = 6.5;
	 X(7,0) = 1; X(7,1) = 1; X(7,2) = 1; Y(7) = 6;	 
	 X(8,0) = 0.25; X(8,1) = 0.25; X(8,2) = 0; Y(8) = 7.5;	
	 X(9,0) = 0; X(9,1) = 0.25; X(9,2) = 0.25; Y(9) = 7.5;	 
	 X(10,0) = 0.25; X(10,1) = 0; X(10,2) = 0.25; Y(10) = 7.5;	 
	 X(11,0) = 0.5; X(11,1) = 0.5; X(11,2) = 0; Y(11) = 7.1;
	 X(12,0) = 0.5; X(12,1) = 0; X(12,2) = 0.5; Y(12) = 7.1;
	 X(13,0) = 0; X(13,1) = 0.5; X(13,2) = 0.5; Y(13) = 7.1;	 
	 X(14,0) = 0.75; X(14,1) = 0; X(14,2) = 0; Y(14) = 7.25;	
	 X(15,0) = 0; X(15,1) = 0.75; X(15,2) = 0; Y(15) = 7.25;	 
	 X(16,0) = 0; X(16,1) = 0; X(16,2) = 0.75; Y(16) = 7.25;	 
	 X(17,0) = 2; X(17,1) = 0; X(17,2) = 0; Y(17) = 5;
	 X(18,0) = 0; X(18,1) = 2; X(18,2) = 0; Y(18) = 5;
	 X(19,0) = 0; X(19,1) = 0; X(19,2) = 2; Y(19) = 5;
	 X(20,0) = 2; X(20,1) = 2; X(20,2) = 0; Y(20) = 3;	 
	 X(21,0) = 2; X(21,1) = 0; X(21,2) = 2; Y(21) = 3;	
	 X(22,0) = 0; X(22,1) = 2; X(22,2) = 2; Y(22) = 3;	 
	 X(23,0) = 2; X(23,1) = 2; X(23,2) = 2; Y(23) = 2.5;	 
	 X(24,0) = 0; X(24,1) = 3; X(24,2) = 0; Y(24) = 4;

	 Vector gp_scale_length = Vector(3) + 1;
 	 Vector gp_scale_length_dic = Vector(3) + 1;
	 real gp_noise_var = 1.0;
	 real gp_sig_var = 1.0;
	 real gp_threshold = 1.0;
	 real delta = 0.5;

	 GaussianProcessBandit gpUcb = GaussianProcessBandit(gp_scale_length,
		gp_scale_length_dic,gp_noise_var, gp_sig_var, 
		gp_threshold, delta);
	 for(int i=1;i<9;i++) {
		 gpUcb.addObservation(X.getRow(i),Y(i));
	 }
	 Vector predictions = Vector();
	 real regret = 0;
	 for(int i=0;i<1000;i++) {
		int bestBandit = gpUcb.predict(X, predictions);
		gpUcb.addObservation(X.getRow(bestBandit),Y(bestBandit));
		regret+=(8 - Y(bestBandit));
	 }
	 std::cout<<"\n -----\n"<<"GPB regret: "<<regret/1000<<"\n -----";
	 std::cout<<"\n";

	 real R = 5;
	 real epsilon = 0.5;

	 LinearBanditTS lbTS = LinearBanditTS(R, delta, epsilon);
	 for(int i=1;i<9;i++) {
		 lbTS.addObservation(X.getRow(i),Y(i));
	 }
	 predictions = Vector();
	 regret = 0;
	 for(int i=0;i<1000;i++) {
		int bestBandit = lbTS.predict(X, predictions);
		lbTS.addObservation(X.getRow(bestBandit),Y(bestBandit));
		regret+=(8 - Y(bestBandit));
	 }
	 std::cout<<"\n -----\n"<<"LBTS regret: "<<regret/1000<<"\n -----";
	 std::cout<<"\n";

 }