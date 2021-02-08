#include "myclasses.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <limits>

using namespace std;
// Reinforcement Functions
void Reinforcements::initialize(int a){
	//  Dynamic pointers were used over vectors as they have superior performance, 
	// 	need to be deconstructed only once per slope profile and for personal practice with pointers
	fill = new float[a];
	Lfree = new float[a - 1];
	Lfixed = new float[a - 1];
	for(int i=0; i<a - 1; i++){
		fill[i] = 0;
		Lfree[i] = 0;
		Lfixed[i] = 0;
	}
	fill[a-1] =0;
}

// Clear dynamicly allocated variables
Reinforcements::~Reinforcements(){
	delete fill;
	delete Lfree;
	delete Lfixed;
}

// SGA Functions
void Individual::initialize(int k, float inMaxPosition[], float inMinPosition[], float inMaxVelocity[], float inMinVelocity[]){
	nVar = k;
	Value = 0;
	PBV = numeric_limits<float>::max();

	for( int i=0; i<k;i++){
		MaxPosition.push_back(inMaxPosition[i]);
		MinPosition.push_back(inMinPosition[i]);
		MaxVelocity.push_back(inMaxVelocity[i]);
		MinVelocity.push_back(inMinVelocity[i]);
		Position.push_back(MinPosition[i]+(float)rand()/(RAND_MAX) * (MaxPosition[i]-MinPosition[i]));
		Velocity.push_back(0);
	}
	PBP = Position;
}

void Individual::reset_velocity_PB(){
	for( int i=0; i<nVar;i++){
		Velocity[i] = 0;
	}
	PBV = Value;
	PBP = Position;
}

void Individual::mutate(){
	for( int i=0; i<nVar;i++){
			Position[i] = MinPosition[i]+(float)rand()/(RAND_MAX) * (MaxPosition[i]-MinPosition[i]);
		}
}

void Individual::update_velocity_position(float GBP[] , int it){
	float personal_velocity[nVar];
	float social_velocity[nVar];

	//Update Velocity
	for( int i=0; i<nVar;i++){
			personal_velocity[i] = c1 * (float)rand()/(RAND_MAX) * ( PBP[i] -Position[i] );
			social_velocity[i] = c2 * (float)rand()/(RAND_MAX) * (GBP[i] - Position[i]);
			Velocity[i] = w * pow(wdamp ,it) * Velocity[i] + personal_velocity[i] + social_velocity[i];

			// Apply Velocity boundary conditions
			Velocity[i] = max( min(Velocity[i], MaxVelocity[i]), MinVelocity[i]);
		}

	//Update Position
	for( int i=0; i<nVar;i++){
		Position[i] = Position[i] + Velocity[i];

		//Apply Position boundary conditions
		Position[i] = max( min(Position[i], MaxPosition[i]), MinPosition[i]);
	}
}

void Individual::update_PB(){
	if(Value < PBV){
		PBV = Value;
		PBP = Position;
	}
}

// Clear dynamicly allocated variables, avoided by the use of vectors as it was computational intensive
// 1.8 times the total execution time to use: delete <pointer>
Individual::~Individual(){	
}

