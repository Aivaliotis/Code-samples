#ifndef MYCLASSES_H_
#define MYCLASSES_H_

#include <vector> 

class Bishop_c {
private:
public:
	int Initial_Strip;
	int Last_Strip;
};

class Reinforcements {
public:
	~Reinforcements();
	void initialize(int a);
	float *fill;
	float anchors;
	float *Lfree;
	float *Lfixed;
	bool fill_exist = false;
	bool anchor_exist = false;
};

class Circle {
public:
	float x0=0;
	float y0=0;
	float R=0;
};

class Individual {
private:
	// float *PBP;				//Personal Best Position
	std::vector<float> PBP;
	float PBV;				//Personal Best Value
	float c1=1;
	float c2=1;
	float w=1;
	float wdamp=0.9;
	float wmin = 0.01;
	std::vector<float> MaxPosition;
	std::vector<float> MaxVelocity;
	std::vector<float> MinPosition;
	std::vector<float> MinVelocity;

public:
	~Individual();
	void mutate();
	void initialize(int k, float inMaxPosition[], float inMinPosition[], float inMaxVelocity[], float inMinVelocity[]);
	void update_velocity_position(float GBP[], int it);
	void update_PB();
	void reset_velocity_PB();
	float Value;
	std::vector<float> Position;
	std::vector<float> Velocity;
	int nVar;
};

class Search_area {
public:
	float x1max;
	float x1min;
	float x3max;
	float x3min;
	float thetamax;
	float thetamin;
};
#endif /* POBLEMPARAMS_H_ */






