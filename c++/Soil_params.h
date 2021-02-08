#ifndef POBLEMPARAMS_H_
#define POBLEMPARAMS_H_

#include <math.h>
#include <vector>
#include "myclasses.h"

class Soil_params {
public:
	// Functions
	~Soil_params();
	void update_params(float H_in, float slope_in, float wsoil_in, float c_in, float fi_deg_in);
	void Strip_creation();
	void update_surface(Reinforcements *a);

	//Soil
	float H=20;
	float slope= 45 * 3.14159265/180;
	float Soil[4][2] = {{-1000,0}, {0,0}, {H/((float)tan(slope)), H}, {10000, H}};
	float wsoil = 20;
	float c = 30;
	float fi_deg = 20 * 3.14159265/180;
	int n_soil = 4;
	int nStrip;
	float *Strip;
	int nstrips;
	int nStrip_foot;
	int nStrip_crown;
	float H_bedrock=-10000000;
	float **Natural_Surface;
	float **Surface;
	std::vector<std::vector<float>> edges;
	int nedge;

	//Fill
	float w_fill = 20;	// KN/m^3
	float H_fill = 10;	// m
	float slope_fill = 45 * 3.14159265/180;	// deg
	float c_fill = 5;	// Kpa
	float fi_fill = 40 * 3.14159265/180;	// deg
	float Fill_cost=0;	// euro/m^3
	float Fill_length=1;

	// Water Horizon
	float Water_horizon=-10000000;
	float wwater = 9.81;
	float nw = sizeof(Water_horizon)/sizeof(float);

	// Additional loads
	float q_start;     //x-coord (m)
	float q_end;       //y-coord (m)
	float q_value;     //Kpa

	// Kinematic equations
	float kv=0;
	float kh=0;

	// Cost Variables
	float cost_C1=200;
	float cost_C2=50; 	//euro/m
	float cost_C3=0;
	float cost_C4=0.025; 	//euro/kNm
	float cost_modifier=0.01; //ratio for obj function

	// Anchor Variables
	float anchor_D=0.2; 	//m
	float anchor_qs=200; 	//kPa
	float Anchor_factor=100;
	float L_precission=0.5;
	float Lfree_max=35;
	float Lfixed_max=5;
	float anchor_spacing=2; //m from center to center
	float anchor_first=0; //
	float anchor_angle=10* 3.14159265/180;
	float **anchor_placement;

	// Initial Population Multiplier
	float p=3;
};

#endif /* POBLEMPARAMS_H_ */


