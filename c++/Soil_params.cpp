
#include "Soil_params.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

void Soil_params::update_params(float H_in, float slope_in, float wsoil_in, float c_in, float fi_deg_in){
	H = H_in;
	slope= slope_in * 3.14159265/180;
	Soil[2][0] = H/((float)tan(slope));
	Soil[2][1] = H;
	Soil[3][0] = fmax(10000, 5*Soil[2][0]);
	Soil[3][1] = H;
	wsoil = wsoil_in;
	c = c_in;
	fi_deg = fi_deg_in * 3.14159265/180;
}

void Soil_params::Strip_creation(){

    int nslices = fmax( ceil(Soil[2][0]-Soil[1][0]) , 20);
    float dxStrip = (Soil[2][0]-Soil[1][0])/nslices;
    float dxStrip_plateau = 1;
    int i;

	float xmin = -ceil(H * 10);
	float xmax = ceil(fmax( 2 * H / tan(0.5 * slope) + 0.1 , 2 * H));

	nStrip_foot = ceil(abs(xmin / dxStrip_plateau));
	nStrip_crown = ceil((xmax-Soil[2][0]) / dxStrip_plateau);
	nStrip = nStrip_foot + nStrip_crown + nslices;

	Natural_Surface = new float*[nStrip];
	for (int i = 0; i < nStrip; i++) {
		Natural_Surface[i] = new float[2];
	}

	Surface = new float*[nStrip];
	for (int i = 0; i < nStrip; i++) {
		Surface[i] = new float[2];
	}

	anchor_placement = new float*[nStrip-1];
	for (int i = 0; i < nStrip-1; i++) {
		anchor_placement[i] = new float[2];
	}

	Natural_Surface[0][0] = xmin;
    Natural_Surface[0][1] = 0;
    i = 1;

    while (i < nStrip){
        if( Natural_Surface[i-1][0] >= 0 && Natural_Surface[i-1][0]<Soil[2][0]){
            Natural_Surface[i][0] = Natural_Surface[i-1][0] + dxStrip;
            Natural_Surface[i][1] = Natural_Surface[i][0]*tan(slope);

            if (Natural_Surface[i][0] == Soil[2][0]){ Natural_Surface[i][1] = H; }
        }
        else{
            Natural_Surface[i][0] = Natural_Surface[i - 1][0] + dxStrip_plateau;

            if (Natural_Surface[i-1][0] < 0){Natural_Surface[i][1] = 0;}
            else {Natural_Surface[i][1] = H;}
        }

        anchor_placement[i-1][0] = (Natural_Surface[i-1][0] + Natural_Surface[i][0])/2;
		anchor_placement[i-1][1] = (Natural_Surface[i-1][1] + Natural_Surface[i][1])/2;
        i++;
        }
}

void Soil_params::update_surface(Reinforcements *a){
	int i;

	// Update of the soil surface height at each strip
	for(i=0; i<nStrip; i++){
		Surface[i][0] = Natural_Surface[i][0];
		Surface[i][1] = Natural_Surface[i][1] + a->fill[i]*H_fill;
	}

	// reset previous allocated edges vector
	edges.clear();
	float b[nStrip];

	// Update the soil edges array
	for(i=0; i < nStrip-2; i++){
	        b[i]=(Surface[i+1][1] - Surface[i][1]) / (Surface[i+1][0] - Surface[i][0]);

	        if(i>0 && (abs(b[i]-b[i-1])>0.001)){
	            edges.push_back({Surface[i][0], Surface[i][1]});
			}
	}

	nedge = edges.size();
}

// Clear dynamicly allocated variables
Soil_params::~Soil_params(){
	delete Natural_Surface;
	delete Surface;
}
