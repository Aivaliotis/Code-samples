#include "utils.h"
#include "Soil_params.h"
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <algorithm>
#include <vector>
#include <numeric>
#include <math.h>
#include "myclasses.h"
#include <omp.h>

using namespace std;

Individual Breed(Individual a, Individual b){
	Individual offspring;
	offspring = a;
	float k;
	for(int i=0; i<a.nVar; i++){
		k=(float)rand()/(RAND_MAX);
		offspring.Position[i] = k*a.Position[i]+(1-k)*b.Position[i];
	}
	return offspring;
}

// Function to return Index positions of sorted array
void sorted_Index(float a[],const int n, int* ptr){
	iota( ptr, ptr + n , 0 ) ;	// fill with 0, 1, ..., n-1
	sort( ptr, ptr + n, [&a] ( int i, int j ) { return a[i] > a[j] ; } ) ;	// indirect (descending) sort on indices
}

// Random permutation
void randperm( int maxNumber, int size, int *ptr){
    iota( ptr, ptr+maxNumber,0);
    random_shuffle(ptr, ptr + maxNumber);
}

float bishop(float *vars, Soil_params *Soil, Reinforcements *reinforcements, bool displayable, Circle *circle){
		float out;
	 	float x0, y0, R, xstart, xend;
	 	int Initial_Strip, Last_Strip;

	 	x0 = circle->x0;
	    y0 = circle->y0;
	    R = circle->R;
	    xstart = *(vars);
	    xend = *(vars+1);
		float Surface[Soil->nStrip + 4][2];
		float Natural_Surface[Soil->nStrip + 4][2];

	    if(reinforcements->fill_exist){
	    	Find_section_points(circle, Soil, xstart, xend, Surface, Natural_Surface);
	    }
	    else{
	    	for(int i = 0; i<Soil->nStrip; i++){
	    		Surface[i][0] = Soil->Surface[i][0];
	    		Surface[i][1] = Soil->Surface[i][1];
	    		Natural_Surface[i][0] = Soil->Natural_Surface[i][0];
	    		Natural_Surface[i][1] = Soil->Natural_Surface[i][1];
	    	}
	    	// creating dummy variables to restrain edges
	    	for(int i = 0; i < 4; i++){
				Surface[Soil->nStrip+i][0] = 3.40282e+038;
				Surface[Soil->nStrip+i][1] = 3.40282e+038;
				Natural_Surface[Soil->nStrip+i][0] = 3.40282e+038;
				Natural_Surface[Soil->nStrip+i][1] = 3.40282e+038;
			}
	    }
	    // Locating first and last strip
	    Initial_Strip = 0;
	    while( xstart >= Surface[Initial_Strip][0]){
	    	Initial_Strip++;
	    }

	    Last_Strip = Soil->nStrip +4 - 1;
		while(xend <= Surface[Last_Strip][0]){
			Last_Strip = Last_Strip - 1;
		}

		int nStrip = Last_Strip - Initial_Strip + 2;

		// Array initialization
		float Q[nStrip], x1[nStrip], x2[nStrip];
		float y1[nStrip], y2[nStrip], y3s[nStrip], y4s[nStrip], y3f[nStrip], y4f[nStrip];
		float heqs[nStrip], As[nStrip], Ws[nStrip],Ls[nStrip], heqf[nStrip], Af[nStrip], Wf[nStrip], Lf[nStrip], U[nStrip], a[nStrip];

		// Slice geometry
		for(int i = 0; i < nStrip; i++ ){
			x1[i] = Surface[Initial_Strip+i-1][0];	// x3 == x1
			x2[i] = Surface[Initial_Strip+i][0];	// x4 == x2
			y2[i] = y0-sqrt(R*R-(x2[i]-x0)*(x2[i]-x0));
			y4s[i] = Natural_Surface[Initial_Strip + i][1];
			if(i>0){
				y3s[i] = y4s[i-1];
				y1[i] = y2[i-1];}
			else{
				y3s[i] = Natural_Surface[Initial_Strip + i - 1][1];
				y1[i] = y0-sqrt(R*R-(x1[0]-x0)*(x1[0]-x0));
			}
		}

		// Correction of initial and last slice
		x1[0] = xstart;
		y1[0] = y0-sqrt(R*R-(x1[0]-x0)*(x1[0]-x0));
		x2[nStrip-1] = xend;
		y2[nStrip-1] = Natural_Surface[Initial_Strip+nStrip-1][1];
		y4s[nStrip-1] = Natural_Surface[Initial_Strip+nStrip-1][1];

		// Slice parameter calculation
		for(int i = 0; i < nStrip; i++ ){
			heqs[i] = ((y3s[i]+y4s[i])/2-(y1[i]+y2[i])/2)/2;
			As[i] = ((y3s[i]-y1[i]+y4s[i]-y2[i])*(x2[i]-x1[i])/2);
			Ws[i] = As[i]*Soil->wsoil;
			Q[i] = 0;
			U[i] = 0;	//Aw.*constants.wwater;
			a[i] = atan((y2[i]-y1[i])/(x2[i]-x1[i]));

			// Sliding surface propertie (natural soil / land fill)
			if(((y1[i] + y2[i])/2 < (y3s[i]+y4s[i])/2)){
				Ls[i] = sqrt((y2[i]-y1[i])*(y2[i]-y1[i])+(x2[i]-x1[i])*(x2[i]-x1[i]));
				Lf[i] = 0;
			}else{
				Lf[i] = sqrt((y2[i]-y1[i])*(y2[i]-y1[i])+(x2[i]-x1[i])*(x2[i]-x1[i]));
				Ls[i] = 0;
			}
		}

		// Fill contribution
		if(reinforcements->fill_exist){
			for(int i = 0; i < nStrip; i++ ){
				y3f[i] = Surface[Initial_Strip+i][1];
				y4f[i] = Surface[Initial_Strip+i+1][1];
				heqf[i] = ((y3f[i]+y4f[i])/2-(y1[i]+y2[i])/2)/2;
				Af[i] = ((y3f[i]-y1[i]+y4f[i]-y2[i])*(x2[i]-x1[i])/2);
				Wf[i] = Af[i]*Soil->w_fill;
			}
		}else{
			for(int i = 0; i < nStrip; i++ ){
				y3f[i] = 0;
				y4f[i] = 0;
				heqf[i] = 0;
				Af[i] = 0;
				Wf[i] = 0;
				Lf[i] = 0;
			}
		}

		// Anchor contribution
		float W_r[nStrip], T_r[nStrip], Lfree[nStrip], Lfixed[nStrip];

		if(reinforcements->anchor_exist){
			for(int i=0; i<nStrip; i++){
				Lfree[i] = reinforcements->Lfree[Initial_Strip+i-1];
				Lfixed[i] = reinforcements->Lfixed[Initial_Strip+i-1];

				if(Lfree[i] > 0 && Lfree[i]*Lfixed[i]>0){
					// section of the anchor line and the slip surface, providing the free length required for full anchor contribution
					float a_l = sin(-Soil->anchor_angle);
					float b_l = Soil->anchor_placement[Initial_Strip+i-1][1]-a_l*Soil->anchor_placement[Initial_Strip+i-1][0];
					float a2 = 1 + pow(a_l,2);
					float a1 = -2*circle->x0+2*a_l*(b_l-circle->y0);
					float a0 = pow(circle->x0,2) + pow((b_l - circle->y0),2)-pow(circle->R,2);
					float D = pow(a1,2)-4*a2*a0;
					float xe = (-a1 + sqrt(D))/(2*a2);
					float ye = a_l*xe + b_l;
					float Lreq = sqrt(pow((Soil->anchor_placement[Initial_Strip+i-1][0]-xe),2) + pow((Soil->anchor_placement[Initial_Strip+i-1][1] - ye),2));

					// anchor force reduction based on the relative position of the anchor bulb to the sliding surface
					float anchor_force=Lfixed[i]*3.14159265*(Soil->anchor_D)*Soil->anchor_qs*fmax(fmin((Lfree[i]+Lfixed[i]-Lreq)/Lfixed[i],1),0);
					W_r[i] = anchor_force*sin(Soil->anchor_angle);
					T_r[i] = anchor_force*cos(Soil->anchor_angle);
				}else{
					W_r[i] = 0;
					T_r[i] = 0;
				}
			}

			// revision of the contribution of the first anchor
			if( xstart>Soil->anchor_placement[Initial_Strip-1][0]){
				W_r[0] = 0;
				T_r[0] = 0;
			}
		}else{
			for(int i = 0; i < nStrip; i++){
				W_r[i] = 0;
				T_r[i] = 0;
			}
		}

		// Simplified Bishop method
		float FS[2], FS_fin, sum_nom1, sum_nom2, sum_den;
		int count;

		FS[0]=0;
		FS[1]=1;
		count=0;
		while( abs(FS[0]-FS[1])>0.00025&&count<25){
			FS[0]=FS[1];
			sum_nom1=0;
			sum_nom2=0;
			sum_den=0;
			for(int i = 0; i < nStrip; i++ ){
				sum_nom1 += ((Soil->c*Ls[i]*cos(a[i])+(Ls[i]>0)*((Ws[i]+Wf[i])*(1-Soil->kv)+W_r[i]+Q[i]-U[i]*cos(a[i]))*tan(Soil->fi_deg))/(cos(a[i])+tan(Soil->fi_deg)*sin(a[i])/FS[0]));
				sum_nom2 += ((Soil->c_fill*Lf[i]*cos(a[i])+(Lf[i]>0)*((Wf[i]+Ws[i])*(1-Soil->kv)+W_r[i]+Q[i]-U[i]*cos(a[i]))*tan(Soil->fi_fill))/(cos(a[i])+tan(Soil->fi_fill)*sin(a[i])/FS[0]));
				sum_den += (((Ws[i]+Wf[i])*(1-Soil->kv)+W_r[i]+Q[i])*sin(a[i])+Ws[i]*Soil->kh*(cos(a[i])-heqs[i]/R)+Wf[i]*Soil->kh*(cos(a[i])-heqf[i]/R)-T_r[i]*((y0-(y4s[i]+y3s[i])/2))/R);
			}
			FS[1]=((sum_nom1+sum_nom2)/sum_den);
			count=count+1;
		}

		FS_fin=FS[1];

		if(FS_fin<0||count==25||sum_den<0||(y0-R<Soil->H_bedrock)) {
			FS_fin = 1000; }	// non convergion solution }

		out = FS_fin;

		return out;
}

float determinant( float matrix[3][3], int n) {
   float det = 0;
   float submatrix[3][3];

   if (n == 2)
      return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
   else {
      for (int x = 0; x < n; x++) {
            int subi = 0;
            for (int i = 1; i < n; i++) {
               int subj = 0;
               for (int j = 0; j < n; j++) {
                  if (j == x)
                  continue;
                  submatrix[subi][subj] = matrix[i][j];
                  subj++;
               }
               subi++;

            }
            det = det + (pow(-1, x) * matrix[0][x] * determinant( submatrix, n - 1 ));
      }
   }
   return det;
}

void correct_circle(float vars[], Circle *circle, Soil_params *Soil){

	float dx, x1, x2, x3, theta, y1, y2, y3;

	bool early_exit = true;
	int i=0;
	x1 = vars[0];

	// adjust the logic parameters
	while(x1 > Soil->edges[i][0]) { i++;}
	int first_edge = i;
	float xnom[3][3], ynom[3][3], denom[3][3];

	dx = 1;
	vars[2] = vars[2]+1*3.14159265/180;	// adjustment for first iteration
	int count =0;

	do {
		count++;

		// adjust points for a logical slip surface (methodology was chosen through simple observations)
		if(early_exit) {
			vars[2] = vars[2]-1*3.14159265/180;
			if (vars[2]<(-60*3.14159265/180))
			{
				vars[2] = -60*3.14159265/180;
				vars[1] += 1;
			}
		}else{
			if(circle->y0<0) { vars[2] = vars[2] - 1*3.14159265/180; }
			else{
				vars[1] = min(vars[1] + 1, Soil->Natural_Surface[Soil->nStrip-1][0]);
				vars[2] = min(vars[2] + 1*3.14159265/180, Soil->slope-1*3.14159265/180);
			}
		}

		x1 = vars[0];
		x3 = vars[1];
		y1 = fmin(fmax(x1*tan(Soil->slope),0),Soil->H);
		y3 = fmin(fmax(x3*tan(Soil->slope),0),Soil->H);
		x2 = x1 + dx;
		theta = vars[2];
		y2 = y1 + tan(theta)*dx;

		// array creation for the circle determination from three points (although vector would have a more readible code, it quadrupled the total execution time)
		xnom[0][0] = x1*x1 + y1*y1;
		xnom[0][1] = y1;
		xnom[0][2] = 1;

		xnom[1][0] = x2*x2 + y2*y2;
		xnom[1][1] = y2;
		xnom[1][2] = 1;

		xnom[2][0] = x3*x3 + y3*y3;
		xnom[2][1] = y3;
		xnom[2][2] = 1;

		ynom[0][0] = x1*x1 + y1*y1;
		ynom[0][1] = x1;
		ynom[0][2] = 1;

		ynom[1][0] = x2*x2 + y2*y2;
		ynom[1][1] = x2;
		ynom[1][2] = 1;

		ynom[2][0] = x3*x3 + y3*y3;
		ynom[2][1] = x3;
		ynom[2][2] = 1;

		denom[0][0] = x1;
		denom[0][1] = y1;
		denom[0][2] = 1;

		denom[1][0] = x2;
		denom[1][1] = y2;
		denom[1][2] = 1;

		denom[2][0] = x3;
		denom[2][1] = y3;
		denom[2][2] = 1;

		circle->x0 = 0.5 * (determinant(xnom,3)/determinant(denom,3));
		circle->y0 = -0.5 * (determinant(ynom,3)/determinant(denom,3));
		circle->R = sqrt(pow((x1-circle->x0),2)+ pow((y1-circle->y0),2));

		early_exit = false;
		for(int i = first_edge; i<Soil->nedge; i++){
			if(circle->R < sqrt(pow((circle->x0 - Soil->edges[i][0]),2) + pow((circle->y0-Soil->edges[i][1]),2))) {
				early_exit = true; 
			}
		}

	} while(early_exit || (circle->y0 < y3));

	vars[0] = x1;
	vars[1] = x3;
	vars[2] = theta;
}

float Critical_FS_SGA(Soil_params *Soil, Reinforcements *reinforcements, Search_area *search_area, bool search_area_signal, float out_vars[], bool verbose){
	// initialization of heuristic enviroment
	float out = 0;	// output
	int nvar=3;	//variable number
	int Pop = 100;	// GA population
	int ElitNum = 0.2 * Pop;	// elite population
	int BreedNum = 0.6 * Pop;	// cross over population
	int nPop = 20;	// PSO population
	int MaxInteractions = 100;	// GA interactions
	int MaxIt = 5;	// PSO interactions

	float inMaxPosition[nvar] = {Soil->Soil[2][0] - 1, Soil->Natural_Surface[Soil->nStrip -1][0], Soil->slope-float(0.1*3.14159265/180)};
	float inMinPosition[nvar] = {Soil->Natural_Surface[0][0], Soil->Soil[2][0], -60 * 3.14159265/180};

	// Initialization of Area space
	if(search_area_signal){
		inMaxPosition[0] = search_area->x1max;
		inMinPosition[0] = search_area->x1min;
		inMaxPosition[1] = search_area->x3max;
		inMinPosition[1] = search_area->x3min;
		inMaxPosition[2] = search_area->thetamax;
		inMinPosition[2] = search_area->thetamin;
	}

	float MaxVelocity[nvar] , MinVelocity[nvar];

	for(int i =0; i<nvar; i++){
		MaxVelocity[i] = 0.2 * (inMaxPosition[i] - inMinPosition[i]);
		MinVelocity[i] = - MaxVelocity[i];
	}

	float GlobalBestValue = 3.40282e+038;
	float It_Value[MaxInteractions];
	float GlobalBestPosition[nvar];

	Individual particle[Pop];
	Individual new_particle[Pop];
	Circle circle[Pop]; // implementation of individual circles in order to be ready for parallelism as the other shared variabes are unaffected by the functions

	for(int i=0; i<Pop;i++){
		particle[i].initialize(nvar, inMaxPosition, inMinPosition, MaxVelocity, MinVelocity);
		particle[i].Value = evaluate_individual_FS(&particle[i], Soil, reinforcements, false, &circle[i]);
		particle[i].update_PB();
	}

	for(int i=0; i<Pop;i++){
		if(particle[i].Value<GlobalBestValue){
			for(int j=0; j<nvar; j++) { GlobalBestPosition[j] = particle[i].Position[j]; }
			GlobalBestValue = particle[i].Value;
		}
	}

	int interaction=0;
	int it = 0;
	int selection[Pop];
	int *selectionptr = selection;

	while (interaction <MaxInteractions){
		interaction++;
		it = 0;
		randperm(Pop, nPop, selectionptr);

		// Swarming process		
		while(it<MaxIt){
			// PSO exploration
			for(int i=0; i<nPop; i++){
				particle[selection[i]].update_velocity_position(GlobalBestPosition, it);
				particle[selection[i]].Value = evaluate_individual_FS(&particle[selection[i]], Soil, reinforcements, false, &circle[selection[i]]);
				particle[selection[i]].update_PB();
			}

			// Update Global Best
			for(int i=0; i<nPop; i++){
				if(particle[selection[i]].Value < GlobalBestValue){
					for(int j=0; j<nvar; j++) { GlobalBestPosition[j] = particle[selection[i]].Position[j]; }
					GlobalBestValue = particle[selection[i]].Value;
				}
			}

			it++;
		}

		for(int i=0; i<nPop; i++){
			particle[selection[i]].reset_velocity_PB();
		}

		//Genetic Process
		float fitness[Pop];
		for(int i=0; i<Pop; i++){
			fitness[i] = particle[i].Value;
		}

		int fitness_Index[Pop];
		int *fitness_Index_ptr = fitness_Index;
		sorted_Index(fitness, Pop, fitness_Index_ptr);

		for(int i=0; i<ElitNum; i++){
				new_particle[i] = Breed(particle[fitness_Index[i]],particle[fitness_Index[(rand()%ElitNum)]]);
			}

		for(int i=ElitNum; i<(BreedNum+ElitNum); i++){
			new_particle[i] = Breed(particle[fitness_Index[(rand()%Pop)]],particle[fitness_Index[(rand()%Pop)]]);
		}

		for(int i=(BreedNum+ElitNum); i<Pop; i++){
			new_particle[i] = particle[i];
			new_particle[i].mutate();
		}

		for(int i=0; i<Pop; i++){
			particle[i] = new_particle[i];
			particle[i].Value = evaluate_individual_FS(&particle[i], Soil, reinforcements, false, &circle[i]);

			// Update Global Best
			if(particle[i].Value<GlobalBestValue){
				for(int j=0; j<nvar; j++){ GlobalBestPosition[j] = particle[i].Position[j]; }
				GlobalBestValue = particle[i].Value;
			}
		}

		It_Value[interaction] = GlobalBestValue;

		if(verbose){ cout << "Interaction : " << interaction << "	Best Value : " << GlobalBestValue <<endl; }

		if(interaction>10 && It_Value[interaction]==It_Value[interaction-10]){
			break;
		}
	}
	
	for(int i=0; i<3; i++){out_vars[i] = GlobalBestPosition[i];}

	out = GlobalBestValue;

	return out;
}

float evaluate_individual_FS(Individual *particle, Soil_params *Soil, Reinforcements *reinforcements, bool displayable, Circle *circle){
	float vars[particle->nVar];

	for(int i=0; i< particle->nVar; i++){ vars[i] = particle->Position[i];	}

	correct_circle(vars, circle, Soil);

	for(int i=0; i< particle->nVar; i++){ particle->Position[i] = vars[i];	}
	return bishop(vars, Soil, reinforcements, displayable, circle);
}

// Find section points of the circle and the natural surface inside the slip surface, as slip sirface properties may differ 
void Find_section_points(Circle *circle, Soil_params *Soil, float xstart, float xend,float Surface[][2],float Natural_Surface[][2]){
	    float a = tan(Soil->slope);
	    float xtest[4], ytest[4];
	    int j=0;

	    // Search between the entry and the foot of the slope
	    if((circle->y0-circle->R)<0){
	    	if((circle->x0 + sqrt(pow(circle->R,2) - pow(circle->y0,2))) < 0){
				xtest[j] = circle->x0 + sqrt(pow(circle->R,2)-pow(circle->y0,2));
				ytest[j] = 0;
				j++;
	    	}
	    	if((circle->x0 - sqrt(pow(circle->R,2) - pow(circle->y0,2))) < xstart){
				xtest[j] = circle->x0 - sqrt(pow(circle->R,2) - pow(circle->y0,2));
				ytest[j] = 0;
				j++;
	    	}
	    }

	    // Search in the inclined section
	    float c1 = 1 + pow(a,2);
	    float c2 = -2 * circle->x0 - 2 * a * circle->y0;
	    float c3 = pow(circle->x0,2)  + pow(circle->y0,2) - pow(circle->R,2);

	    float D = pow(c2,2) - (4*c1*c3);
	    float sol[2] = {(-c2-sqrt(D)) / (2*c1), (-c2+sqrt(D)) / (2*c1) };

	    for(int i=0; i<2; i++){
	    	if((sol[i] > 0) && sol[i] < Soil->Soil[2][0]){
		    	xtest[j] = sol[i];
		    	ytest[j] = a*sol[i] ;
		    	j++;
		    }
	    }

	    // Introduce these points as new strips
	    for(int i=0; i < j; i++ ){
	    	Surface[i][0] = xtest[i];
	    	Surface[i][1] = ytest[i];
	    	Natural_Surface[i][0] = xtest[i];
	    	Natural_Surface[i][1] = fmax(xtest[i]*tan(Soil->slope),0);
	    }

	    for(int i=0; i < Soil->nStrip; i++ ){
			Surface[j+i][0] = Soil->Surface[i][0];
			Surface[j+i][1] = Soil->Surface[i][1];
			Natural_Surface[j+i][0] = Soil->Natural_Surface[i][0];
			Natural_Surface[j+i][1] = Soil->Natural_Surface[i][1];
		}

		// Fill empy sections points with illogical value
	    for(int i=j; i < 4; i++ ){
			Surface[Soil->nStrip+i][0] = 3.40282e+038;
			Surface[Soil->nStrip+i][1] = 3.40282e+038;
			Natural_Surface[Soil->nStrip+i][0] = 3.40282e+038;
			Natural_Surface[Soil->nStrip+i][1] = 3.40282e+038;
		}

		// bubble sort the strips 
	    for (int i = 0; i < j; i++){
	    	for (int ij = 0; ij < Soil->nStrip + 4 - i-1; ij++) {
	    		if (Surface[ij][0] > Surface[ij+1][0]){
	    			swap(Surface[ij], Surface[ij+1]);
	    			swap(Natural_Surface[ij], Natural_Surface[ij+1]);
	    		}
	    	}
	    }
}