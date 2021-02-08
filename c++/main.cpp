#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <ctime>
#include "utils.h"
#include "Soil_params.h"
#include "Soil_params.cpp"
#include "myclasses.h"
#include "myclasses.cpp"
#include "utils.h"
#include "utils.cpp"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>
#include <fstream>

using namespace std;

int	main() {
	
	clock_t tStart = clock();
	int counter1=0;
	int maxiteration = 100;
	// Open or Create new txt file to save simplified bishop results
 	ofstream file;
 
    file.open("Unreinforced Soil dataset.txt",ios_base::app);
    if(!file)
    {
      cout<<"Error in creating file.."<<endl;
      return 0;
    }
    cout<<"\nFile created successfully."<<endl;
	
	// multithreaded process spread 8 threads to increase efficiency, setting every variable inside the pragma (private) to avoid memory overhead except a shared counter to follow the progress
    omp_set_num_threads(8);
	#pragma omp parallel shared(counter1)
	{
		Soil_params Soil;
		srand(int(time(NULL)) ^ omp_get_thread_num());	// setting a distinct random seed for each thread
		float H,slope,wsoil,c,fi_deg;

	#pragma omp for //schedule(dynamic)
		for(int counter=0; counter<maxiteration; counter++){

			// setting random slope properties to create a wide spectrum of solutions
			H=(float)rand()/(RAND_MAX)*30+10;	//m
			slope=(float)rand()/(RAND_MAX)*40+10;	//deg 
			wsoil=(float)rand()/(RAND_MAX)*20+10;	//KN/m3
			c=(float)rand()/(RAND_MAX)*30;	//KPa
			fi_deg=(float)rand()/(RAND_MAX)*30+5;	//deg


			Soil.update_params(H, slope, wsoil, c, fi_deg);
			Reinforcements reinforcements;
			Circle test_circle;
	
			Soil.Strip_creation();
	
			reinforcements.initialize(Soil.nStrip);
			Soil.update_surface(&reinforcements);
	
			float test_vars[3];
			// Slip surface search area based on entry point, exit point and entry angle
			Search_area search_area;
			search_area.x1max = Soil.Soil[2][0] - 1;
			search_area.x1min = Soil.Natural_Surface[0][0];
			search_area.x3max = Soil.Natural_Surface[Soil.nStrip -1][0];
			search_area.x3min = Soil.Soil[2][0];
			search_area.thetamax = Soil.slope-(0.1*3.14159265/180);
			search_area.thetamin = -60 * 3.14159265/180;
				
			float FS;
	
			FS = Critical_FS_SGA(&Soil, &reinforcements, &search_area, true, test_vars,false);

			// creating critical section for threads as only one thread must execute this command at the time
			#pragma omp critical
			{
				file << H <<" "<<slope<<" "<<wsoil<<" "<<c<<" "<<fi_deg<< " " << test_vars[0] << " " << test_vars[1] << " " << test_vars[2] <<" "<<FS << endl; // output in txt file
				counter1++;
			}

			cout << counter1 << "/" << maxiteration << endl;	// display counter to verify that it runs properly
		}
	}
	cout << float(clock() - tStart)/1000 << endl;
}