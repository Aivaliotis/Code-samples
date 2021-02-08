#ifndef UTILS_H_
#define UTILS_H_
#include "myclasses.h"
#include "Soil_params.h"

Individual Breed(Individual a, Individual b);

void randperm(int a, int b, int* ptr);

void sorted_Index(float a[], int n, int* ptr);

float determinant( float matrix[10][10], int n);

void correct_circle(float vars[], Circle *circle, Soil_params *Soil);

float bishop(float *vars, Soil_params *Soil, Reinforcements *reinforcements, bool displayable, Circle *circle);

float determinant( float matrix[3][3], int n);

float Critical_FS_SGA(Soil_params *Soil, Reinforcements *reinforcements, Search_area *search_area, bool search_area_signal, float out_vars[], bool verbose);

float evaluate_individual_FS(Individual *particle, Soil_params *Soil, Reinforcements *reinforcements, bool displayable, Circle *circle);

void Find_section_points(Circle *circle, Soil_params *Soil, float xstart, float xend,float Surface[][2],float Natural_Surface[][2]);

void Fill_opt_SD_zone_angle(Soil_params *Soil, Reinforcements* reinforcements, Search_area *search_area);

#endif /* UTILS_H_ */
