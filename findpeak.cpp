#include<iostream>
#include<fstream>
#include<cstdlib>
#include<math.h>
#include"setup.h"
using namespace std;
double peak(double * x, double * y, double x_min, double x_max, int dim)
{
	int dim_new=0;
	int index;
	double y_max,freq;
	for (int i=0;i<dim;i++){
		if (x[i] >= x_min && x[i] <= x_max) dim_new++;
	}
	double * x_new = new double [dim_new];
	double * y_new = new double [dim_new];
	int j=0;
	for (int i=0;i<dim;i++){
		if (x[i] >= x_min && x[i] <= x_max){
		   x_new[j]=x[i];
		   y_new[j]=y[i];
		   j++;
		}
	}
	//sorting of x.
	for (int i=1;i<dim_new-1;i++){
		if (x_new[dim_new-1] >= x_new[i-1] && x_new[dim_new-1] <= x_new[i]){
			xxxxxxx;
		}
	}
	for (int i=0;i<dim_new-1;i++){
		y_max=y_new[0];
		index=0;
		if (y_max <= y_new[i+1]){
			y_max=y_new[i+1];
			index=i+1;
		}
	}
	freq=(x_new[index-1]+x_new[index+1])/2.0;
	if (fabs(freq-x_new[index]) < x_error) cout << "  Warning! x" << endl;
	return freq;
}
