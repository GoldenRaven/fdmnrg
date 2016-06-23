#include<iomanip>
#include<iostream>
#include<fstream>
#include<math.h>
#include<omp.h>
#include"setup.h"
using namespace std;
double Z;
void func_wn(void)
{
	cout << "func_wn():" << endl;
	double * wn=new double [N_max+2];
	double * Zn=new double [N_max+2];
	for (int n=0;n<N_max+2;n++){
		Zn[n]=0;
		wn[n]=0;
	}
	for (int n=n0;n<N_max+1;n++){//partition function: \sum_{n_0}^{N_max} Zn
		for (int s=num_kept;s<num_basis[n];s++){
			Zn[n]=Zn[n]+exp(-1.0*pow(Lambda,-(n-1-1)/2.0)*Beta*eigen[n][s].eig_val);
		    //cout << "Beta*eigen[n][s].eig_val= " << Beta*eigen[n][s].eig_val << endl;
		}
		//cout << "Zn= " << Zn[n] << endl;
		Z=Z+Zn[n]*pow(dim_dot,N_max-(n-1));
	}
	for (int s=0;s<num_basis[N_max+1];s++){//Z_{N_max+1}
		Zn[N_max+1]=Zn[N_max+1]+exp(-1.0*pow(Lambda,-(N_max+1-1-1)/2.0)*Beta*eigen[N_max+1][s].eig_val);
	    //cout << "sum in N_max+1= " << sum << endl;
	}
	Z=Z+Zn[N_max+1];
	cout << "  " << "partition Z= " << Z << endl;
	double sum=0;
	ofstream fwn("wn.dat");
	for (int i=0;i<N_max+2;i++){
		wn[i]=pow(dim_dot,N_max+1-i)*Zn[i]/Z;
		fwn << setw(3) << i << setw(20) << wn[i] << endl;
		sum=sum+wn[i];
		cout << "  sum of wn= " << setprecision(16) << scientific << sum << endl;
	}
	delete [] wn;
	delete [] Zn;
	cout << "  ";cout << "Time leaved:    ";date_time();
}
