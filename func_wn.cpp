#include<iomanip>
#include<iostream>
#include<fstream>
#include<math.h>
#include<omp.h>
#include<cstdlib>
#include"setup.h"
using namespace std;
double Z;
void func_wn(void)
{
	cout << "func_wn():  " << endl;date_time();
	double * wn=new double [N_max];
	double * Zn=new double [N_max];
	for (int n=0;n<N_max;n++){
		Zn[n]=0;
		wn[n]=0;
	}
	for (int n=n0;n<N_max-1;n++){//partition function: \sum_{n_0}^{N_max} Zn
		//if (n==6) exit(0);
		for (int s=num_kept;s<num_basis[n];s++){
			//cout << setw(20) << eigen[n][s].eig_val << setw(20) << exp(-1.0*pow(Lambda,-(n-1-1)/2.0)*Beta*eigen[n][s].eig_val) << endl;
			Zn[n]=Zn[n]+exp(-1.0*pow(Lambda,-(n-1-1)/2.0)*Beta*eigen[n][s].eig_val);
		    //cout << "Beta*eigen[n][s].eig_val= " << Beta*eigen[n][s].eig_val << endl;
		}
		//cout << "n=  " << n << "   sum exp=  " << Zn[n] << endl;
		Zn[n]=Zn[n]*pow(dim_dot,N_max-1-n);
		//cout << "n=  " << n << "exp*pow=  " << Zn[n] << endl;
		Z=Z+Zn[n];
	}
	for (int s=0;s<num_basis[N_max-1];s++){//Z_{N_max-1}
		Zn[N_max-1]=Zn[N_max-1]+exp(-1.0*pow(Lambda,-(N_max-1-1-1)/2.0)*Beta*eigen[N_max-1][s].eig_val);
	}
		//cout << "n=  " << N_max-1 << "sum exp=  " << Zn[N_max-1] << endl;
		//cout << "n=  " << N_max-1 << "exp*pow=  " << Zn[N_max-1] << endl;
	Z=Z+Zn[N_max-1];
	cout << "  " << "partition Z= " << Z << endl;
	double sum=0;
	ofstream fwn("wn.dat");
	for (int i=n0;i<N_max;i++){
		wn[i]=Zn[i]/Z;
		fwn << setw(3) << i << setw(20) << wn[i] << endl;
		sum=sum+wn[i];
		cout << "  sum of wn= " << setprecision(16) << scientific << sum << endl;
	}
	delete [] wn;
	delete [] Zn;
	cout << "  ";cout << "Time leaved:    ";date_time();
}
