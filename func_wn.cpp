#include<iomanip>
#include<iostream>
#include<fstream>
#include<math.h>
#include<omp.h>
#include<cstdlib>
#include"setup.h"
using namespace std;
double ** exp_z;
void func_wn(void)
{
	cout << "func_wn():  ";date_time();cout << endl;
	double exp_Z(int n,int l);
	exp_z=new double * [N_max];
	for (int n=0;n<n0;n++){
		exp_z[n]=new double [1];
	}
	for (int n=n0;n<N_max-1;n++){
		exp_z[n]=new double [(dim_dot-1)*num_kept];
	}
	for (int n=N_max-1;n<N_max;n++){
		exp_z[n]=new double [num_basis[n]];
	}
	double * wn=new double [N_max];
	double sum=0;
	ofstream fwn("wn.dat");
	for (int n=n0;n<N_max-1;n++){
	    double suml=0;
		for (int l=num_kept;l<num_basis[n];l++){
			exp_z[n][l-num_kept]=exp_Z(n,l);//!!!
		    suml=suml+exp_z[n][l-num_kept];//!!!
		}
		wn[n]=pow(4,N_max-1-n)*suml;
		fwn << setw(5) << n << setw(25) << scientific << setprecision(15) << wn[n] << endl;
		sum=sum+wn[n];
        cout << "  " << setw(5) << n << "  sum of wn= " << setprecision(16) << scientific << sum << endl;
	}
	for (int n=N_max-1;n<N_max;n++){
	    double suml=0;
		for (int l=0;l<num_basis[n];l++){
			exp_z[n][l]=exp_Z(n,l);
		    suml=suml+exp_z[n][l];
		}
		wn[n]=pow(4,N_max-1-n)*suml;
		fwn << setw(5) << n << setw(25) << scientific << setprecision(15) << wn[n] << endl;
		sum=sum+wn[n];
        cout << "  " << setw(5) << n << "  sum of wn= " << setprecision(16) << scientific << sum << endl;
	}
	delete [] wn;
	cout << "Time leaved:    ";date_time();cout << endl;
}

double exp_Z(int n,int l)
{
	double ans;
	double sumn=0;
	for (int n1=n0;n1<N_max-1;n1++){
		if (-1.0*Beta*(E_GS[n1]*pow(Lambda,-1.0*(n1-1-1)/2.0)-E_GS[n]*pow(Lambda,-1.0*(n-1-1)/2.0)) > 500){
			ans=0;
			return ans;
		}else{
	        double suml1=0;
		    for (int l1=num_kept;l1<num_basis[n1];l1++){
		    	suml1=suml1+exp(-1.0*Beta*pow(Lambda,-1.0*(n1-1-1)/2.0)*eigen[n1][l1].eig_val_relat);
		    }
		    sumn=sumn+pow(4,N_max-1-n1)*exp(-1.0*Beta*(E_GS[n1]*pow(Lambda,-1.0*(n1-1-1)/2.0)-E_GS[n]*pow(Lambda,-1.0*(n-1-1)/2.0)))*suml1;
		}
	}
	for (int n1=N_max-1;n1<N_max;n1++){
		if (-1.0*Beta*(E_GS[n1]*pow(Lambda,-1.0*(n1-1-1)/2.0)-E_GS[n]*pow(Lambda,-1.0*(n-1-1)/2.0)) > 500){
			//cout << n << "  " << n1 << "  " << -1.0*Beta*(E_GS[n1]*pow(Lambda,-1.0*(n1-1-1)/2.0)-E_GS[n]*pow(Lambda,-1.0*(n-1-1)/2.0)) << endl;
			ans=0;
			return ans;
		}else{
	        double suml1=0;
		    for (int l1=0;l1<num_basis[n1];l1++){
		    	suml1=suml1+exp(-1.0*Beta*pow(Lambda,-1.0*(n1-1-1)/2.0)*eigen[n1][l1].eig_val_relat);
		    }
		    sumn=sumn+pow(4,N_max-1-n1)*exp(-1.0*Beta*(E_GS[n1]*pow(Lambda,-1.0*(n1-1-1)/2.0)-E_GS[n]*pow(Lambda,-1.0*(n-1-1)/2.0)))*suml1;
		}
	}
	ans=exp(-1.0*Beta*pow(Lambda,-1.0*(n-1-1)/2.0)*eigen[n][l].eig_val_relat)/sumn;
	return ans;
}
