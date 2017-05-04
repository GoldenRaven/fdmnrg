#include<iomanip>
#include<iostream>
#include<fstream>
#include<math.h>
#include<omp.h>
#include<cstdlib>
#include"setup.h"
using namespace std;
double ** exp_z;
void func_wn(int diff)
{
	cout << "func_wn():  ";date_time();cout << endl;
	double exp_Z(int n,int l, int diff);
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
	ofstream fwn("wn.dat",ios_base::out | ios_base::app);
	for (int n=n0;n<N_max-1;n++){
	    double suml=0;
		for (int l=num_kept;l<num_basis[n];l++){
			//cout << n << "  " << l << "  " << endl;
			//cout << exp_Z(n,l) << "    " << n << "    " << l << endl;
			exp_z[n][l-num_kept]=exp_Z(n,l,diff);//!!!
			//cout << "xxxxxxxxxxxx" << endl;
		    suml=suml+exp_z[n][l-num_kept];//!!!
			//cout << "xxxxxxxxxxxx" << endl;
		}
		wn[n]=pow(dim_dot,N_max-1-n)*suml;
		fwn << setw(5) << n << setw(25) << scientific << setprecision(15) << wn[n] << endl;
		sum=sum+wn[n];
        cout << "  " << setw(5) << n << "  sum of wn= " << setprecision(16) << scientific << sum << endl;
	}
	for (int n=N_max-1;n<N_max;n++){
	    double suml=0;
		for (int l=0;l<num_basis[n];l++){
			exp_z[n][l]=exp_Z(n,l,diff);
		    suml=suml+exp_z[n][l];
		}
		wn[n]=pow(dim_dot,N_max-1-n)*suml;
		fwn << setw(5) << n << setw(25) << scientific << setprecision(15) << wn[n] << endl;
		sum=sum+wn[n];
        cout << "  " << setw(5) << n << "  sum of wn= " << setprecision(16) << scientific << sum << endl;
	}
	fwn.close();
	delete [] wn;
	cout << "Time leaved:    ";date_time();cout << endl;
}

double exp_Z(int n,int l,int diff)
{
	double ans;
	double sum=0;
	for (int n1=n0;n1<N_max-1;n1++){
		for (int l1=num_kept;l1<num_basis[n1];l1++){
			double expo=0;
			expo=-1.0*Beta*((E_GS[n1]+eigen[n1][l1].eig_val_relat)*pow(Lambda,-1.0*(n1-1-diff)/2.0)-(E_GS[n]+eigen[n][l].eig_val_relat)*pow(Lambda,-1.0*(n-1-diff)/2.0));
			if (expo > 709) return 0;
			sum=sum+pow(dim_dot,N_max-1-n1)*exp(expo);
		}
	}
	for (int n1=N_max-1;n1<N_max;n1++){
		for (int l1=num_kept;l1<num_basis[n1];l1++){
			double expo=0;
			expo=-1.0*Beta*((E_GS[n1]+eigen[n1][l1].eig_val_relat)*pow(Lambda,-1.0*(n1-1-diff)/2.0)-(E_GS[n]+eigen[n][l].eig_val_relat)*pow(Lambda,-1.0*(n-1-diff)/2.0));
			if (expo > 709) return 0;
			sum=sum+pow(dim_dot,N_max-1-n1)*exp(expo);
		}
	}
	ans=1.0/sum;
	return ans;
}
