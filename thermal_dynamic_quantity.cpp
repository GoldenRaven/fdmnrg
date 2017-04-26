#include<iostream>
#include<iomanip>
#include<fstream>
#include<math.h>
#include"setup.h"

double sum_up,sum_down,sumn_up,sumn_down,sum_E,sumn_E,sum_Z,sumn_Z,lnZ;
using namespace std;
void occu_imp(void)
{
	cout << "Occupation of impurity:" << endl;
	sumn_up=0;
	sumn_down=0;
	for (int n=n0;n<N_max-1;n++){
		sum_up=0;
		sum_down=0;
		for (int l=num_kept;l<num_basis[n];l++){
			sum_up=sum_up+exp_z[n][l-num_kept]*occu_imp_up_eigen[n][l][l];
			sum_down=sum_down+exp_z[n][l-num_kept]*occu_imp_down_eigen[n][l][l];
		}
		sumn_up=sumn_up+pow(4,N_max-1-n)*sum_up;
		sumn_down=sumn_down+pow(4,N_max-1-n)*sum_down;
	}
	for (int n=N_max-1;n<N_max;n++){
		sum_up=0;
		sum_down=0;
		for (int l=num_kept;l<num_basis[n];l++){
			sum_up=sum_up+exp_z[n][l-num_kept]*occu_imp_up_eigen[n][l][l];
			sum_down=sum_down+exp_z[n][l-num_kept]*occu_imp_down_eigen[n][l][l];
		}
		sumn_up=sumn_up+pow(4,N_max-1-n)*sum_up;
		sumn_down=sumn_down+pow(4,N_max-1-n)*sum_down;
	}
	ofstream f_occu("occu.dat", ios_base::out | ios_base::app);
	//up_occupation    down_occupation    up_occupation+down_occupation
	f_occu << scientific << setprecision(15) << setw(25) << sumn_up << setw(25) << sumn_down << setw(25) << sumn_up+sumn_down << endl;
	cout << setw(35) << "  occupation_up:   " << scientific << setprecision(15) << setw(25) << sumn_up << endl;
	cout << setw(35) << "  occupation_down: " << scientific << setprecision(15) << setw(25) << sumn_down << endl;
	cout << setw(35) << "  occupation     : " << scientific << setprecision(15) << setw(25) << sumn_up+sumn_down << endl;
	cout << "Time leaved:    ";date_time();cout << endl;
	f_occu.close();
}

double inner_energy(int diff)
{
	cout << "inner_energy: " << endl;
	sumn_E=0;
	for (int n=n0;n<N_max-1;n++){
		sum_E=0;
		for (int l=num_kept;l<num_basis[n];l++){
			sum_E=sum_E+exp_z[n][l-num_kept]*(eigen[n][l].eig_val_relat+E_GS[n])*pow(Lambda,-1.0*(n-1-diff)/2.0);
		}
		sumn_E=sumn_E+pow(4,N_max-1-n)*sum_E;
	}
	for (int n=N_max-1;n<N_max;n++){
		sum_E=0;
		for (int l=num_kept;l<num_basis[n];l++){
			sum_E=sum_E+exp_z[n][l-num_kept]*(eigen[n][l].eig_val_relat+E_GS[n])*pow(Lambda,-1.0*(n-1-diff)/2.0);
		}
		sumn_E=sumn_E+pow(4,N_max-1-n)*sum_E;
	}
	cout << setw(35) << "  inner energy : " << scientific << setprecision(15) << setw(25) << sumn_E << endl;
	cout << "Time leaved:    ";date_time();cout << endl;
	return sumn_E;
}

double entropy_kB(int diff)
{
	//cout << "entropy: " << endl;
	sumn_E=0;
	sumn_Z=0;
	for (int n=n0;n<N_max-1;n++){
		sum_E=0;
		sum_Z=0;
		for (int l=num_kept;l<num_basis[n];l++){
			sum_E=sum_E+exp_z[n][l-num_kept]*(eigen[n][l].eig_val_relat+E_GS[n])*pow(Lambda,-1.0*(n-1-diff)/2.0);
			sum_Z=sum_Z+exp(-1.0*Beta*((eigen[n][l].eig_val_relat+E_GS[n])*pow(Lambda,-1.0*(n-1-diff)/2.0)-E_GS[N_max-1]*pow(Lambda,-1.0*(N_max-1-1-diff)/2.0)));
			//cout << n << "   " << l << "  " << eigen[n][l].eig_val*pow(Lambda,-1.0*(n-1-diff)/2.0) << "    " << (eigen[n][l].eig_val_relat+E_GS[n])*pow(Lambda,-1.0*(n-1-diff)/2.0) << endl;
			//cout << n << "   " << l << "   "  << eigen[n][l].eig_val*pow(Lambda,-1.0*(n-1-diff)/2.0) << "     " << E_GS[N_max-1]*pow(Lambda,-1.0*(N_max-1-1-diff)/2.0) << "   1xxxxx    " << exp(-1.0*Beta*(eigen[n][l].eig_val*pow(Lambda,-1.0*(n-1-diff)/2.0)-E_GS[N_max-1]*pow(Lambda,-1.0*(N_max-1-1-diff)/2.0))) << endl;
			//cout << "Dot  " << n << "  N_up_" << left << setw(5) << scientific << eigen[n][l].quant_num_upnum << "  N_down_" << left << setw(5) << scientific << eigen[n][l].quant_num_downnum << setw(25) << setprecision(15) << eigen[n][l].eig_val_relat << (eigen[n][l].eig_val_relat+E_GS[n])*pow(Lambda,-1.0*(n-1-diff)/2.0) << "   " << eigen[n][l].k << endl;
		}
		sumn_E=sumn_E+pow(4,N_max-1-n)*sum_E;
		sumn_Z=sumn_Z+pow(4,N_max-1-n)*sum_Z;
	}
	for (int n=N_max-1;n<N_max;n++){
		sum_E=0;
		sum_Z=0;
		for (int l=num_kept;l<num_basis[n];l++){
			sum_E=sum_E+exp_z[n][l-num_kept]*(eigen[n][l].eig_val_relat+E_GS[n])*pow(Lambda,-1.0*(n-1-diff)/2.0);
			sum_Z=sum_Z+exp(-1.0*Beta*((eigen[n][l].eig_val_relat+E_GS[n])*pow(Lambda,-1.0*(n-1-diff)/2.0)-E_GS[N_max-1]*pow(Lambda,-1.0*(N_max-1-1-diff)/2.0)));
			//cout << n << "   " << l << "   "  << eigen[n][l].eig_val*pow(Lambda,-1.0*(n-1-diff)/2.0) << "     " << E_GS[N_max-1]*pow(Lambda,-1.0*(N_max-1-1-diff)/2.0) << "   1xxxxx    " << exp(-1.0*Beta*(eigen[n][l].eig_val*pow(Lambda,-1.0*(n-1-diff)/2.0)-E_GS[N_max-1]*pow(Lambda,-1.0*(N_max-1-1-diff)/2.0))) << endl;
		}
		sumn_E=sumn_E+pow(4,N_max-1-n)*sum_E;
		sumn_Z=sumn_Z+pow(4,N_max-1-n)*sum_Z;
	}
	//cout << -1.0*Beta*E_GS[N_max-1]*pow(Lambda,-1.0*(N_max-1-1-diff)/2.0) << endl;
	//cout << " xxxxxxxxx    " << sumn_Z << endl;
	//cout << log(sumn_Z) << " xxxxxxxxx" << endl;
	lnZ=-1.0*Beta*E_GS[N_max-1]*pow(Lambda,-1.0*(N_max-1-1-diff)/2.0)+log(sumn_Z);
	//cout << setw(35) << "  entropy : " << scientific << setprecision(15) << setw(25) << Beta*sumn_E + lnZ << endl;
	//cout << "Time leaved:    ";date_time();cout << endl;
	return Beta*sumn_E + lnZ;
}
