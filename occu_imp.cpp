#include<iostream>
#include<iomanip>
#include<fstream>
#include<math.h>
#include"setup.h"

using namespace std;
void occu_imp(void)
{
	cout << "Occupation of impurity:" << endl;
	double sum_up,sum_down,sumn_up=0,sumn_down=0;
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
	cout << "  occupation_up:   " << scientific << setprecision(15) << setw(25) << sumn_up << endl;
	cout << "  occupation_down: " << scientific << setprecision(15) << setw(25) << sumn_down << endl;
	cout << "  occupation     : " << scientific << setprecision(15) << setw(25) << sumn_up+sumn_down << endl;
	cout << "Time leaved:    ";date_time();cout << endl;
	f_occu.close();
}
