#include<fstream>
#include"setup.h"

int dim_imp=4;
int dim_dot=4;
int num_kept=512;
int N;
int n0;//number of step starting discarding states.num_basis[n0]>num_kept
int quant_tmp=1000;
double Beta_bar=0.6;//unit same as Beta!
double U=0.1;// in unit of D
double Ed;//in unite of D
double Lambda=2.0;
double alpha=0.69;
double temperature=1.0e-9;//k_B*T/D
double omega0=0.5*temperature;
bool smear=false;
bool unsmear=true;
bool Q=false;
bool Q_Sz=false;
bool N_up_N_down=true;
