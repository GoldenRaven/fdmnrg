#include"setup.h"

int dim_imp;
int dim_dot;
int num_kept;
int N;
int n0;//number of step starting discarding states.num_basis[n0]>num_kept
int quant_tmp=1000;
double Beta_bar;//unit same as Beta!
double U;// in unit of D
double Ed;//in unite of D
double Lambda;
double alpha;
double temperature;//k_B*T/D
double omega0=temperature;
bool smear;
bool unsmear;
bool Q;
bool Q_Sz;
bool N_up_N_down;
