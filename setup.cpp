#include"struct.h"
int job_id;
int dim_imp;
int dim_dot;
int num_kept;
int N;
int N_max;
int n0;//number of step starting discarding states.num_basis[n0]>num_kept
int quant_tmp=1000;
double Beta_bar;//unit same as Beta!
double U;// in unit of D
double Ed_up;//in unite of D
double Ed_down;//in unite of D
double Lambda;
double alpha;
double temperature;//k_B*T/D
double Beta;
double omega0=temperature;
double tc;
double td;
double *** c_up_eigen;
double *** c_down_eigen;
double *** occu_imp_up_eigen;
double *** occu_imp_down_eigen;
double *** stm_f_up_eigen;
double *** stm_f_down_eigen;
double * E_GS;
int * num_basis;
int * num_eigen_kept;//[N_max+2]
BASIS ** basis_ordered;
NEW_BASIS *** basis_kj;
EIGEN_STATE ** eigen;
bool smear;
bool unsmear;
bool Q;
bool Q_Sz;
bool N_up_N_down;
bool occupation;
bool imp_dos;
bool stm_dos;
