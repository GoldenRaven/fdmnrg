#include"struct.h"
int job_id=1;
int dim_imp=4;
int dim_dot=4;
int num_kept=512;
int N;
int N_max=0;
int n0;//number of step starting discarding states.num_basis[n0]>num_kept
int quant_tmp=1000;
double Beta_bar=0.69;//unit same as Beta!
double U;// in unit of D
double Ed_up;//in unite of D
double Ed_down;//in unite of D
double Lambda=1.8;
double alpha=0.4;
double temperature;//k_B*T/D
double Beta=0.69;
double omega0=temperature;
double tc=0;
double td=1;
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
bool smear=0;
bool unsmear=1;
bool Q=0;
bool Q_Sz=0;
bool N_up_N_down=1;
bool occupation=1;
bool imp_dos=1;
bool stm_dos=0;
bool eig=0;
double Omega=0;
char * smooth="newsc";
