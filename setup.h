#ifndef SETUP_H_
#define SETUP_H_
#include"struct.h"

extern int job_id;
extern int dim_imp;//=4;
extern int dim_dot;//=4;
extern int num_kept;//=16;
extern int N;
extern int N_max;//=3;//200; No. of dots added.
extern int n0;//number of step starting discarding states.num_basis_new[n0]>num_kept
extern int quant_tmp;
extern double Ed_up;//=-2.0;
extern double Ed_down;//=-2.0;
extern double U;//=1.0;
extern double Lambda;//=5.0;
extern double Beta;
extern double Beta_bar;
extern double temperature;
extern double alpha;//width in log-Gaussian
extern double omega0;//width in Gaussian
extern double tc;//stm coupling with c electron
extern double td;//stm coupling with d electron
extern int * num_basis;//[N_max+2].
extern int * num_eigen_kept;//[N_max+2]
extern double * E_GS;
extern double ** exp_z;
extern double *** vect;
extern double *** c_up_eigen;
extern double *** c_down_eigen;
extern double *** occu_imp_up_eigen;
extern double *** occu_imp_down_eigen;
extern double *** stm_f_up_eigen;
extern double *** stm_f_down_eigen;
extern bool smear;
extern bool unsmear;
extern bool Q;
extern bool Q_Sz;
extern bool N_up_N_down;
extern bool occupation;
extern bool imp_dos;
extern bool stm_dos;
extern bool eig;
extern double Omega;
extern char * smooth;
extern EIGEN_STATE ** eigen;
extern BASIS ** basis_ordered;
extern NEW_BASIS *** basis_kj;

void date_time(void);
#endif
