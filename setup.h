#ifndef SETUP_H_
#define SETUP_H_

extern int dim_imp;//=4;
extern int dim_dot;//=4;
extern int num_kept;//=16;
extern int N;
extern int N_max;//=3;//200; No. of dots added.
extern int n0;//number of step starting discarding states.num_basis_new[n0]>num_kept
extern int quant_tmp;
extern double Ed;//=-2.0;
extern double U;//=1.0;
extern double Lambda;//=5.0;
extern double Beta;
extern double Beta_bar;
extern double temperature;
extern double alpha;//width in log-Gaussian
extern double omega0;//width in Gaussian
extern int * num_basis;//[N_max+2].
extern int * num_eigen_kept;//[N_max+2]
extern double * E_GS;
extern double ** exp_z;
extern double *** vect;
extern double *** c_dag_up_basis;
extern double *** c_dag_down_basis;
extern double *** c_up_basis;
extern double *** c_down_basis;
extern double *** c_dag_up_eigen;
extern double *** c_dag_down_eigen;
extern double *** c_up_eigen;
extern double *** c_down_eigen;
extern bool smear;
extern bool Q;
extern bool Q_Sz;
extern bool N_up_N_down;

struct EIGEN_STATE{
	int sort,k,n,quant_num_totalnum,quant_num_upnum,quant_num_downnum;//dim_imp;
	double eig_val,eig_val_relat;
	double * eigen_vect;
};
struct BASIS{
	int k,j,n,quant_num_totalnum,quant_num_upnum,quant_num_downnum;
};
struct NEW_BASIS{
	int sort,k,j,n,quant_num_totalnum,quant_num_upnum,quant_num_downnum;
};

extern EIGEN_STATE ** eigen;
extern BASIS ** basis_ordered;
extern NEW_BASIS *** basis_kj;

void date_time(void);
#endif
