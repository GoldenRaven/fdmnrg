#ifndef STRUCT_H
#define STRUCT_H
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
#endif
