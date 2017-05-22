#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include<math.h>
#include<mkl.h>
#include<omp.h>
#include<string.h>
#include"setup.h"
#include"struct.h"
using namespace std;

void genoutput(void);
void occu_imp_operator(int);
int func_delta(int a, int b){
	if (a==b){
		return 1;
	}else{
		return 0;
	}
}

void iterative_dia_band(void)
{
	int num_block;
	int * num_basis_block;
	BASIS ** block;
	double *** H_bij;
	double ** temp1;
	double ** temp2;
	double ** temp3;
	double ** temp4;
	double *** vect;
	NEW_BASIS *** basis_kj;
	cout << "  ";cout << "Iteriterative_dia_band():    ";date_time();cout << endl;

	ifstream f_input("input_band");
	f_input >> temperature ;
	f_input >> Lambda      ;
	f_input >> num_kept    ;
	f_input >> dim_imp     ;
	f_input >> dim_dot     ;
	f_input >> Beta_bar    ;
	f_input >> Q           ;
	f_input >> Q_Sz        ;
	f_input >> N_up_N_down ;

	Beta=1.0/temperature;//1.0/(k_B*T/D)
	N = int(1-(2.0*log(Beta_bar/Beta)/log(Lambda)));//site -1,0,...,N_max
	N_max = N+15;//!!!!!!!!!!!!!!!!!!!!

    double coupling_imp_dot_up;
    double coupling_imp_dot_down;
	double * pe_up = new double [N_max+1];
	double * pe_down = new double [N_max+1];
	double * ptn_up = new double [N_max+1];
	double * ptn_down = new double [N_max+1];
	E_GS=new double [N_max];
	num_basis=new int [N_max+1];
    num_eigen_kept=new int [N_max+1];
	ifstream fin_chain("chain_band.dat");
	cout << "  chain parameter:" << endl;
	cout << setw(13) << "  #site" << setw(29) << "   tn_up" << setw(29) << "  tn_down" << setw(29) << "    en_up" << setw(29) << "      en_down" << endl;
	for (int i=-1;i<N_max+1;i++){//from second line: coupling, t0, t1,...
		double temp;
		fin_chain >> scientific;
		if (i==-1){
		    fin_chain >> temp >> coupling_imp_dot_up >> coupling_imp_dot_down;
			cout << "  " << scientific << setw(5) << i << setw(30) << setprecision(20) << coupling_imp_dot_up << setw(30) << setprecision(20) << coupling_imp_dot_down << endl;
		}else{
			//fin_chain >> setprecision(20) >> temp >> setprecision(20) >> ptn_up[i] >> setprecision(20) >> ptn_down[i] >> setprecision(20) >> pe_up[i] >> pe_down[i] >> setprecision(20) >> temp;
			fin_chain >> setprecision(20) >> temp >> setprecision(20) >> ptn_up[i] >> setprecision(20) >> ptn_down[i] >> setprecision(20) >> pe_up[i] >> pe_down[i];
			cout << "  " << scientific << setw(5) << i << setw(30) << setprecision(20) << ptn_up[i] << setw(30) << setprecision(20) << ptn_down[i] << setw(30) << setprecision(20) << pe_up[i] << setw(30) << setprecision(20) << pe_down[i] << endl;
		}
	}
	N_max=N_max-1;
	Ed_up=pe_up[0];
	Ed_down=pe_down[0];
	U=0;
	genoutput();
	cout << "   " << " n  " << "  num_basis[n]  " << "  num_eigen_kept[n]  " << "  n0 " << endl;
	for (int n=0;n<N_max+1;n++){
		if (n<10){
	        if (dim_imp*pow(dim_dot,n)<num_kept){
	            num_basis[n]=dim_imp*pow(dim_dot,n); //n=0,1,...
			    num_eigen_kept[n]=dim_imp*pow(dim_dot,n);
	        }else{
	            if (dim_imp*pow(dim_dot,n) >= dim_dot*num_kept){
	            	num_basis[n]=dim_dot*num_kept;
	            }else{
	                num_basis[n]=dim_imp*pow(dim_dot,n); //n=0,1,...
				}
	        	num_eigen_kept[n]=num_kept;
			}
			if (num_basis[n]<=num_kept){
				n0=n+1;
			}
		}else{
	        num_basis[n]=num_basis[9]; //                 9 means big enough.
	        num_eigen_kept[n]=num_eigen_kept[9];//n=0,1,...    9 means big enough.
		}
	}
	//for (int n=0;n<N_max+1;n++){
	for (int n=0;n<8;n++){
		cout << "    " << setw(3) << n << "     " << setw(7) << num_basis[n] << "             " << setw(7) << num_eigen_kept[n] << "       " << setw(3) << n0 << endl;
	}
	cout << "  step starting to discard state: " << n0 << endl;
    int quant_num_totalnum_dot[4]={0,1,1,2};//{0,up,down,updown}
    int quant_num_upnum_dot[4]={0,1,0,1};
    int quant_num_downnum_dot[4]={0,0,1,1};
	double c_up_imp[4][4]={{0,1,0,0},{0,0,0,0},{0,0,0,1},{0,0,0,0}};
	double c_down_imp[4][4]={{0,0,1,0},{0,0,0,-1},{0,0,0,0},{0,0,0,0}};
	double c_dag_up_imp[4][4]={{0,0,0,0},{1,0,0,0},{0,0,0,0},{0,0,1,0}};
	double c_dag_down_imp[4][4]={{0,0,0,0},{0,0,0,0},{1,0,0,0},{0,-1,0,0}};
	double c_up_dot[4][4]={{0,1,0,0},{0,0,0,0},{0,0,0,1},{0,0,0,0}};
	double c_down_dot[4][4]={{0,0,1,0},{0,0,0,-1},{0,0,0,0},{0,0,0,0}};
	double c_dag_up_dot[4][4]={{0,0,0,0},{1,0,0,0},{0,0,0,0},{0,0,1,0}};
	double c_dag_down_dot[4][4]={{0,0,0,0},{0,0,0,0},{1,0,0,0},{0,-1,0,0}};
	double c_up_eigen_0[4][4];
	double c_down_eigen_0[4][4];
	double c_dag_up_eigen_0[4][4];
	double c_dag_down_eigen_0[4][4];
	for (int i=0;i<4;i++){
		for (int j=0;j<4;j++){
			c_up_eigen_0[i][j]=0;
			c_down_eigen_0[i][j]=0;
			c_dag_up_eigen_0[i][j]=0;
			c_dag_down_eigen_0[i][j]=0;
		}
	}
    basis_ordered=new BASIS * [N_max];
	basis_kj=new NEW_BASIS ** [N_max];
    eigen=new EIGEN_STATE * [N_max];
	vect=new double ** [N_max];
	for (int n=0;n<1;n++){//basis_kj[0] shouldn't be used!
	    basis_kj[n]=new NEW_BASIS * [dim_dot];
	    for (int i=0;i<dim_dot;i++){
	        basis_kj[n][i]=new NEW_BASIS [dim_dot];
	    }
	}
	for (int n=1;n<N_max;n++){//basis_kj[n]
	    basis_kj[n]=new NEW_BASIS * [num_eigen_kept[n-1]];
		for (int i=0;i<num_eigen_kept[n-1];i++){
		    basis_kj[n][i]=new NEW_BASIS [dim_dot];
		}
	}
	for (int n=0;n<N_max;n++){
        basis_ordered[n]=new BASIS [num_basis[n]];// 二维可变数组.储存基矢排序.4,16,64,
        eigen[n]=new EIGEN_STATE [num_basis[n]];// eigen[0][k] --> |k;-1>
		vect[n]=new double * [num_basis[n]];
		for (int i=0;i<num_basis[n];i++){
			vect[n][i]=new double [num_basis[n]];
			for (int j=0;j<num_basis[n];j++){
				vect[n][i][j]=0;
			}
		}
	}
	BASIS basis_temp0={0,1,-1,0,0,0};basis_ordered[0][0]=basis_temp0;//|0>=|1,0;-1>
	BASIS basis_temp1={0,2,-1,1,1,0};basis_ordered[0][1]=basis_temp1;//|up>=|2,0;-1>
	BASIS basis_temp2={0,3,-1,1,0,1};basis_ordered[0][2]=basis_temp2;//|down>=|3,0;-1>
	BASIS basis_temp3={0,4,-1,2,1,1};basis_ordered[0][3]=basis_temp3;//|up down>=|4,0;-1>
	vect[0][0][0]=1;vect[0][0][1]=0;vect[0][0][2]=0;vect[0][0][3]=0;
	vect[0][1][0]=0;vect[0][1][1]=1;vect[0][1][2]=0;vect[0][1][3]=0;
    vect[0][2][0]=0;vect[0][2][1]=0;vect[0][2][2]=1;vect[0][2][3]=0;
    vect[0][3][0]=0;vect[0][3][1]=0;vect[0][3][2]=0;vect[0][3][3]=1;
	double temp[4]={0,Ed_up,Ed_down,Ed_up+Ed_down+U};//sorting of eigenvalue
	for (int j=0;j<num_basis[0]-1;j++){
		if (temp[j]<temp[j+1]){
			double tempp=temp[j];
			temp[j]=temp[j+1];
			temp[j+1]=tempp;
		}
	}
	EIGEN_STATE eig_temp0={1,1,0,1,1,0,(Ed_up)/sqrt(Lambda),(Ed_up-temp[3])/sqrt(Lambda),vect[0][1]};
	EIGEN_STATE eig_temp1={2,2,0,1,0,1,(Ed_down)/sqrt(Lambda),(Ed_down-temp[3])/sqrt(Lambda),vect[0][2]};
	EIGEN_STATE eig_temp2={3,3,0,0,0,0,0,(0-temp[3])/sqrt(Lambda),vect[0][0]};
	EIGEN_STATE eig_temp3={4,4,0,2,1,1,(Ed_up+Ed_down+U)/sqrt(Lambda),(Ed_up+Ed_down+U-temp[3])/sqrt(Lambda),vect[0][3]};
	eigen[0][0]=eig_temp0;//| up     > = |1,0>
	eigen[0][1]=eig_temp1;//| down   > = |2,0>
	eigen[0][2]=eig_temp2;//| 0      > = |3,0>
	eigen[0][3]=eig_temp3;//| updown > = |4,0>
	E_GS[0]=temp[3]/sqrt(Lambda);
	for (int i=0;i<4;i++){
		cout << eigen[0][i].eig_val*sqrt(Lambda) << endl;
	}
	//exit(0);
	temp1=new double * [dim_dot];
	temp2=new double * [dim_dot];
	temp3=new double * [dim_dot];
	temp4=new double * [dim_dot];
	for (int i=0;i<dim_dot;i++){
		temp1[i]=new double [dim_dot];
		temp2[i]=new double [dim_dot];
		temp3[i]=new double [dim_dot];
		temp4[i]=new double [dim_dot];
		for (int j=0;j<dim_dot;j++){
			temp1[i][j]=0;
			temp2[i][j]=0;
			temp3[i][j]=0;
			temp4[i][j]=0;
		}
	}
	for (int i=0;i<num_basis[0];i++){
        for (int j=0;j<num_basis[0];j++){
			for (int k=0;k<num_basis[0];k++){
				temp1[i][j]=temp1[i][j]+eigen[0][i].eigen_vect[k]*c_up_imp[k][j];
				temp2[i][j]=temp2[i][j]+eigen[0][i].eigen_vect[k]*c_down_imp[k][j];
				temp3[i][j]=temp3[i][j]+eigen[0][i].eigen_vect[k]*c_dag_up_imp[k][j];
				temp4[i][j]=temp4[i][j]+eigen[0][i].eigen_vect[k]*c_dag_down_imp[k][j];
			}
		}
	}
	for (int i=0;i<num_basis[0];i++){
        for (int j=0;j<num_basis[0];j++){
			for (int k=0;k<num_basis[0];k++){
				c_up_eigen_0[i][j]=c_up_eigen_0[i][j]+temp1[i][k]*eigen[0][j].eigen_vect[k];
				c_down_eigen_0[i][j]=c_down_eigen_0[i][j]+temp2[i][k]*eigen[0][j].eigen_vect[k];
				c_dag_up_eigen_0[i][j]=c_dag_up_eigen_0[i][j]+temp3[i][k]*eigen[0][j].eigen_vect[k];
				c_dag_down_eigen_0[i][j]=c_dag_down_eigen_0[i][j]+temp4[i][k]*eigen[0][j].eigen_vect[k];
			}
		}
	}
	for (int i=0;i<dim_dot;i++){
		delete [] temp1[i];
		delete [] temp2[i];
		delete [] temp3[i];
		delete [] temp4[i];
	}
	delete [] temp1;
	delete [] temp2;
	delete [] temp3;
	delete [] temp4;
	ofstream f_E_GS("E_GS.dat",ios_base::out|ios_base::app);
	for (int n=1;n<N_max;n++){
		char str0[20],str1[15];
		sprintf(str0,"%d",n);
		strcpy(str1,"_U.dat");
		strcat(str0,str1);
		//ofstream f_U(str0);
		char str02[20],str2[15];
		sprintf(str02,"%d",n);
		strcpy(str2,"_eigenvalue.dat");
		strcat(str02,str2);
		ofstream f_eig_val(str02);
		char str03[20],str3[15];
		sprintf(str03,"%d",n);
		strcpy(str3,"_d_up.dat");
		strcat(str03,str3);
		//ofstream f_d_up(str03);
		char str04[20],str4[15];
		sprintf(str04,"%d",n);
		strcpy(str4,"_d_down.dat");
		strcat(str04,str4);
		//ofstream f_d_down(str04);
		{int i=0;
	    for (int k=0;k<num_eigen_kept[n-1];k++){
	        for (int j=0;j<dim_dot;j++){
	    		basis_ordered[n][i].k=k+1;
	    		basis_ordered[n][i].j=j+1;
	    		basis_ordered[n][i].n=n-1;
	    		basis_ordered[n][i].quant_num_totalnum=eigen[n-1][k].quant_num_totalnum + quant_num_totalnum_dot[j];
	    		basis_ordered[n][i].quant_num_upnum=eigen[n-1][k].quant_num_upnum + quant_num_upnum_dot[j];
	    		basis_ordered[n][i].quant_num_downnum=eigen[n-1][k].quant_num_downnum+quant_num_downnum_dot[j];
	    	    basis_kj[n][k][j].k=k+1;
	    	    basis_kj[n][k][j].j=j+1;
	    	    basis_kj[n][k][j].n=n-1;
				basis_kj[n][k][j].quant_num_totalnum=eigen[n-1][k].quant_num_totalnum+quant_num_totalnum_dot[j];
	    		basis_kj[n][k][j].quant_num_upnum=eigen[n-1][k].quant_num_upnum+quant_num_upnum_dot[j];
	    		basis_kj[n][k][j].quant_num_downnum=eigen[n-1][k].quant_num_downnum+quant_num_downnum_dot[j];
				i++;
	        }
	    }}
		if (Q){
		    for (int i=1;i<num_basis[n];i++){//sorting of quantum number: Q.
		    	for (int ii=1;ii<=num_basis[n]-i;ii++){
		    		if (basis_ordered[n][ii-1].quant_num_totalnum > basis_ordered[n][ii].quant_num_totalnum){
		    			BASIS temp = basis_ordered[n][ii-1];
		    			basis_ordered[n][ii-1]=basis_ordered[n][ii];
		    			basis_ordered[n][ii]=temp;
		    		}
		    	}
		    }
	        num_block=0;
	        for (int i=1;i<=num_basis[n];i++){
		    	basis_kj[n][basis_ordered[n][i-1].k-1][basis_ordered[n][i-1].j-1].sort=i-1;
		    	if (basis_ordered[n][i-1].quant_num_totalnum != basis_ordered[n][i].quant_num_totalnum){
	        		num_block=num_block+1;
	        	}
	        }
	        num_basis_block=new int [num_block];
	        for (int i=0;i<num_block;i++){
	        	num_basis_block[i]=0;
	        }
	        {int ii=0;
	        for (int i=1;i<=num_basis[n];i++){
	        	if (basis_ordered[n][i-1].quant_num_totalnum != basis_ordered[n][i].quant_num_totalnum){
	        		num_basis_block[ii]=i;//-num_basis_block[ii-2];
	        	    ii++;
	        	}
	        }}
		}else if(Q_Sz){
		    for (int i=1;i<num_basis[n];i++){//sorting of quantum number: Q, S_z.
		    	for (int ii=1;ii<=num_basis[n]-i;ii++){
		    		if (quant_tmp*basis_ordered[n][ii-1].quant_num_totalnum+(basis_ordered[n][ii-1].quant_num_upnum-basis_ordered[n][ii-1].quant_num_downnum) > quant_tmp*basis_ordered[n][ii].quant_num_totalnum+(basis_ordered[n][ii].quant_num_upnum-basis_ordered[n][ii].quant_num_downnum)){
		    			BASIS temp = basis_ordered[n][ii-1];
		    			basis_ordered[n][ii-1]=basis_ordered[n][ii];
		    			basis_ordered[n][ii]=temp;
		    		}
		    	}
		    }
	        num_block=0;
	        for (int i=1;i<=num_basis[n];i++){
		    	//cout << basis_ordered[n][i-1].quant_num_totalnum << "  N  " << i-1 << endl;
		    	basis_kj[n][basis_ordered[n][i-1].k-1][basis_ordered[n][i-1].j-1].sort=i-1;
		    	if (quant_tmp*basis_ordered[n][i-1].quant_num_totalnum+(basis_ordered[n][i-1].quant_num_upnum-basis_ordered[n][i-1].quant_num_downnum) != quant_tmp*basis_ordered[n][i].quant_num_totalnum+(basis_ordered[n][i].quant_num_upnum-basis_ordered[n][i].quant_num_downnum)){// Q, Sz
	        		num_block=num_block+1;
	        	}
	        }
	        num_basis_block=new int [num_block];
	        for (int i=0;i<num_block;i++){
	        	num_basis_block[i]=0;
	        }
	        {int ii=0;
	        for (int i=1;i<=num_basis[n];i++){
	        	if (quant_tmp*basis_ordered[n][i-1].quant_num_totalnum+(basis_ordered[n][i-1].quant_num_upnum-basis_ordered[n][i-1].quant_num_downnum) != quant_tmp*basis_ordered[n][i].quant_num_totalnum+(basis_ordered[n][i].quant_num_upnum-basis_ordered[n][i].quant_num_downnum)){
	        		num_basis_block[ii]=i;//-num_basis_block[ii-2];
	        	    ii++;
	        	}
	        }}
		}else if(N_up_N_down){
		    for (int i=1;i<num_basis[n];i++){//sorting of quantum number: up, down.
		    	for (int ii=1;ii<=num_basis[n]-i;ii++){
		    		if (quant_tmp*basis_ordered[n][ii-1].quant_num_upnum+(basis_ordered[n][ii-1].quant_num_downnum) > quant_tmp*basis_ordered[n][ii].quant_num_upnum+(basis_ordered[n][ii].quant_num_downnum)){
		    			BASIS temp = basis_ordered[n][ii-1];
		    			basis_ordered[n][ii-1]=basis_ordered[n][ii];
		    			basis_ordered[n][ii]=temp;
		    		}
		    	}
		    }
	        num_block=0;
	        for (int i=1;i<=num_basis[n];i++){
		    	//cout << basis_ordered[n][i-1].quant_num_totalnum << "  N  " << i-1 << endl;
		    	basis_kj[n][basis_ordered[n][i-1].k-1][basis_ordered[n][i-1].j-1].sort=i-1;
		    	if (quant_tmp*basis_ordered[n][i-1].quant_num_upnum+(basis_ordered[n][i-1].quant_num_downnum) != quant_tmp*basis_ordered[n][i].quant_num_upnum+(basis_ordered[n][i].quant_num_downnum)){// N_up, N_down
	        		num_block=num_block+1;
	        	}
	        }
	        num_basis_block=new int [num_block];
	        for (int i=0;i<num_block;i++){
	        	num_basis_block[i]=0;
	        }
	        {int ii=0;
	        for (int i=1;i<=num_basis[n];i++){
	        	if (quant_tmp*basis_ordered[n][i-1].quant_num_upnum+(basis_ordered[n][i-1].quant_num_downnum) != quant_tmp*basis_ordered[n][i].quant_num_upnum+(basis_ordered[n][i].quant_num_downnum)){
	        		num_basis_block[ii]=i;//-num_basis_block[ii-2];
	        	    ii++;
	        	}
	        }}
		}
		/*
		 *for (int i=0;i<num_basis[n];i++){
		 *    cout << "basis  " << basis_ordered[n][i].quant_num_totalnum << "  " << basis_ordered[n][i].quant_num_upnum << "  " << basis_ordered[n][i].quant_num_downnum << endl;
		 *}
		 */
		for (int i=num_block-1;i>0;i--){
		    num_basis_block[i]=num_basis_block[i]-num_basis_block[i-1];
		}//1,4,6,4,1.
		/*
		 *for (int i=1;i<=num_block;i++){
		 *    cout << num_basis_block[i-1] << "  xxx" << endl;
		 *}
		 */
	    block=new BASIS * [num_block];
	    for (int i=0;i<num_block;i++){
	    	block[i]=new BASIS [num_basis_block[i]];
	    }
	    {int ii=0;//initialization of block.
	    for (int i=0;i<num_block;i++){
	    	for (int j=0;j<num_basis_block[i];j++){
	    		block[i][j]=basis_ordered[n][ii];
				//cout << block[i][j].quant_num_totalnum << "  xxxx" << endl;
	    		ii++;
	    	}
	    }}
	    H_bij=new double ** [num_block];
	    for (int b=0;b<num_block;b++){
	        H_bij[b]=new double * [num_basis_block[b]];
	        for (int i=0;i<num_basis_block[b];i++){
	            H_bij[b][i]=new double [num_basis_block[b]];
				for (int j=0;j<num_basis_block[b];j++){
					H_bij[b][i][j]=0;
				}
	        }
	    }
		//cout << "  Dot "<< n;cout << "; Number of blocks: " << setw(4) << num_block << ";    Total number of basis: " << num_basis[n] << endl;
		if (n==1){
			int sum=0;
		    int kk=0;
	        for (int b=0;b<num_block;b++){
				//cout << "    ";cout << "Basis number in block "<< setw(4) << left << b << " : " << setw(4) << num_basis_block[b] << ";    Total electon number: " << setw(4) << block[b][0].quant_num_totalnum << endl;
	            for (int i=0;i<num_basis_block[b];i++){
					//cout << "| " << block[b][i].k << " , " << block[b][i].j << " ; " << n-1 << " > " << block[b][i].quant_num_totalnum << "  " << eigen[n-1][block[b][i].k-1].quant_num_totalnum << "  " << quant_num_totalnum_dot[block[b][i].j-1] << "  " << eigen[n-1][block[b][i].k-1].eig_val_relat << endl;
			        //cout<<"i= "<< i <<"  "<< block[b][i].quant_num_totalnum << "  " << block[b][i].k << "  " << block[b][i].j << "  " << block[b][i].n << endl;
	                for (int j=0;j<num_basis_block[b];j++){
	        		    double sum_up=0;double sum_down=0;
	        			sum_up=sum_up+c_up_dot[block[b][i].j-1][block[b][j].j-1]*c_dag_up_eigen_0[block[b][i].k-1][block[b][j].k-1]*pow(-1,eigen[n-1][block[b][j].k-1].quant_num_totalnum);
	        			sum_up=sum_up+c_dag_up_dot[block[b][i].j-1][block[b][j].j-1]*c_up_eigen_0[block[b][i].k-1][block[b][j].k-1]*pow(-1,eigen[n-1][block[b][i].k-1].quant_num_totalnum);
	        			sum_down=sum_down+c_down_dot[block[b][i].j-1][block[b][j].j-1]*c_dag_down_eigen_0[block[b][i].k-1][block[b][j].k-1]*pow(-1,eigen[n-1][block[b][j].k-1].quant_num_totalnum);
	        			sum_down=sum_down+c_dag_down_dot[block[b][i].j-1][block[b][j].j-1]*c_down_eigen_0[block[b][i].k-1][block[b][j].k-1]*pow(-1,eigen[n-1][block[b][i].k-1].quant_num_totalnum);
	        			H_bij[b][i][j]=(sqrt(Lambda)*eigen[n-1][block[b][i].k-1].eig_val_relat + pow(Lambda,(n-1)/2.0)*(pe_up[n]*quant_num_upnum_dot[block[b][i].j-1] + pe_down[n]*quant_num_downnum_dot[block[b][i].j-1]))*func_delta(block[b][i].k,block[b][j].k)*func_delta(block[b][i].j,block[b][j].j) + pow(Lambda,(n-1)/2.0)*(ptn_up[0]*sum_up+ptn_down[0]*sum_down); //  H_bij=<block[b][i]|H|block[b][j]>. Attention! eigen[n+1][block[b][i].k-1].eigen_value?截断后重新连续排序!
						if (fabs(H_bij[b][i][j]) < 1e-30){
							H_bij[b][i][j]=0;
						}
						//cout << setw(8) << H_bij[b][i][j] << " | " << block[b][i].k << "  " << block[b][i].j << "  " << block[b][i].n << " | " << block[b][j].k << "  " << block[b][j].j << "  " << block[b][j].n << "  | " << block[b][0].quant_num_totalnum << endl;
	        		}
					//cout << " eigen_value: " << eigen[n-1][block[b][i].k-1].eig_val_relat << endl;
	        	}
		        int lda=num_basis_block[b];
				//double * vect_temp=new double [num_basis_block[b]*num_basis_block[b]];
				MKL_Complex16 * vect_temp=new MKL_Complex16 [num_basis_block[b]*num_basis_block[b]];
	            double * value=new double [num_basis_block[b]];
		        {int k=0;
		        for (int i=0;i<num_basis_block[b];i++){
		        	for (int j=0;j<num_basis_block[b];j++){
		        		vect_temp[k].real = H_bij[b][i][j];
		        		vect_temp[k].imag = 0;
		        		k++;
		        	}
		        }}
				double abstol=1e-10;
				int * isuppz=new int [2*num_basis_block[b]];
				double vl,vu;
				int il,iu;
				int * m;
				m=&num_basis_block[b];
		        int info=LAPACKE_zheevr(LAPACK_ROW_MAJOR,'V','A','U',num_basis_block[b],vect_temp,num_basis_block[b],vl,vu,il,iu,abstol,m,value,vect_temp,num_basis_block[b],isuppz);
				//LAPACKE_dsyevd(LAPACK_ROW_MAJOR,'V','U',num_basis_block[b],vect_temp,lda,value);
		        for (int i=0;i<num_basis_block[b];i++){
		            eigen[n][kk].eig_val=value[i];
					//eigen_ordered[n+1][kk].eigen_value=value[i];
			        //eigen_ordered[n+1][kk].quant_num_totalnum=block[b][0].quant_num_totalnum;
					int k=0;
					//for (int j=i*num_basis_block[b];j<(i+1)*num_basis_block[b];j++){
					for (int j=0;j<num_basis_block[b];j++){
						vect[n][kk][sum+k]=vect_temp[i+j*num_basis_block[b]].real;
						if (fabs(vect_temp[i+j*num_basis_block[b]].imag) > 1e-30){
							cout << vect_temp[i+j*num_basis_block[b]].imag << endl;
						}
						//vect[n][kk][sum+k]=vect_temp[j+i*num_basis_block[b]];
						//vect[n][kk][sum+k]=vect_temp[j];
						k++;
		        	}
		            eigen[n][kk].eigen_vect=vect[n][kk];
		            //eigen_ordered[n+1][kk].eigen_vect=vect[n+2][kk];
		            kk++;
					//cout << "    ";cout << "    ";cout << "eigenvalue: " << setw(9) << value[i] << endl;
		        }
				sum=sum+num_basis_block[b];
				//cout << "    ";cout << "Basis number in block "<< setw(4) << left << b << " : " << setw(4) << num_basis_block[b] << ";    Total electon number: " << setw(4) << block[b][0].quant_num_totalnum << scientific << "; min_eigenvalue: " << setw(15) << setprecision(10) << value[0] << "; max_eigenvalue: " << setw(15) << value[num_basis_block[b]-1] << endl;
		        delete [] vect_temp;
				delete [] value;
			}
		}else{
			double ** c_up_dot_basis=new double * [num_basis[n-1]];
			double ** c_down_dot_basis=new double * [num_basis[n-1]];
			double ** c_dag_up_dot_basis=new double * [num_basis[n-1]];
			double ** c_dag_down_dot_basis=new double * [num_basis[n-1]];
			double ** c_up_dot_eigen=new double * [num_eigen_kept[n-1]];
			double ** c_down_dot_eigen=new double * [num_eigen_kept[n-1]];
			double ** c_dag_up_dot_eigen=new double * [num_eigen_kept[n-1]];
			double ** c_dag_down_dot_eigen=new double * [num_eigen_kept[n-1]];
			for (int i=0;i<num_basis[n-1];i++){
				c_up_dot_basis[i]=new double  [num_basis[n-1]];
				c_down_dot_basis[i]=new double  [num_basis[n-1]];
				c_dag_up_dot_basis[i]=new double  [num_basis[n-1]];
				c_dag_down_dot_basis[i]=new double  [num_basis[n-1]];
				for (int j=0;j<num_basis[n-1];j++){
					c_up_dot_basis[i][j]=0;
					c_down_dot_basis[i][j]=0;
					c_dag_up_dot_basis[i][j]=0;
					c_dag_down_dot_basis[i][j]=0;
				}
			}
			for (int i=0;i<num_eigen_kept[n-1];i++){
				c_up_dot_eigen[i]=new double  [num_eigen_kept[n-1]];
				c_down_dot_eigen[i]=new double  [num_eigen_kept[n-1]];
				c_dag_up_dot_eigen[i]=new double  [num_eigen_kept[n-1]];
				c_dag_down_dot_eigen[i]=new double  [num_eigen_kept[n-1]];
				for (int j=0;j<num_eigen_kept[n-1];j++){
					c_up_dot_eigen[i][j]=0;
					c_down_dot_eigen[i][j]=0;
					c_dag_up_dot_eigen[i][j]=0;
					c_dag_down_dot_eigen[i][j]=0;
				}
			}
#pragma omp parallel for 
			for (int i=0;i<num_basis[n-1];i++){
				for (int j=0;j<num_basis[n-1];j++){
					c_up_dot_basis[i][j]=func_delta(basis_ordered[n-1][i].k,basis_ordered[n-1][j].k)*c_up_dot[basis_ordered[n-1][i].j-1][basis_ordered[n-1][j].j-1]*pow(-1,basis_ordered[n-1][j].quant_num_totalnum-quant_num_totalnum_dot[basis_ordered[n-1][j].j-1]);//eigen[n-2][basis_ordered[n-1][j].k-1].quant_num_totalnum);
					c_down_dot_basis[i][j]=func_delta(basis_ordered[n-1][i].k,basis_ordered[n-1][j].k)*c_down_dot[basis_ordered[n-1][i].j-1][basis_ordered[n-1][j].j-1]*pow(-1,basis_ordered[n-1][j].quant_num_totalnum-quant_num_totalnum_dot[basis_ordered[n-1][j].j-1]);//eigen[n-2][basis_ordered[n-1][j].k-1].quant_num_totalnum);
					c_dag_up_dot_basis[i][j]=func_delta(basis_ordered[n-1][i].k,basis_ordered[n-1][j].k)*c_dag_up_dot[basis_ordered[n-1][i].j-1][basis_ordered[n-1][j].j-1]*pow(-1,basis_ordered[n-1][j].quant_num_totalnum-quant_num_totalnum_dot[basis_ordered[n-1][j].j-1]);//eigen[n-2][basis_ordered[n-1][j].k-1].quant_num_totalnum);
					c_dag_down_dot_basis[i][j]=func_delta(basis_ordered[n-1][i].k,basis_ordered[n-1][j].k)*c_dag_down_dot[basis_ordered[n-1][i].j-1][basis_ordered[n-1][j].j-1]*pow(-1,basis_ordered[n-1][j].quant_num_totalnum-quant_num_totalnum_dot[basis_ordered[n-1][j].j-1]);//eigen[n-2][basis_ordered[n-1][j].k-1].quant_num_totalnum);
				}
			}
	        double * matrix_U=new double [num_basis[n-1]*num_eigen_kept[n-1]];
	        double * matrix_c=new double [num_basis[n-1]*num_basis[n-1]];
	        double * matrix_cc=new double [num_eigen_kept[n-1]*num_basis[n-1]];
	        double * matrix_ccc=new double [num_eigen_kept[n-1]*num_eigen_kept[n-1]];
	        {int k=0;
	        for (int i=0;i<num_eigen_kept[n-1];i++){
	        	for (int j=0;j<num_basis[n-1];j++){
	        		matrix_U[k]=eigen[n-1][i].eigen_vect[j];
	        		k++;
	        	}
	        }}
	        {int k=0;
	        for (int i=0;i<num_basis[n-1];i++){
	        	for (int j=0;j<num_basis[n-1];j++){
	        		matrix_c[k]=c_up_dot_basis[j][i];
	        		k++;
	        	}
	        }}
	        cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,num_eigen_kept[n-1],num_basis[n-1],num_basis[n-1],1,matrix_U,num_basis[n-1],matrix_c,num_basis[n-1],0,matrix_cc,num_eigen_kept[n-1]);
	        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,num_eigen_kept[n-1],num_eigen_kept[n-1],num_basis[n-1],1,matrix_cc,num_eigen_kept[n-1],matrix_U,num_basis[n-1],0,matrix_ccc,num_eigen_kept[n-1]);
	        for (int i=0;i<num_eigen_kept[n-1];i++){
	        	for (int j=0;j<num_eigen_kept[n-1];j++){
	        		c_up_dot_eigen[j][i]=matrix_ccc[num_eigen_kept[n-1]*i+j];
	        	}
	        }
	        {int k=0;
	        for (int i=0;i<num_basis[n-1];i++){
	        	for (int j=0;j<num_basis[n-1];j++){
	        		matrix_c[k]=c_down_dot_basis[j][i];
	        		k++;
	        	}
	        }}
	        cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,num_eigen_kept[n-1],num_basis[n-1],num_basis[n-1],1,matrix_U,num_basis[n-1],matrix_c,num_basis[n-1],0,matrix_cc,num_eigen_kept[n-1]);
	        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,num_eigen_kept[n-1],num_eigen_kept[n-1],num_basis[n-1],1,matrix_cc,num_eigen_kept[n-1],matrix_U,num_basis[n-1],0,matrix_ccc,num_eigen_kept[n-1]);
	        for (int i=0;i<num_eigen_kept[n-1];i++){
	        	for (int j=0;j<num_eigen_kept[n-1];j++){
	        		c_down_dot_eigen[j][i]=matrix_ccc[num_eigen_kept[n-1]*i+j];
	        	}
	        }
	        {int k=0;
	        for (int i=0;i<num_basis[n-1];i++){
	        	for (int j=0;j<num_basis[n-1];j++){
	        		matrix_c[k]=c_dag_up_dot_basis[j][i];
	        		k++;
	        	}
	        }}
	        cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,num_eigen_kept[n-1],num_basis[n-1],num_basis[n-1],1,matrix_U,num_basis[n-1],matrix_c,num_basis[n-1],0,matrix_cc,num_eigen_kept[n-1]);
	        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,num_eigen_kept[n-1],num_eigen_kept[n-1],num_basis[n-1],1,matrix_cc,num_eigen_kept[n-1],matrix_U,num_basis[n-1],0,matrix_ccc,num_eigen_kept[n-1]);
	        for (int i=0;i<num_eigen_kept[n-1];i++){
	        	for (int j=0;j<num_eigen_kept[n-1];j++){
	        		c_dag_up_dot_eigen[j][i]=matrix_ccc[num_eigen_kept[n-1]*i+j];
	        	}
	        }
	        {int k=0;
	        for (int i=0;i<num_basis[n-1];i++){
	        	for (int j=0;j<num_basis[n-1];j++){
	        		matrix_c[k]=c_dag_down_dot_basis[j][i];
	        		k++;
	        	}
	        }}
	        cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,num_eigen_kept[n-1],num_basis[n-1],num_basis[n-1],1,matrix_U,num_basis[n-1],matrix_c,num_basis[n-1],0,matrix_cc,num_eigen_kept[n-1]);
	        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,num_eigen_kept[n-1],num_eigen_kept[n-1],num_basis[n-1],1,matrix_cc,num_eigen_kept[n-1],matrix_U,num_basis[n-1],0,matrix_ccc,num_eigen_kept[n-1]);
	        for (int i=0;i<num_eigen_kept[n-1];i++){
	        	for (int j=0;j<num_eigen_kept[n-1];j++){
	        		c_dag_down_dot_eigen[j][i]=matrix_ccc[num_eigen_kept[n-1]*i+j];
	        	}
	        }
	        delete [] matrix_U;
	        delete [] matrix_c;
	        delete [] matrix_cc;
	        delete [] matrix_ccc;
			/*
			 *for (int i=0;i<num_basis[n-1];i++){
			 *    for (int j=0;j<num_basis[n-1];j++){
			 *        cout << setw(10) << c_up_dot_basis[i][j] << " basis" << " |  " << setw(3) << basis_ordered[n-1][i].k-1 << "  " << setw(3) << basis_ordered[n-1][i].j-1 << " |  " << setw(3) << basis_ordered[n-1][j].k-1 << "  " << setw(3) << basis_ordered[n-1][j].j-1 << endl;
			 *        //cout << setw(10) << c_up_dot_basis[i][j] << setw(10) << c_dag_up_dot_basis[j][i] << "  up  " << i << "  " << j << "  " << n << endl;
			 *        //cout << setw(10) << c_up_dot_basis[i][j] << "  up  " << i << "  " << j << " |  " << basis_ordered[n-1][i].k-1 << "  " << basis_ordered[n-1][i].j-1 << " |  " << basis_ordered[n-1][j].k-1 << "  " << basis_ordered[n-1][j].j-1 << endl;
			 *    }
			 *}
			 */
			int sum=0;
			int kk=0;
			BASIS * basis_old=new BASIS [num_basis[n-1]];//n=0
			for (int i=0;i<num_basis[n-1];i++){
				basis_old[i]=basis_ordered[n-1][i];
			}
		    for (int b=0;b<num_block;b++){
#pragma omp parallel for 
		        for (int i=0;i<num_basis_block[b];i++){
					//cout << "| " << block[b][i].k << " , " << block[b][i].j << " ; " << n-1 << " > " << block[b][i].quant_num_totalnum << "  " << eigen[n-1][block[b][i].k-1].quant_num_totalnum << "  " << quant_num_totalnum_dot[block[b][i].j-1] << "  " << eigen[n-1][block[b][i].k-1].eig_val_relat << endl;
		            for (int j=0;j<num_basis_block[b];j++){
						double sum_up2=0;double sum_down2=0;double sum_up=0;double sum_down=0;
						sum_up=pow(-1,eigen[n-1][block[b][j].k-1].quant_num_totalnum)*c_dag_up_dot_eigen[block[b][i].k-1][block[b][j].k-1]*c_up_dot[block[b][i].j-1][block[b][j].j-1] + pow(-1,eigen[n-1][block[b][i].k-1].quant_num_totalnum)*c_dag_up_dot[block[b][i].j-1][block[b][j].j-1]*c_up_dot_eigen[block[b][i].k-1][block[b][j].k-1];
						sum_down=pow(-1,eigen[n-1][block[b][j].k-1].quant_num_totalnum)*c_dag_down_dot_eigen[block[b][i].k-1][block[b][j].k-1]*c_down_dot[block[b][i].j-1][block[b][j].j-1] + pow(-1,eigen[n-1][block[b][i].k-1].quant_num_totalnum)*c_dag_down_dot[block[b][i].j-1][block[b][j].j-1]*c_down_dot_eigen[block[b][i].k-1][block[b][j].k-1];
						//H_bij[b][i][j]=(sqrt(Lambda)*eigen[n-1][block[b][i].k-1].eig_val_relat )*func_delta(block[b][i].k,block[b][j].k)*func_delta(block[b][i].j,block[b][j].j) + pow(Lambda,(n-2)/2.0)*(ptn_up[n-2]*sum_up+ptn_down[n-2]*sum_down);   //  H_bij=<block[b][i]|H|block[b][j]>. Attention! eigen[n+1][block[b][i].k-1].eigen_value?截断后重新连续排序!
						H_bij[b][i][j]=(sqrt(Lambda)*eigen[n-1][block[b][i].k-1].eig_val_relat + pow(Lambda,(n-1)/2.0)*(pe_up[n]*quant_num_upnum_dot[block[b][i].j-1] +pe_down[n]*quant_num_downnum_dot[block[b][i].j-1]))*func_delta(block[b][i].k,block[b][j].k)*func_delta(block[b][i].j,block[b][j].j) + pow(Lambda,(n-1)/2.0)*(ptn_up[n-1]*sum_up+ptn_down[n-1]*sum_down);   //  H_bij=<block[b][i]|H|block[b][j]>. Attention! eigen[n+1][block[b][i].k-1].eigen_value?截断后重新连续排序!
						if (fabs(H_bij[b][i][j]) < 1e-30){
							H_bij[b][i][j]=0;
						}
						//cout << setw(12) << setprecision(9) << H_bij[b][i][j] << " | " << setw(4) << block[b][i].k << setw(4) << block[b][i].j << setw(4) << block[b][i].n << " | " << setw(4) << block[b][j].k << setw(4) << block[b][j].j << setw(4) << block[b][j].n << " | " << setw(4) << block[b][i].quant_num_totalnum << setw(4) << block[b][i].quant_num_totalnum << " | c_N_" << block[b][i].k-1 << "_" << block[b][j].k-1 << setw(15) << right << c_dag_up_dot_eigen[block[b][i].k-1][block[b][j].k-1] << setw(3) << i+1 << setw(3) << j+1 << left << endl;
		    		}
					//cout << " eigen_value: " << eigen[n-1][block[b][i].k-1].eig_val << endl;
		    	}
		        int lda=num_basis_block[b];
				//double * vect_temp=new double [num_basis_block[b]*num_basis_block[b]];
				MKL_Complex16 * vect_temp=new MKL_Complex16 [num_basis_block[b]*num_basis_block[b]];
	            double * value=new double [num_basis_block[b]];
		        {int k=0;
		        for (int i=0;i<num_basis_block[b];i++){
		        	for (int j=0;j<num_basis_block[b];j++){
		        		vect_temp[k].real = H_bij[b][i][j];
		        		vect_temp[k].imag = 0;
		        		k++;
		        	}
		        }}
				double abstol=1e-10;
				int * isuppz=new int [2*num_basis_block[b]];
				double vl,vu;
				int il,iu;
				int * m;
				m=&num_basis_block[b];
		        int info=LAPACKE_zheevr(LAPACK_ROW_MAJOR,'V','A','U',num_basis_block[b],vect_temp,num_basis_block[b],vl,vu,il,iu,abstol,m,value,vect_temp,num_basis_block[b],isuppz);
				//LAPACKE_dsyevd(LAPACK_ROW_MAJOR,'V','U',num_basis_block[b],vect_temp,lda,value);
		        for (int i=0;i<num_basis_block[b];i++){
		            eigen[n][kk].eig_val=value[i];
					//cout << "eigenvalue: " << setw(3) << i << setw(9) << value[i] << endl;
					//cout << "eigval  " << setw(3) << i << setw(9) << eigen[n][kk].eig_val << endl;
					int k=0;
					for (int j=0;j<num_basis_block[b];j++){//row-wise
						vect[n][kk][sum+k]=vect_temp[i+j*num_basis_block[b]].real;
						if (fabs(vect_temp[i+j*num_basis_block[b]].imag) > 1e-30){
							cout << vect_temp[i+j*num_basis_block[b]].imag << endl;
						}
						//vect[n][kk][sum+k]=vect_temp[j];
						k++;
		        	}
		            eigen[n][kk].eigen_vect=vect[n][kk];
		            kk++;
		        }
				sum=sum+num_basis_block[b];
				//cout << "    ";cout << "Basis number in block "<< setw(4) << left << b << " : " << setw(4) << num_basis_block[b] << ";    Total electon number: " << setw(4) << block[b][0].quant_num_totalnum << scientific << "; min_eigenvalue: " << setw(15) << setprecision(10) << value[0] << "; max_eigenvalue: " << setw(15) << value[num_basis_block[b]-1] << endl;
			    delete [] vect_temp;
				delete [] value;
		    }
			for (int i=0;i<num_eigen_kept[n-1];i++){
				delete [] c_dag_up_dot_eigen[i];
				delete [] c_dag_down_dot_eigen[i];
				delete [] c_up_dot_eigen[i];
				delete [] c_down_dot_eigen[i];
			}
			for (int i=0;i<num_basis[n-1];i++){
				delete [] c_dag_up_dot_basis[i];
				delete [] c_dag_down_dot_basis[i];
				delete [] c_up_dot_basis[i];
				delete [] c_down_dot_basis[i];
			}
			delete [] c_dag_up_dot_eigen;
			delete [] c_dag_down_dot_eigen;
			delete [] c_up_dot_eigen;
			delete [] c_down_dot_eigen;
			delete [] basis_old;
			delete [] c_dag_up_dot_basis;
			delete [] c_dag_down_dot_basis;
			delete [] c_up_dot_basis;
			delete [] c_down_dot_basis;
		}
		for (int k=0;k<num_basis[n];k++){//initionialization of sort,k,n.
			eigen[n][k].sort=k+1;
			eigen[n][k].k=k+1;
			eigen[n][k].n=n;
		}
		{int sum1=0,sum2=0;
		for(int b=0;b<1;b++){//initionialization of quantum number.
			sum2=sum2+num_basis_block[b];
			for (int i=0;i<num_basis_block[b];i++){
			    eigen[n][i].quant_num_totalnum=block[b][0].quant_num_totalnum;
			    eigen[n][i].quant_num_upnum=block[b][0].quant_num_upnum;
			    eigen[n][i].quant_num_downnum=block[b][0].quant_num_downnum;
				//cout << "basi  " << block[b][0].quant_num_totalnum << "  " << block[b][0].quant_num_upnum << "  " << block[b][0].quant_num_downnum << endl;
				//cout << "basi  " << eigen[n][i].quant_num_totalnum << "  " << eigen[n][i].quant_num_upnum << "  " << eigen[n][i].quant_num_downnum << endl;
				//cout << " before sort " << i+1 << "   " << eigen[n][i].quant_num_totalnum << "  " << eigen[n][i].quant_num_upnum-eigen[n][i].quant_num_downnum  << "  " << setprecision(12) << eigen[n][i].eig_val << endl;
			}
		}
		for(int b=1;b<num_block;b++){//initionialization of quantum number.
			sum2=sum2+num_basis_block[b];
			sum1=sum1+num_basis_block[b-1];
			for (int i=sum1;i<sum2;i++){
			    eigen[n][i].quant_num_totalnum=block[b][0].quant_num_totalnum;
			    eigen[n][i].quant_num_upnum=block[b][0].quant_num_upnum;
			    eigen[n][i].quant_num_downnum=block[b][0].quant_num_downnum;
				//cout << "basis  " << block[b][0].quant_num_totalnum << "  " << block[b][0].quant_num_upnum << "  " << block[b][0].quant_num_downnum << endl;
				//cout << "basis  " << eigen[n][i].quant_num_totalnum << "  " << eigen[n][i].quant_num_upnum << "  " << eigen[n][i].quant_num_downnum << endl;
				//cout << " before sort " << i+1 << "   " << eigen[n][i].quant_num_totalnum << "  " << eigen[n][i].quant_num_upnum-eigen[n][i].quant_num_downnum  << "  " << setprecision(12) << eigen[n][i].eig_val << endl;
			}
		}}
	    for (int k=1;k<=num_basis[n]-1;k++){//sort of eigen_value.
	    	for (int ii=1;ii<=num_basis[n]-k;ii++){
	    		if (eigen[n][ii-1].eig_val > eigen[n][ii].eig_val){
	    			EIGEN_STATE temp = eigen[n][ii-1];//.quant_num_totalnum;
	    			eigen[n][ii-1]=eigen[n][ii];
	    			eigen[n][ii]=temp;
	    		}
	    	}
			eigen[n][num_basis[n]-k].k=num_basis[n]-k+1;
	    }
		//cout << eigen[n][0].k << endl;
		eigen[n][0].k=1;
		for (int i=0;i<num_basis[n];i++){
		    eigen[n][i].eig_val_relat=eigen[n][i].eig_val-eigen[n][0].eig_val;
			//cout << "in eigen  " << eigen[n][i].k-1 << endl;
			//cout << n << "  " << basis_ordered[n][i].quant_num_totalnum <<" | "<<  basis_ordered[n][i].k <<"  "<< basis_ordered[n][i].j << " | " << eigen[n][i].sort << endl;
		}
		//ofstream cup("cup",ios::binary);
		E_GS[n]=sqrt(Lambda)*E_GS[n-1]+eigen[n][0].eig_val;
		f_E_GS << "Dot  " << n << scientific << setw(25) << setprecision(15) << E_GS[n] << setw(25) <<  E_GS[n]*pow(Lambda,-1.0*(n-1)/2.0) << endl;
		for (int i=0;i<num_basis[n];i++){
			if (Q){
			    f_eig_val << "Dot  " << n << "  Q_" << left << setw(5) << scientific << eigen[n][i].quant_num_totalnum << setw(25) << setprecision(15) << eigen[n][i].eig_val << setw(25) << setprecision(15) << eigen[n][i].eig_val_relat << setw(25) << (eigen[n][i].eig_val_relat+E_GS[n])*pow(Lambda,-1.0*(n-1)/2.0) << "   " << eigen[n][i].k << endl;
			}else if(Q_Sz){
			    f_eig_val << "Dot  " << n << "  Q_" << left << setw(5) << scientific << eigen[n][i].quant_num_totalnum << "  Sz_" << left << setw(5) << scientific << (eigen[n][i].quant_num_upnum-eigen[n][i].quant_num_downnum)/2.0  << setw(25) << setprecision(15) << eigen[n][i].eig_val << setw(25) << setprecision(15) << eigen[n][i].eig_val_relat << setw(25) << (eigen[n][i].eig_val_relat+E_GS[n])*pow(Lambda,-1.0*(n-1)/2.0) << "   " << eigen[n][i].k << endl;
			}else if(N_up_N_down){
			    f_eig_val << "Dot  " << n << "  N_up_" << left << setw(5) << scientific << eigen[n][i].quant_num_upnum << "  N_down_" << left << setw(5) << scientific << eigen[n][i].quant_num_downnum << setw(25) << setprecision(15) << eigen[n][i].eig_val << setw(25) << setprecision(15) << eigen[n][i].eig_val_relat << (eigen[n][i].eig_val_relat+E_GS[n])*pow(Lambda,-1.0*(n-1)/2.0) << "   " << eigen[n][i].k << endl;
			}
			/*
			 *for (int j=0;j<num_basis[n];j++){
			 *    f_U << i << "    " << j << "    " << eigen[n][j].eigen_vect[i] << endl;
			 *    f_d_up << i << "    " << j << "    " << c_up_eigen[n][i][j] << endl;
			 *    f_d_down << i << "    " << j << "    " << c_down_eigen[n][i][j] << endl;
			 *}
			 */
		}
		//cout << "    ";cout << "Time leaved:    ";date_time();cout << endl;
		for (int i=0;i<num_block;i++){
			delete [] block[i];
		}
		delete [] block;
		for (int b=0;b<num_block;b++){
		    for (int i=0;i<num_basis_block[b];i++){
		        delete [] H_bij[b][i];
		    }
		    delete [] H_bij[b];
		}
		delete [] H_bij;
		delete [] num_basis_block;
		//if (n==6) exit(0);
	}
	cout << "Time leaved:    ";date_time();cout << endl;
}
