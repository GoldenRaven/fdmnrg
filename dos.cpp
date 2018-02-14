#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include<math.h>
#include<mkl.h>
#include<omp.h>
#include"setup.h"
using namespace std;
double Pi=3.141592653589793238;
double *** rho_red_temp;
double (*pf)(double,double);
double P_K(double,double);
double **** eigen_sigma;
double ** eigen_value;
double **** rho_red;
void density_of_state(void)
{
    cout << "density_of_state: "; date_time();cout << endl;
    double dos1_up(double);
    double dos1_down(double);
    double dos2_up(double);
    double dos2_down(double);
    double dos3_up(double);
    double dos3_down(double);
    void   rho_red_operator(int);
    double P_LG(double,double);
    double P_G(double,double);
    double P_lorentz(double,double);
    double P_sum_rule(double,double);
    eigen_sigma=new double *** [N_max+2];
    for (int n=0;n<1;n++){//eigen_sigma[0][sigma][k][j] shouldn't be used!
        eigen_sigma[n]=new double ** [1];
        for (int sigma=0;sigma<1;sigma++){
            eigen_sigma[n][sigma]=new double * [1];
            for (int k=0;k<1;k++){
                eigen_sigma[n][sigma][k]=new double [1];
                for (int j=0;j<1;j++){
                    eigen_sigma[n][sigma][k][j]=0;
                }
            }
        }
    }
    for (int n=1;n<N_max+2;n++){
        eigen_sigma[n]=new double ** [dim_dot];
        for (int sigma=0;sigma<dim_dot;sigma++){
            eigen_sigma[n][sigma]=new double * [num_eigen_kept[n-1]];
            for (int k=0;k<num_eigen_kept[n-1];k++){
                eigen_sigma[n][sigma][k]=new double [num_basis[n]];
            }
        }
    }
    for (int n=1;n<N_max+2;n++){//initialization of eigen_sigma[n][sigma][k][j].
        for (int sigma=0;sigma<dim_dot;sigma++){
            for (int k=0;k<num_eigen_kept[n-1];k++){
                for (int j=0;j<num_basis[n];j++){
                    //transformation matrix after sorting.
                    eigen_sigma[n][sigma][k][j]=eigen[n][j].eigen_vect[basis_kj[n][k][sigma].sort];
                    //eigen_sigma[n][sigma][k][j]=vect[n][j][basis_kj[n][k][sigma].sort];
                }
            }
        }
    }
    eigen_value = new double * [N_max+2];
    for (int n=0;n<N_max+2;n++){
        eigen_value[n] = new double [num_basis[n]];
        for (int i=0;i<num_basis[n];i++){
            eigen_value[n][i] = eigen[n][i].eig_val;
        }
    }
    for (int n=0;n<N_max+2;n++){
        delete [] eigen[n];
    }
    delete [] eigen;
    rho_red      = new double **  [N_max+1];//reduced density matrix.
    rho_red_temp = new double *** [N_max+1];//summing term in each rho_red.
    for (int n=0;n<N_max+1;n++){
        rho_red[n]=new double * [num_kept];
        for (int i=0;i<num_kept;i++){
            rho_red[n][i]=new double [num_kept];
            for (int j=0;j<num_kept;j++){
                rho_red[n][i][j]=0;
            }
        }
    }
    for (int n=0;n<n0;n++){//n < n0 not used!
        rho_red_temp[n] = new double ** [1];
        for(int n2=0;n2<1;n2++){
            rho_red_temp[n][n2] = new double * [1];
            for (int i=0;i<1;i++){
                rho_red_temp[n][n2][i] = new double [1];
            }
        }
    }
    for (int n1=N_max;n1>=n0;n1--){
        rho_red_operator(n1);
    }
    for (int n=0;n<n0;n++){
        for(int n2=0;n2<1;n2++){
            for (int i=0;i<1;i++){
                delete [] rho[n][n2][i];
            }
            delete [] rho[n][n2];
        }
        delete [] rho[n];
    }
    for (int n1=N_max-2;n1>=n0;n1--){
        for (int n2=0;n2<N_max-2-n1+1;n2++){
            for (int i=0;i<num_kept;i++){
                delete [] rho[n1][n2][i];
            }
            delete [] rho[n1][n2];
        }
        delete [] rho[n1];
    }
    delete [] rho;
    cout << "  Time finished all reduced density matrice:    ";date_time();cout << endl;
    cout << scientific << left << "            ";
    cout << setw(14) << "freqency"<< setw(14) << "dos1_up"<< setw(14) << "dos2_up";
    cout << setw(14) << "dos3_up";
    cout << setw(14) << "dos1_down"<< setw(14) << "dos2_down";
    cout << setw(14) << "dos3_down";
    cout << endl;
    double DOS1_UP,DOS2_UP,DOS3_UP,DOS1_DOWN,DOS2_DOWN,DOS3_DOWN;
    //test sum rule
    double freqency;
    freqency=1.0;
    pf=P_sum_rule;
    DOS1_UP=dos1_up(freqency);
    DOS2_UP=dos2_up(freqency);
    DOS3_UP=dos3_up(freqency);
    DOS1_DOWN=dos1_down(freqency);
    DOS2_DOWN=dos2_down(freqency);
    DOS3_DOWN=dos3_down(freqency);
    cout << "  sum_rule  ";
    cout << setw(14) << setprecision(5) << freqency;
    cout << setw(14) << setprecision(5) << DOS1_UP;
    cout << setw(14) << setprecision(5) << DOS2_UP;
    cout << setw(14) << setprecision(5) << DOS3_UP;
    cout << setw(14) << setprecision(5) << DOS1_DOWN;
    cout << setw(14) << setprecision(5) << DOS2_DOWN;
    cout << setw(14) << setprecision(5) << DOS3_DOWN;
    cout << "  sum_up= ";
    cout << setw(20) << setprecision(15) << DOS1_UP+DOS2_UP+DOS3_UP;
    cout << "  sum_down= ";
    cout << setw(20) << setprecision(15) << DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
    cout << endl;
    pf=P_K;
    if (smear) {
        ifstream f_freqency("freqency");
        ofstream f_dos_smeared("freq_dos_smeared.dat");
        while (!f_freqency.eof()){//???????????
            f_freqency >> freqency;
            //if (fabs(freqency) < omega0){
            DOS1_UP=dos1_up(freqency);
            DOS2_UP=dos2_up(freqency);
            DOS3_UP=dos3_up(freqency);
            DOS1_DOWN=dos1_down(freqency);
            DOS2_DOWN=dos2_down(freqency);
            DOS3_DOWN=dos3_down(freqency);
            cout << "            ";
            cout << setw(14) << setprecision(5) << freqency;
            cout << setw(14) << setprecision(5) << DOS1_UP;
            cout << setw(14) << setprecision(5) << DOS2_UP;
            cout << setw(14) << setprecision(5) << DOS3_UP;
            cout << setw(14) << setprecision(5) << DOS1_DOWN;
            cout << setw(14) << setprecision(5) << DOS2_DOWN;
            cout << setw(14) << setprecision(5) << DOS3_DOWN;
            date_time();cout << endl;
            //}
            f_dos_smeared << left;
            f_dos_smeared << setw(20) << setprecision(10) << freqency;
            f_dos_smeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP;
            f_dos_smeared << setw(20) << setprecision(10) << DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
            f_dos_smeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP + DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
            f_dos_smeared << endl;
        }
        f_freqency.close();
        f_dos_smeared.close();
    }
    if (unsmear){
        bool smear_tmp=smear;
        smear=false;
        ifstream f_freqency("freqency");
        ofstream f_dos_unsmeared("freq_dos_unsmeared.dat");
        while (!f_freqency.eof()){//???????????
            f_freqency >> freqency;
            DOS1_UP=dos1_up(freqency);
            DOS2_UP=dos2_up(freqency);
            DOS3_UP=dos3_up(freqency);
            DOS1_DOWN=dos1_down(freqency);
            DOS2_DOWN=dos2_down(freqency);
            DOS3_DOWN=dos3_down(freqency);
            cout << "           ";
            //cout << "  unsmeared ";
            cout << setw(14) << setprecision(5) << freqency;
            cout << setw(14) << setprecision(5) << DOS1_UP;
            cout << setw(14) << setprecision(5) << DOS2_UP;
            cout << setw(14) << setprecision(5) << DOS3_UP;
            cout << setw(14) << setprecision(5) << DOS1_DOWN;
            cout << setw(14) << setprecision(5) << DOS2_DOWN;
            cout << setw(14) << setprecision(5) << DOS3_DOWN;
            date_time();cout << endl;
            f_dos_unsmeared << left;
            f_dos_unsmeared << setw(20) << setprecision(10) << freqency;
            f_dos_unsmeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP;
            f_dos_unsmeared << setw(20) << setprecision(10) << DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
            f_dos_unsmeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP + DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
            f_dos_unsmeared << endl;
        }
        smear=smear_tmp;
        //exit(0);
        f_freqency.close();
        f_dos_unsmeared.close();
    }
    cout << "Time leaved:    ";date_time();cout << endl;
}

void rho_red_operator(int n1)
{
    rho_red_temp[n1] = new double ** [N_max+1-n1];
    for (int n2=0;n2<N_max+1-n1;n2++){
        rho_red_temp[n1][n2] = new double * [num_kept];
        for (int i=0;i<num_kept;i++){
            rho_red_temp[n1][n2][i] = new double [num_kept];
            for (int j=0;j<num_kept;j++){
                rho_red_temp[n1][n2][i][j] = 0;
            }
        }
    }
    double * matrix_U=new double [num_basis[n]*num_basis[n]];
    double * matrix_c=new double [num_basis[n]*num_basis[n]];
    double * matrix_cc=new double [num_basis[n]*num_basis[n]];
    double * matrix_ccc=new double [num_basis[n]*num_basis[n]];

    for (int n2=0;n2<N_max+1-n1;n2++){
        if (n2 == 0){
            for (int sigma=0;sigma<dim_dot;sigma++){
                {int k=0;//colume-wise
                    for (int i=0;i<num_basis[n];i++){
                        for (int j=0;j<num_basis[n];j++){
                            matrix_U[k]=eigen[n][i].eigen_vect[j];
                            k++;
                        }
                    }}
                {int k=0;
                    for (int i=0;i<num_basis[n];i++){
                        for (int j=0;j<num_basis[n];j++){
                            matrix_c[k]=occu_imp_up_basis[n][j][i];
                            k++;
                        }
                    }}
                cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,num_basis[n],num_basis[n],num_basis[n],1,matrix_U,num_basis[n],matrix_c,num_basis[n],0,matrix_cc,num_basis[n]);
                cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,num_basis[n],num_basis[n],num_basis[n],1,matrix_cc,num_basis[n],matrix_U,num_basis[n],0,matrix_ccc,num_basis[n]);
                for (int i=0;i<num_basis[n];i++){
                    for (int j=0;j<num_basis[n];j++){
                        occu_imp_up_eigen[n][j][i]=matrix_ccc[num_basis[n]*i+j];
                    }
                }

            }
        }
    }
}

double dos1_up(double freqency)
{
    double sum=0;
    for (int n=n0;n<N_max+1;n++){
        double sum_n=0;
#pragma omp parallel for reduction(+:sum_n)
        for (int i=num_kept;i<num_basis[n];i++){
            for (int j=num_kept;j<num_basis[n];j++){
                double omegan=0;
                omegan=pow(Lambda,-1.0*(eigen[n][0].n-1)/2.0)*(eigen_value[n][i]-eigen_value[n][j]);
                sum_n=sum_n+c_up_eigen[n][i][j]*c_up_eigen[n][i][j]*(exp_z[n][i-num_kept]+exp_z[n][j-num_kept])*(*pf)(freqency,-1.0*omegan);
            }
        }
        sum=sum+pow(4,eigen[N_max+1][0].n-eigen[n][0].n)*sum_n;
    }
    for (int n=N_max+1;n<N_max+2;n++){
        double sum_n=0;
#pragma omp parallel for reduction(+:sum_n)
        for (int i=0;i<num_basis[n];i++){
            for (int j=0;j<num_basis[n];j++){
                double omegan=0;
                omegan=pow(Lambda,-1.0*(n-1-1)/2.0)*(eigen_value[n][i]-eigen_value[n][j]);
                sum_n=sum_n+c_up_eigen[n][i][j]*c_dag_up_eigen[n][j][i]*(exp_z[n][i]+exp_z[n][j])*(*pf)(freqency,-1.0*omegan);
                //sum_n=sum_n+c_dag_up_eigen[n][i][j]*c_up_eigen[n][j][i]*(exp(-1.0*Beta*pow(Lambda,-1.0*(n-1-1)/2.0)*eigen[n][i].eig_val)+exp(-1.0*Beta*pow(Lambda,-1.0*(n-1-1)/2.0)*eigen[n][j].eig_val))*P_K(freqency,-1.0*omegan);
            }
        }
        sum=sum+pow(4,N_max-1-n)*sum_n;
    }
    return sum;
}
double dos1_down(double freqency)
{
    double sum=0;
    double P_h_unsmeared(double freqency);
    double P_h_smeared(double freqency);
    for (int n=n0;n<N_max-1;n++){
        double sum_n=0;
#pragma omp parallel for reduction(+:sum_n)
        for (int i=num_kept;i<num_basis[n];i++){
            for (int j=num_kept;j<num_basis[n];j++){
                double omegan=0;
                omegan=pow(Lambda,-1.0*(n-1-1)/2.0)*(eigen_value[n][i]-eigen_value[n][j]);
                sum_n=sum_n+c_down_eigen[n][i][j]*c_dag_down_eigen[n][j][i]*(exp_z[n][i-num_kept]+exp_z[n][j-num_kept])*(*pf)(freqency,-1.0*omegan);
                //sum_n=sum_n+c_dag_down_eigen[n][i][j]*c_down_eigen[n][j][i]*(exp(-1.0*Beta*pow(Lambda,-1.0*(n-1-1)/2.0)*eigen[n][i].eig_val)+exp(-1.0*Beta*pow(Lambda,-1.0*(n-1-1)/2.0)*eigen[n][j].eig_val))*P_K(freqency,-1.0*omegan);
            }
        }
        sum=sum+pow(4,N_max-1-n)*sum_n;
    }
    for (int n=N_max-1;n<N_max;n++){
        double sum_n=0;
#pragma omp parallel for reduction(+:sum_n)
        for (int i=0;i<num_basis[n];i++){
            for (int j=0;j<num_basis[n];j++){
                double omegan=0;
                omegan=pow(Lambda,-1.0*(n-1-1)/2.0)*(eigen_value[n][i]-eigen_value[n][j]);
                sum_n=sum_n+c_down_eigen[n][i][j]*c_dag_down_eigen[n][j][i]*(exp_z[n][i]+exp_z[n][j])*(*pf)(freqency,-1.0*omegan);
            }
        }
        sum=sum+pow(4,N_max-1-n)*sum_n;
    }
    return sum;
}
double dos2_up(double freqency)
{
    double sum=0;
    double P_h_unsmeared(double freqency);
    double P_h_smeared(double freqency);
    for (int n=n0;n<N_max-1;n++){
        double sum_n=0;
#pragma omp parallel for reduction(+:sum_n)
        for (int i=num_kept;i<num_basis[n];i++){
            for (int j=0;j<num_kept;j++){
                double omegan=0;
                omegan=(eigen_value[n][i]-eigen_value[n][j])*pow(Lambda,-1.0*(n-1-1)/2.0);
                sum_n=sum_n+c_up_eigen[n][i][j]*c_dag_up_eigen[n][j][i]*exp_z[n][i-num_kept]*(*pf)(freqency,-1.0*omegan)+c_dag_up_eigen[n][i][j]*c_up_eigen[n][j][i]*exp_z[n][i-num_kept]*(*pf)(freqency,omegan);
                //sum_n=sum_n+c_dag_up_eigen[n][i][j]*c_up_eigen[n][j][i]*exp(-1.0*Beta*eigen[n][i].eig_val*pow(Lambda,-1.0*(n-1-1)/2.0))*P_K(freqency,-1.0*omegan)+c_up_eigen[n][i][j]*c_dag_up_eigen[n][j][i]*exp(-1.0*Beta*eigen[n][i].eig_val*pow(Lambda,-1.0*(n-1-1)/2.0))*P_K(freqency,omegan);
            }
        }
        sum=sum+sum_n*pow(4,N_max-1-n);
    }
    return sum;
}
double dos2_down(double freqency)
{
    double sum=0;
    double P_h_unsmeared(double freqency);
    double P_h_smeared(double freqency);
    for (int n=n0;n<N_max-1;n++){
        double sum_n=0;
#pragma omp parallel for reduction(+:sum_n)
        for (int i=num_kept;i<num_basis[n];i++){
            for (int j=0;j<num_kept;j++){
                double omegan=0;
                omegan=(eigen_value[n][i]-eigen_value[n][j])*pow(Lambda,-1.0*(n-1-1)/2.0);
                sum_n=sum_n+c_down_eigen[n][i][j]*c_dag_down_eigen[n][j][i]*exp_z[n][i-num_kept]*(*pf)(freqency,-1.0*omegan)+c_dag_down_eigen[n][i][j]*c_down_eigen[n][j][i]*exp_z[n][i-num_kept]*(*pf)(freqency,omegan);
                //sum_n=sum_n+c_dag_down_eigen[n][i][j]*c_down_eigen[n][j][i]*exp(-1.0*Beta*eigen[n][i].eig_val*pow(Lambda,-1.0*(n-1-1)/2.0))*P_K(freqency,-1.0*omegan)+c_down_eigen[n][i][j]*c_dag_down_eigen[n][j][i]*exp(-1.0*Beta*eigen[n][i].eig_val*pow(Lambda,-1.0*(n-1-1)/2.0))*P_K(freqency,omegan);
            }
        }
        sum=sum+sum_n*pow(4,N_max-1-n);
    }
    return sum;
}
double dos3_up(double freqency)
{
    double sum=0;
    double P_h_unsmeared(double freqency);
    double P_h_smeared(double freqency);
    //ofstream f_up("up.dat");
    for (int n=n0;n<N_max-1;n++){
        double sum_n=0;
        double * matrix_rho_KK=new double [num_eigen_kept[n]*num_eigen_kept[n]];
        double * matrix_B=new double [(num_basis[n]-num_eigen_kept[n])*num_eigen_kept[n]];
        double * matrix_prod1_DK=new double [(num_basis[n]-num_eigen_kept[n])*num_eigen_kept[n]];
        double * matrix_prod2_KD=new double [num_eigen_kept[n]*(num_basis[n]-num_eigen_kept[n])];
        {int k=0;
            for (int i=0;i<num_eigen_kept[n];i++){
                for (int j=0;j<num_eigen_kept[n];j++){
                    matrix_rho_KK[k]=rho_red_temp[n][j][i];
                    k++;
                }
            }}
        {int k=0;
            for (int i=0;i<num_eigen_kept[n];i++){
                for (int j=num_eigen_kept[n];j<num_basis[n];j++){
                    matrix_B[k]=c_dag_up_eigen[n][j][i];
                    k++;
                }
            }}
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,num_basis[n]-num_eigen_kept[n],num_eigen_kept[n],num_eigen_kept[n],1,matrix_B,num_basis[n]-num_eigen_kept[n],matrix_rho_KK,num_eigen_kept[n],0,matrix_prod1_DK,num_basis[n]-num_eigen_kept[n]);
        {int k=0;
            for (int i=num_kept;i<num_basis[n];i++){
                for (int j=0;j<num_kept;j++){
                    matrix_B[k]=c_dag_up_eigen[n][j][i];
                    k++;
                }
            }}
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,num_kept,num_basis[n]-num_kept,num_kept,1,matrix_rho_KK,num_kept,matrix_B,num_kept,0,matrix_prod2_KD,num_kept);
#pragma omp parallel for reduction(+:sum_n)
        for (int i=0;i<num_kept;i++){
            for (int j=num_kept;j<num_basis[n];j++){
                double omegan=0;
                omegan=(eigen_value[n][i]-eigen_value[n][j])*pow(Lambda,-1.0*(n-1-1)/2.0);
                //omegan=pow(Lambda,-1.0*(n-1-1)/2.0)*(eigen[n][i].eig_val-eigen[n][j].eig_val);
                //if (n < 6 ){
                ////cout << n << "  " << i << "  " << j <<  "  wn=  " << omegan << " A_nij=  " << c_dag_up_eigen[n][i][j]*c_up_eigen[n][j][i]*(exp(-1.0*Beta*pow(Lambda,-1.0*(n-1-1)/2.0)*eigen[n][i].eig_val_relat)+exp(-1.0*Beta*pow(Lambda,-1.0*(n-1-1)/2.0)*eigen[n][j].eig_val_relat)) <<  "  P= " << P_lorentz(freqency,-1.0*omegan) << endl;
                //cout << n << "  " << i << "  " << j <<  "  wn=  " << omegan << "  f=  " << freqency << "  P= " << (*pf)(freqency,-1.0*omegan) << endl;
                //}
                //sum_n=sum_n+c_up_eigen[n][i][j]*c_dag_up_eigen[n][j][k]*rho_red_temp[n][k][i]*(*pf)(freqency,-1.0*omegan)+rho_red_temp[n][i][k]*c_dag_up_eigen[n][k][j]*c_up_eigen[n][j][i]*(*pf)(freqency,omegan);
                //cout << "n= "  << n << "  i= " << i << "  j= " << j << "  " << (j-num_kept)+(dim_dot-1)*num_kept*i << "   " << matrix_prod1_DK[(j-num_kept)+(num_basis[n]-num_kept)*i] << "    " << i+num_kept*(j-num_kept) << "   " << matrix_prod2_KD[i+num_kept*(j-num_kept)] << endl;
                sum_n=sum_n+c_up_eigen[n][i][j]*matrix_prod1_DK[(j-num_kept)+(num_basis[n]-num_kept)*i]*(*pf)(freqency,-1.0*omegan)+matrix_prod2_KD[i+num_kept*(j-num_kept)]*c_up_eigen[n][j][i]*(*pf)(freqency,omegan);
                //cout << "x" << endl;
                //f_up << freqency << "  " << omegan << "  " << (*pf)(freqency,-1.0*omegan) << endl;
            }
        }
        sum=sum+sum_n;
        delete [] matrix_rho_KK;
        delete [] matrix_B;
        delete [] matrix_prod1_DK;
        delete [] matrix_prod2_KD;
    }
    return sum;
}
double dos3_down(double freqency)
{
    double sum=0;
    double P_h_unsmeared(double freqency);
    double P_h_smeared(double freqency);
    //ofstream f_down("down.dat");
    for (int n=n0;n<N_max-1;n++){
        double sum_n=0;
        double * matrix_rho_KK=new double [num_eigen_kept[n]*num_eigen_kept[n]];
        double * matrix_B=new double [(num_basis[n]-num_eigen_kept[n])*num_eigen_kept[n]];
        double * matrix_prod1_DK=new double [(num_basis[n]-num_eigen_kept[n])*num_eigen_kept[n]];
        double * matrix_prod2_KD=new double [num_eigen_kept[n]*(num_basis[n]-num_eigen_kept[n])];
        {int k=0;
            for (int i=0;i<num_eigen_kept[n];i++){
                for (int j=0;j<num_eigen_kept[n];j++){
                    matrix_rho_KK[k]=rho_red_temp[n][j][i];
                    k++;
                }
            }}
        {int k=0;
            for (int i=0;i<num_eigen_kept[n];i++){
                for (int j=num_eigen_kept[n];j<num_basis[n];j++){
                    matrix_B[k]=c_dag_down_eigen[n][j][i];
                    k++;
                }
            }}
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,num_basis[n]-num_eigen_kept[n],num_eigen_kept[n],num_eigen_kept[n],1,matrix_B,num_basis[n]-num_eigen_kept[n],matrix_rho_KK,num_eigen_kept[n],0,matrix_prod1_DK,num_basis[n]-num_eigen_kept[n]);
        {int k=0;
            for (int i=num_eigen_kept[n];i<num_basis[n];i++){
                for (int j=0;j<num_eigen_kept[n];j++){
                    matrix_B[k]=c_dag_down_eigen[n][j][i];
                    k++;
                }
            }}
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,num_eigen_kept[n],num_basis[n]-num_eigen_kept[n],num_eigen_kept[n],1,matrix_rho_KK,num_eigen_kept[n],matrix_B,num_eigen_kept[n],0,matrix_prod2_KD,num_eigen_kept[n]);
#pragma omp parallel for reduction(+:sum_n)
        for (int i=0;i<num_kept;i++){
            for (int j=num_kept;j<num_basis[n];j++){
                double omegan=0;
                omegan=(eigen_value[n][i]-eigen_value[n][j])*pow(Lambda,-1.0*(n-1-1)/2.0);
                sum_n=sum_n+c_down_eigen[n][i][j]*matrix_prod1_DK[(j-num_kept)+(num_basis[n]-num_kept)*i]*(*pf)(freqency,-1.0*omegan)+matrix_prod2_KD[i+num_kept*(j-num_kept)]*c_down_eigen[n][j][i]*(*pf)(freqency,omegan);
                /*
                 *for (int k=0;k<num_kept;k++){
                 *    sum_n=sum_n+c_down_eigen[n][i][j]*c_dag_down_eigen[n][j][k]*rho_red_temp[n][k][i]*(*pf)(freqency,-1.0*omegan)+  rho_red_temp[n][i][k]*c_dag_down_eigen[n][k][j]*c_down_eigen[n][j][i]*(*pf)(freqency,omegan);
                 *}
                 */
                //f_down << freqency << "  " << omegan << "  " << (*pf)(freqency,-1.0*omegan) << endl;
            }
        }
        sum=sum+sum_n;
        delete [] matrix_rho_KK;
        delete [] matrix_B;
        delete [] matrix_prod1_DK;
        delete [] matrix_prod2_KD;
    }
    return sum;
}

double P_K(double freqency, double omegan)//centered at omegan! not -1.0*omegan.
{
    double P_L(double,double);
    double P_G(double,double);
    double P_h(double);
    double ans;
    ans=P_L(freqency,omegan)*P_h(omegan)+P_G(freqency,omegan)*(1-P_h(omegan));// P_h(omegan)
    //ans=P_L(freqency,omegan)*P_h(freqency)+P_G(freqency,omegan)*(1-P_h(freqency));// P_h(freqency)
    return ans;
}
double P_L(double freqency,double omegan)
{
    double theta(double);
    double gamma=alpha/4.0;
    double ans;
    //ans=1.0/(sqrt(Pi)*alpha*fabs(omegan))*exp(-pow(log10(fabs(omegan/freqency))/alpha+gamma-alpha/2.0,2))*exp(-alpha*(gamma-alpha/4.0));//Log-Gaussian
    if (fabs(omegan)<=1e-13 || fabs(freqency)<=1e-13){
        ans=0;
    }else{
        ans=theta(freqency*omegan)/(sqrt(Pi)*alpha*fabs(omegan))*exp(-1.0*pow(log(fabs(omegan/freqency))/alpha+gamma-alpha/2.0,2.0))*exp(-1.0*alpha*(gamma-alpha/4.0));//eq. 1b final
        //ans=theta(freqency*omegan)/(sqrt(Pi)*alpha*fabs(freqency))*exp(-1.0*pow(log10(fabs(freqency/omegan))/alpha-gamma,2.0));//eq. 1b middle
    }
    return ans;
}
double P_G(double freqency,double omegan)
{
    double ans;
    ans=1.0/(omega0*sqrt(Pi))*exp(-1.0*pow((freqency-omegan)/omega0,2.0));
    return ans;
}
double theta(double freqency)
{
    double ans;
    if (freqency > 1e-13 ){
        ans=1.0;
    }else if (fabs(freqency)<=1e-13){
        ans=1.0/2.0;
        //}else if (freqency > 0 && freqency < 1e-20){
        //ans=1.0/2.0;
    }else {
        ans=0.0;
    }
    return ans;
}
double P_h(double omegan)
{
    double ans;
    if (smear && (fabs(omegan) < omega0)){
        ans=exp(-1.0*pow(log(fabs(omegan/omega0))/alpha,2.0));
    }else{
        ans=1.0;
    }
    return ans;
}
double P_sum_rule(double freqency, double omegan)
{
    double ans;
    ans=1.0;
    return ans;
}
