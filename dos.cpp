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
double **** rho;
void density_of_state(void)
{
    cout << "density_of_state: "; date_time();cout << endl;
    double dos1(double ***, double ***, double);
    double dos2(double ***, double ***, double);
    double dos3(double ***, double ***, double);
    void   rho_red(int);
    double P_LG(double,double);
    double P_G(double,double);
    double P_lorentz(double,double);
    double P_sum_rule(double,double);
    eigen_sigma=new double *** [N_max];
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
    for (int n=1;n<N_max;n++){
        eigen_sigma[n]=new double ** [dim_dot];
        for (int sigma=0;sigma<dim_dot;sigma++){
            eigen_sigma[n][sigma]=new double * [num_eigen_kept[n-1]];
            for (int k=0;k<num_eigen_kept[n-1];k++){
                eigen_sigma[n][sigma][k]=new double [num_basis[n]];
                for (int j=0;j<num_basis[n];j++){
                    eigen_sigma[n][sigma][k][j]=0;
                }
            }
        }
    }
    for (int n=1;n<N_max;n++){//initialization of eigen_sigma[n][sigma][k][j].
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
    eigen_value = new double * [N_max];
    for (int n=0;n<N_max;n++){
        eigen_value[n] = new double [num_basis[n]];
        for (int i=0;i<num_basis[n];i++){
            eigen_value[n][i] = eigen[n][i].eig_val;
        }
    }
    for (int n=0;n<N_max;n++){
        delete [] eigen[n];
    }
    delete [] eigen;
    rho_red_temp=new double ** [N_max-1];
    for (int n=0;n<N_max-1;n++){
        rho_red_temp[n]=new double * [num_kept];
        for (int i=0;i<num_kept;i++){
            rho_red_temp[n][i]=new double [num_kept];
            for (int j=0;j<num_kept;j++){
                rho_red_temp[n][i][j]=0;
            }
        }
    }
    rho = new double *** [N_max-1];
    for (int n=0;n<n0;n++){
        rho[n] = new double ** [1];
        for(int n2=0;n2<1;n2++){
            rho[n][n2] = new double * [1];
            for (int i=0;i<1;i++){
                rho[n][n2][i] = new double [1];
            }
        }
    }
    for (int n1=N_max-2;n1>=n0;n1--){
        rho_red(n1);
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
    cout << scientific << left << "                     ";
    cout << setw(14) << "freqency";
    cout << setw(14) << "dos1_up";
    cout << setw(14) << "dos2_up";
    cout << setw(14) << "dos3_up";
    cout << setw(14) << "dos1_down";
    cout << setw(14) << "dos2_down";
    cout << setw(14) << "dos3_down";
    cout << endl;
    double DOS1_UP,DOS2_UP,DOS3_UP,DOS1_DOWN,DOS2_DOWN,DOS3_DOWN;
    //test sum rule
    double freqency;
    freqency=1.0;
    pf=P_sum_rule;
    DOS1_UP=dos1(c_up_eigen, c_up_eigen, freqency);
    DOS2_UP=dos2(c_up_eigen, c_up_eigen, freqency);
    DOS3_UP=dos3(c_up_eigen, c_up_eigen, freqency);
    DOS1_DOWN=dos1(c_down_eigen, c_down_eigen, freqency);
    DOS2_DOWN=dos2(c_down_eigen, c_down_eigen, freqency);
    DOS3_DOWN=dos3(c_down_eigen, c_down_eigen, freqency);
    cout << "imp_dos sum_rule  ";
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
        ofstream f_imp_dos_smeared("imp_dos_smeared.dat");
        ofstream f_stm_dos_smeared("stm_dos_smeared.dat");
        while (!f_freqency.eof()){//???????????
            f_freqency >> freqency;
            if (imp_dos) {
                DOS1_UP=dos1(c_up_eigen, c_up_eigen, freqency);
                DOS2_UP=dos2(c_up_eigen, c_up_eigen, freqency);
                DOS3_UP=dos3(c_up_eigen, c_up_eigen, freqency);
                DOS1_DOWN=dos1(c_down_eigen, c_down_eigen, freqency);
                DOS2_DOWN=dos2(c_down_eigen, c_down_eigen, freqency);
                DOS3_DOWN=dos3(c_down_eigen, c_down_eigen, freqency);
                cout << "    imp          ";
                cout << setw(14) << setprecision(5) << freqency;
                cout << setw(14) << setprecision(5) << DOS1_UP;
                cout << setw(14) << setprecision(5) << DOS2_UP;
                cout << setw(14) << setprecision(5) << DOS3_UP;
                cout << setw(14) << setprecision(5) << DOS1_DOWN;
                cout << setw(14) << setprecision(5) << DOS2_DOWN;
                cout << setw(14) << setprecision(5) << DOS3_DOWN;
                date_time();cout << endl;
                //}
                f_imp_dos_smeared << left;
                f_imp_dos_smeared << setw(20) << setprecision(10) << freqency;
                f_imp_dos_smeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP;
                f_imp_dos_smeared << setw(20) << setprecision(10) << DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
                f_imp_dos_smeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP + DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
                f_imp_dos_smeared << endl;
            }
            if (stm_dos) {
                DOS1_UP=dos1(stm_f_up_eigen, stm_f_up_eigen, freqency);
                DOS2_UP=dos2(stm_f_up_eigen, stm_f_up_eigen, freqency);
                DOS3_UP=dos3(stm_f_up_eigen, stm_f_up_eigen, freqency);
                DOS1_DOWN=dos1(stm_f_down_eigen, stm_f_down_eigen, freqency);
                DOS2_DOWN=dos2(stm_f_down_eigen, stm_f_down_eigen, freqency);
                DOS3_DOWN=dos3(stm_f_down_eigen, stm_f_down_eigen, freqency);
                cout << "    stm          ";
                cout << setw(14) << setprecision(5) << freqency;
                cout << setw(14) << setprecision(5) << DOS1_UP;
                cout << setw(14) << setprecision(5) << DOS2_UP;
                cout << setw(14) << setprecision(5) << DOS3_UP;
                cout << setw(14) << setprecision(5) << DOS1_DOWN;
                cout << setw(14) << setprecision(5) << DOS2_DOWN;
                cout << setw(14) << setprecision(5) << DOS3_DOWN;
                date_time();cout << endl;
                //}
                f_stm_dos_smeared << left;
                f_stm_dos_smeared << setw(20) << setprecision(10) << freqency;
                f_stm_dos_smeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP;
                f_stm_dos_smeared << setw(20) << setprecision(10) << DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
                f_stm_dos_smeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP + DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
                f_stm_dos_smeared << endl;
            }
        }
        f_freqency.close();
        f_imp_dos_smeared.close();
        f_stm_dos_smeared.close();
    }
    if (unsmear){
        bool smear_tmp=smear;
        smear=false;
        ifstream f_freqency("freqency");
        ofstream f_imp_dos_unsmeared("imp_dos_unsmeared.dat");
        ofstream f_stm_dos_unsmeared("stm_dos_unsmeared.dat");
        while (!f_freqency.eof()){
            f_freqency >> freqency;
            if (imp_dos) {
                DOS1_UP=dos1(c_up_eigen, c_up_eigen, freqency);
                DOS2_UP=dos2(c_up_eigen, c_up_eigen, freqency);
                DOS3_UP=dos3(c_up_eigen, c_up_eigen, freqency);
                DOS1_DOWN=dos1(c_down_eigen, c_down_eigen, freqency);
                DOS2_DOWN=dos2(c_down_eigen, c_down_eigen, freqency);
                DOS3_DOWN=dos3(c_down_eigen, c_down_eigen, freqency);
                cout << "    imp          ";
                cout << setw(14) << setprecision(5) << freqency;
                cout << setw(14) << setprecision(5) << DOS1_UP;
                cout << setw(14) << setprecision(5) << DOS2_UP;
                cout << setw(14) << setprecision(5) << DOS3_UP;
                cout << setw(14) << setprecision(5) << DOS1_DOWN;
                cout << setw(14) << setprecision(5) << DOS2_DOWN;
                cout << setw(14) << setprecision(5) << DOS3_DOWN;
                date_time();cout << endl;
                f_imp_dos_unsmeared << left;
                f_imp_dos_unsmeared << setw(20) << setprecision(10) << freqency;
                f_imp_dos_unsmeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP;
                f_imp_dos_unsmeared << setw(20) << setprecision(10) << DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
                f_imp_dos_unsmeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP + DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
                f_imp_dos_unsmeared << endl;
            }
            if (stm_dos) {
                DOS1_UP=dos1(stm_f_up_eigen, stm_f_up_eigen, freqency);
                DOS2_UP=dos2(stm_f_up_eigen, stm_f_up_eigen, freqency);
                DOS3_UP=dos3(stm_f_up_eigen, stm_f_up_eigen, freqency);
                DOS1_DOWN=dos1(stm_f_down_eigen, stm_f_down_eigen, freqency);
                DOS2_DOWN=dos2(stm_f_down_eigen, stm_f_down_eigen, freqency);
                DOS3_DOWN=dos3(stm_f_down_eigen, stm_f_down_eigen, freqency);
                cout << "    stm          ";
                cout << setw(14) << setprecision(5) << freqency;
                cout << setw(14) << setprecision(5) << DOS1_UP;
                cout << setw(14) << setprecision(5) << DOS2_UP;
                cout << setw(14) << setprecision(5) << DOS3_UP;
                cout << setw(14) << setprecision(5) << DOS1_DOWN;
                cout << setw(14) << setprecision(5) << DOS2_DOWN;
                cout << setw(14) << setprecision(5) << DOS3_DOWN;
                date_time();cout << endl;
                //}
                f_stm_dos_unsmeared << left;
                f_stm_dos_unsmeared << setw(20) << setprecision(10) << freqency;
                f_stm_dos_unsmeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP;
                f_stm_dos_unsmeared << setw(20) << setprecision(10) << DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
                f_stm_dos_unsmeared << setw(20) << setprecision(10) << DOS1_UP+DOS2_UP+DOS3_UP + DOS1_DOWN+DOS2_DOWN+DOS3_DOWN;
                f_stm_dos_unsmeared << endl;
            }
        }
        smear=smear_tmp;
        //exit(0);
        f_freqency.close();
        f_imp_dos_unsmeared.close();
        f_stm_dos_unsmeared.close();
    }
    cout << "Time leaved:    ";date_time();cout << endl;
}

void rho_red(int n1)
{
    rho[n1] = new double ** [N_max-2-n1+1];
    for (int n2=0;n2<N_max-2-n1+1;n2++){
        rho[n1][n2] = new double * [num_kept];
        for (int i=0;i<num_kept;i++){
            rho[n1][n2][i] = new double [num_kept];
            for (int j=0;j<num_kept;j++){
                rho[n1][n2][i][j] = 0;
            }
        }
    }
    for (int n2=0;n2<N_max-2-n1;n2++){
        double ** temp1 = new double * [num_kept];
        double ** temp2 = new double * [num_kept];
        for (int i=0;i<num_kept;i++){
            temp1[i] = new double [num_kept];
            temp2[i] = new double [num_kept];
        }
        for (int sigma=0;sigma<dim_dot;sigma++){
            for (int i=0;i<num_kept;i++){
                for (int j=0;j<num_kept;j++){
                    temp1[i][j]=0;
                    temp2[i][j]=0;
                }
            }
#pragma omp parallel for
            for (int i=0;i<num_kept;i++){
                for (int j=0;j<num_kept;j++){
                    for (int k=0;k<num_kept;k++){
                        temp1[i][j]=temp1[i][j]+eigen_sigma[n1+1][sigma][i][k]*rho[n1+1][n2][k][j];
                    }
                }
            }
#pragma omp parallel for
            for (int i=0;i<num_kept;i++){
                for (int j=0;j<num_kept;j++){
                    for (int k=0;k<num_kept;k++){
                        temp2[i][j]=temp2[i][j]+temp1[i][k]*eigen_sigma[n1+1][sigma][j][k];
                    }
                }
            }
            for (int i=0;i<num_kept;i++){
                for (int j=0;j<num_kept;j++){
                    rho[n1][n2][i][j]=rho[n1][n2][i][j]+temp2[i][j];
                }
            }
        }
        for (int i=0;i<num_kept;i++){
            delete [] temp1[i];
            delete [] temp2[i];
        }
        delete [] temp1;
        delete [] temp2;
    }
    double ** temp1 = new double * [num_kept];
    for (int i=0;i<num_kept;i++){
        temp1[i] = new double [num_kept];
    }
    for (int sigma=0;sigma<dim_dot;sigma++){
        for (int i=0;i<num_kept;i++){
            for (int j=0;j<num_kept;j++){
                temp1[i][j]=0;
            }
        }
        int k_start;
        if (n1==N_max-2){
            k_start=0;
        }else{
            k_start=num_kept;
        }
#pragma omp parallel for
        for (int i=0;i<num_kept;i++){
            for (int j=0;j<num_kept;j++){
                for (int k=k_start;k<num_basis[n1];k++){
                    temp1[i][j]=temp1[i][j]+eigen_sigma[n1+1][sigma][i][k]*eigen_sigma[n1+1][sigma][j][k]*exp_z[n1+1][k-k_start];
                }
            }
        }
        for (int i=0;i<num_kept;i++){
            for (int j=0;j<num_kept;j++){
                rho[n1][N_max-2-n1][i][j]=rho[n1][N_max-2-n1][i][j]+temp1[i][j];
            }
        }
    }
    for (int i=0;i<num_kept;i++){
        delete [] temp1[i];
    }
    delete [] temp1;
    for (int n2=0;n2<N_max-2-n1+1;n2++){
        for (int i=0;i<num_kept;i++){
            for (int j=0;j<num_kept;j++){
                rho_red_temp[n1][i][j]=rho_red_temp[n1][i][j]+pow(dim_dot,n2)*rho[n1][n2][i][j];
            }
        }
    }
}
double dos1(double *** operatorA, double *** operatorB, double freqency)
{
    double sum=0;
    for (int n=n0;n<N_max-1;n++){
        double sum_n=0;
#pragma omp parallel for reduction(+:sum_n)
        for (int i=num_kept;i<num_basis[n];i++){
            for (int j=num_kept;j<num_basis[n];j++){
                double omegan=0;
                omegan=pow(Lambda,-1.0*(n-1-1)/2.0)*(eigen_value[n][i]-eigen_value[n][j]);
                sum_n=sum_n+operatorA[n][i][j]*operatorB[n][i][j]*(exp_z[n][i-num_kept]+exp_z[n][j-num_kept])*(*pf)(freqency,-1.0*omegan);
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
                sum_n=sum_n+operatorA[n][i][j]*operatorB[n][i][j]*(exp_z[n][i]+exp_z[n][j])*(*pf)(freqency,-1.0*omegan);
            }
        }
        sum=sum+pow(4,N_max-1-n)*sum_n;
    }
    return sum;
}
double dos2(double *** operatorA, double *** operatorB, double freqency)
{
    double sum=0;
    for (int n=n0;n<N_max-1;n++){
        double sum_n=0;
#pragma omp parallel for reduction(+:sum_n)
        for (int i=num_kept;i<num_basis[n];i++){
            for (int j=0;j<num_kept;j++){
                double omegan=0;
                omegan=(eigen_value[n][i]-eigen_value[n][j])*pow(Lambda,-1.0*(n-1-1)/2.0);
                sum_n=sum_n+operatorA[n][i][j]*operatorB[n][i][j]*exp_z[n][i-num_kept]*(*pf)(freqency,-1.0*omegan)+operatorB[n][j][i]*operatorA[n][j][i]*exp_z[n][i-num_kept]*(*pf)(freqency,omegan);
            }
        }
        sum=sum+sum_n*pow(4,N_max-1-n);
    }
    return sum;
}
double dos3(double *** operatorA, double *** operatorB, double freqency)
{
    double sum=0;
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
                    matrix_B[k]=operatorB[n][i][j];
                    k++;
                }
            }}
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,num_basis[n]-num_eigen_kept[n],num_eigen_kept[n],num_eigen_kept[n],1,matrix_B,num_basis[n]-num_eigen_kept[n],matrix_rho_KK,num_eigen_kept[n],0,matrix_prod1_DK,num_basis[n]-num_eigen_kept[n]);
        {int k=0;
            for (int i=num_kept;i<num_basis[n];i++){
                for (int j=0;j<num_kept;j++){
                    matrix_B[k]=operatorB[n][i][j];
                    k++;
                }
            }}
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,num_kept,num_basis[n]-num_kept,num_kept,1,matrix_rho_KK,num_kept,matrix_B,num_kept,0,matrix_prod2_KD,num_kept);
#pragma omp parallel for reduction(+:sum_n)
        for (int i=0;i<num_kept;i++){
            for (int j=num_kept;j<num_basis[n];j++){
                double omegan=0;
                omegan=(eigen_value[n][i]-eigen_value[n][j])*pow(Lambda,-1.0*(n-1-1)/2.0);
                sum_n=sum_n+operatorA[n][i][j]*matrix_prod1_DK[(j-num_kept)+(num_basis[n]-num_kept)*i]*(*pf)(freqency,-1.0*omegan)+matrix_prod2_KD[i+num_kept*(j-num_kept)]*operatorA[n][j][i]*(*pf)(freqency,omegan);
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
    if (freqency == 0 ){
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
    if (freqency > 0 ){
        ans=1.0;
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
