#include<iomanip>
#include<iostream>
#include<fstream>
#include<math.h>
#include<omp.h>
#include"setup.h"
using namespace std;
void func_wn(void)
{
    cout << "func_wn():" << endl;
	double exp_Z(int n,int l);
	double * wn = new double [N_max+2];
	double sum=0;
	ofstream fwn("wn.dat");
	for (int n=n0;n<N_max+1;n++){
	    double suml=0;
		for (int l=num_kept;l<num_basis[n];l++){
		    suml=suml+exp_Z(n,l);
		}
		wn[n]=pow(4,N_max+1-n)*suml;
		fwn << setw(5) << n << setw(25) << scientific << setprecision(15) << wn[n] << endl;
		sum=sum+wn[n];
        cout << "  sum of wn= " << setprecision(16) << scientific << sum << endl;
	}
	for (int n=N_max+1;n<N_max+2;n++){
	    double suml=0;
		for (int l=0;l<num_basis[n];l++){
		    suml=suml+exp_Z(n,l);
		}
		wn[n]=pow(4,N_max+1-n)*suml;
		fwn << setw(5) << n << setw(25) << scientific << setprecision(15) << wn[n] << endl;
		sum=sum+wn[n];
        cout << "  sum of wn= " << setprecision(16) << scientific << sum << endl;
	}
    delete [] wn;
    cout << "  ";cout << "Time leaved:    ";date_time();
}

double exp_Z(int n,int l)
{
	double ans;
	double sumn=0;
	for (int n1=n0;n1<N_max+1;n1++){
		if (-1.0*Beta*(E_GS[n1]*pow(Lambda,-1.0*(n1-1-1)/2.0)-E_GS[n]*pow(Lambda,-1.0*(n-1-1)/2.0)) > 500){
			cout << n << "  " << n1 << "  " << scientific << -1.0*Beta*(E_GS[n1]*pow(Lambda,-1.0*(n1-1-1)/2.0)-E_GS[n]*pow(Lambda,-1.0*(n-1-1)/2.0)) << endl;
			ans=0;
			return ans;
		}else{
	        double suml1=0;
		    for (int l1=num_kept;l1<num_basis[n1];l1++){
		    	suml1=suml1+exp(-1.0*Beta*pow(Lambda,-1.0*(n1-1-1)/2.0)*eigen[n1][l1].eig_val_relat);
		    }
		    sumn=sumn+pow(4,N_max+1-n1)*exp(-1.0*Beta*(E_GS[n1]*pow(Lambda,-1.0*(n1-1-1)/2.0)-E_GS[n]*pow(Lambda,-1.0*(n-1-1)/2.0)))*suml1;
		}
	}
	for (int n1=N_max+1;n1<N_max+2;n1++){
		if (-1.0*Beta*(E_GS[n1]*pow(Lambda,-1.0*(n1-1-1)/2.0)-E_GS[n]*pow(Lambda,-1.0*(n-1-1)/2.0)) > 500){
			cout << n << "  " << n1 << "  " << -1.0*Beta*(E_GS[n1]*pow(Lambda,-1.0*(n1-1-1)/2.0)-E_GS[n]*pow(Lambda,-1.0*(n-1-1)/2.0)) << endl;
			ans=0;
			return ans;
		}else{
	        double suml1=0;
		    for (int l1=0;l1<num_basis[n1];l1++){
		    	suml1=suml1+exp(-1.0*Beta*pow(Lambda,-1.0*(n1-1-1)/2.0)*eigen[n1][l1].eig_val_relat);
		    }
		    sumn=sumn+pow(4,N_max+1-n1)*exp(-1.0*Beta*(E_GS[n1]*pow(Lambda,-1.0*(n1-1-1)/2.0)-E_GS[n]*pow(Lambda,-1.0*(n-1-1)/2.0)))*suml1;
		}
	}
	ans=exp(-1.0*Beta*pow(Lambda,-1.0*(n-1-1)/2.0)*eigen[n][l].eig_val_relat)/sumn;
	return ans;
}

/*
 *void func_wn(void)
 *{
 *    cout << "func_wn():" << endl;
 *    double * wn=new double [N_max+2];
 *    double * Zn=new double [N_max+2];
 *    for (int n=0;n<N_max+2;n++){
 *        Zn[n]=0;
 *        wn[n]=0;
 *    }
 *    for (int n=n0;n<N_max+1;n++){//partition function: \sum_{n_0}^{N_max} Zn
 *        for (int s=num_kept;s<num_basis[n];s++){
 *            Zn[n]=Zn[n]+exp(-1.0*pow(Lambda,-(n-1-1)/2.0)*Beta*(eigen[n][s].eig_val_relat));
 *            //Zn[n]=Zn[n]+exp(-1.0*Beta*(eigen[n][s].eig_val_relat));
 *            //Zn[n]=Zn[n]+exp(-1.0*pow(Lambda,-(n-1-1)/2.0)*Beta*(eigen[n][s].eig_val_relat+E_GS[n]));
 *            //cout << "Beta*eigen[n][s].eig_val= " << Beta*eigen[n][s].eig_val << endl;
 *        }
 *        cout << "Z" << n << " = " << Zn[n] << endl;
 *        Z=Z+Zn[n]*pow(dim_dot,N_max-(n-1));
 *    }
 *    for (int s=0;s<num_basis[N_max+1];s++){//Z_{N_max+1}
 *        Zn[N_max+1]=Zn[N_max+1]+exp(-1.0*pow(Lambda,-(N_max+1-1-1)/2.0)*Beta*(eigen[N_max+1][s].eig_val_relat));
 *        //Zn[N_max+1]=Zn[N_max+1]+exp(-1.0*pow(Lambda,-(N_max+1-1-1)/2.0)*Beta*(eigen[N_max+1][s].eig_val_relat+E_GS[N_max+1]));
 *        //cout << "sum in N_max+1= " << sum << endl;
 *    }
 *    cout << "Z" << N_max+1 << " = " << Zn[N_max+1] << endl;
 *    Z=Z+Zn[N_max+1];
 *    cout << "  " << "partition Z= " << Z << endl;
 *    double sum=0;
 *    ofstream fwn("wn.dat");
 *    for (int i=0;i<N_max+2;i++){
 *        wn[i]=pow(dim_dot,N_max+1-i)*Zn[i]/Z;
 *        fwn << setw(3) << i << setw(20) << wn[i] << endl;
 *        sum=sum+wn[i];
 *        cout << "  sum of wn= " << setprecision(16) << scientific << sum << endl;
 *    }
 *    delete [] wn;
 *    delete [] Zn;
 *    cout << "  ";cout << "Time leaved:    ";date_time();
 *}
 */
