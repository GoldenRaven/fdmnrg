#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include"setup.h"
//void impurity(void);
void iterative_dia_band(void);
void iterative_dia_total(void);
void func_wn(int);
void date_time(void);
void density_of_state(void);
//void deallocate(void);
void occu_imp(void);
double inner_energy(int);
double entropy_kB(int);
int main()
{
    using namespace std;
    cout << "Job started on: "; date_time();cout << endl;
    double S_band,S_total;
    clock_t start,finish;
    start=clock();
    ifstream f_id("job_id");
    f_id >> job_id;
    f_id.close();
    if (job_id==0){//0 for (occupation) calculation.
        cout << "chain of total: " << endl;
        iterative_dia_total();
        func_wn(1);
        if (occupation==1){
	    occu_imp();
	}
    }else if(job_id==1){//1 for (entropy) calculation.
        cout << "chain of band: " << endl;
        iterative_dia_band();//chain of band
        func_wn(0);
        S_band=entropy_kB(0);
        cout << "chain of total: " << endl;
        iterative_dia_total();
        func_wn(1);
        S_total=entropy_kB(1);
        cout << "Entropy: " << endl;
        cout << setw(30) << "  temperature: " << setprecision(15) << setw(30) << temperature << endl;
        cout << setw(30) << "  S_total:     " << setprecision(15) << setw(30) << S_total << endl;
        cout << setw(30) << "  S_band:      " << setprecision(15) << setw(30) << S_band  << endl;
        cout << setw(30) << "  S_imp:       " << setprecision(15) << setw(30) << S_total-S_band << endl;
        cout << "  Time leaved:    ";date_time();cout << endl;
        if (occupation==1){
	    occu_imp();
	}
    }else if(job_id==2){//2 for (imp_DOS, stm_DOS) calculation.
	cout << "chain of total: " << endl;
	iterative_dia_total();
	func_wn(1);
	density_of_state();
        if (occupation==1){
	    occu_imp();
	}
    }
    //deallocate();
    finish=clock();
    std::cout << "CPU time used:    " << (double)(finish - start)/CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "Job finished on:    ";date_time();
    return 0;
}
