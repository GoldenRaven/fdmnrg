//2016/02/27
//2016/03/03只剩对角化与文件读入tn与en.
//2016/03/09完成对角化.
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
void deallocate(void);
void occu_imp(void);
double inner_energy(int);
double entropy_kB(int);
//void entroy(void);
int main()
{
    using namespace std;
    cout << "Job started on: "; date_time();
    double S_band,S_total,E_band=0,E_total=0;
    clock_t start,finish;
    start=clock();
    ifstream f_id("job_id");
    f_id >> job_id;
    f_id.close();
    if (job_id==1){
        //impurity();
        cout << "chain of band: " << endl;
        iterative_dia_band();
        func_wn(0);
        S_band=entropy_kB(0);
        iterative_dia_total();
        func_wn(1);
        S_total=entropy_kB(1);
        cout << "entropy: " << setprecision(15) << setw(30) << temperature << setw(30) << S_total << "    " <<  setw(30) << S_band << "    " << setw(30) << S_total-S_band << endl;
        occu_imp();
    }else if(job_id==2){
        cout << "chain of band: " << endl;
        iterative_dia_total();
        func_wn(1);
        S_total=entropy_kB(1);
        cout << "entropy: " << setprecision(15) << setw(30) << temperature << setw(30) << S_total << "    " <<  setw(30) << S_band << "    " << setw(30) << S_total-S_band << endl;
        occu_imp();
        density_of_state();
    }
    //deallocate();
    finish=clock();
    std::cout << "CPU time used:    " << (double)(finish - start)/CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "Job finished on:    ";date_time();
    return 0;
}
