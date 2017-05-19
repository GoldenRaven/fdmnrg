//2016/02/27
//2016/03/03只剩对角化与文件读入tn与en.
//2016/03/09完成对角化.
#include<iostream>
#include<iomanip>
#include<cstdlib>
#include"setup.h"
//void impurity(void);
void iterative_dia(int);
void func_wn();
void date_time(void);
void density_of_state(void);
void deallocate(void);
void occu_imp(void);
double  inner_energy();
double entropy_kB();
//void entroy(void);
int main()
{
	using namespace std;
	cout << "Job started on: "; date_time();
	double S_band=0,S_total=0,E_band=0,E_total=0;
	clock_t start,finish;
	start=clock();
	//impurity();
	cout << "chain of band: " << endl;
	iterative_dia(0);
	func_wn();
	E_band=inner_energy();
	S_band=entropy_kB();
	//iterative_dia(1);
	//func_wn();
	//E_total=inner_energy();
	//S_total=entropy_kB();
	cout << "inner energy: " << setprecision(15) << setw(30) << temperature << setw(30) << E_total << "    " <<  setw(30) << E_band << "    " << setw(30) << E_total-E_band << endl;
	cout << "entropy: " << setprecision(15) << setw(30) << temperature << setw(30) << S_total << "    " <<  setw(30) << S_band << "    " << setw(30) << S_total-S_band << endl;
	//occu_imp();
	//density_of_state();
	//deallocate();
	finish=clock();
	std::cout << "CPU time used:    " << (double)(finish - start)/CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "Job finished on:    ";date_time();
	return 0;
}
