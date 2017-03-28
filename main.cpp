//2016/02/27
//2016/03/03只剩对角化与文件读入tn与en.
//2016/03/09完成对角化.
#include<iostream>
#include<fstream>
#include<cstdlib>
#include"setup.h"
//void impurity(void);
void iterative_dia(void);
void func_wn(void);
void date_time(void);
void density_of_state(void);
void deallocate(void);
int main()
{
	using namespace std;
	cout << "Job started on: "; date_time();
	clock_t start,finish;
	start=clock();

	ifstream f_input("input");
	f_input >> findpeak    ;
	f_input >> x_error     ;
	f_input >> U           ;
	f_input >> Ed          ;
	f_input >> temperature ;
	f_input >> Lambda      ;
	f_input >> alpha       ;
	f_input >> num_kept    ;
	f_input >> smear       ;
	f_input >> unsmear     ;
	f_input >> omega0      ;
	f_input >> dim_imp     ;
	f_input >> dim_dot     ;
	f_input >> Beta_bar    ;
	f_input >> Q           ;
	f_input >> Q_Sz        ;
	f_input >> N_up_N_down ;

	//impurity();
	iterative_dia();
	func_wn();
	density_of_state();
	//deallocate();
	finish=clock();
	std::cout << "CPU time used:    " << (double)(finish - start)/CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "Job finished on:    ";date_time();
	return 0;
}
