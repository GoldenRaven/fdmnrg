//2016/02/27
//2016/03/03只剩对角化与文件读入tn与en.
//2016/03/09完成对角化.
#include<iostream>
//void impurity(void);
void iterative_dia(void);
void func_wn(void);
void date_time(void);
void density_of_state(void);
void deallocate(void);
int main()
{
	using namespace std;
	cout << " Job started on: "; date_time();
	clock_t start,finish;
	start=clock();
	//impurity();
	iterative_dia();
	func_wn();
	density_of_state();
	deallocate();
	finish=clock();
	std::cout << "CPU time used:    " << (double)(finish - start)/CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "Job finished on:    ";date_time();
	return 0;
}
