#include<fstream>
using namespace std;
int main(void)
{
    ofstream f_freq("freqency");
	double freqency;
	double a,b,c,d;
	int N1,N2,N3;
	double dx1,dx2,dx3;
	a=-1.0;
	b=-0.08;
	c=0.08;
	d=1.0;
	N1=100;
	N2=300;
	N3=100;
	/*
	 *a=-3.0;
	 *b=-0.1;
	 *c=0.1;
	 *d=3.0;
	 *N1=200;
	 *N2=400;
	 *N3=200;
	 */
	/*
	 *a=-7.42e-4;
	 *b=-3.71e-4;
	 *c=3.71e-4;
	 *d=7.42e-4;
	 *N1=15;
	 *N2=60;
	 *N3=15;
	 */
	dx1=(b-a)/N1;
	dx2=(c-b)/N2;
	dx3=(d-c)/N3;
	for (int i=0;i<=N1;i++){
		f_freq << a+i*dx1 << endl;
	}
	for (int i=1;i<=N2;i++){
		f_freq << b+i*dx2 << endl;
	}
	for (int i=1;i<=N3;i++){
		f_freq << c+i*dx3 << endl;
	}
	return 0;
}

