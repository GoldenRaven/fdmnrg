#include<fstream>
#include<iomanip>
#include<math.h>
using namespace std;
void uniform(void);
void logarithm(double,double, int,double);
void userdefined1(double a,double b,double c,double d, int N1, int N2, int N3, double lambda);
void userdefined2(double, int,double);
int main(void)
{
    //uniform();
    //userdefined2(1,300,0.92);
    logarithm(-1,1,560,0.92);
    //userdefined1(-0.5*0.02,0.5*0.02,60,0.7);
    //userdefined1(-1,-0.5,0.5,1,10,380,10,0.9);
    return 0;
}

void uniform(void)
{
    double a,b,c,d;
    int N1,N2,N3;
    double dx1,dx2,dx3;
    ofstream f_freq("freqency");
    a=-1.0;
    b=-0.08;
    c=0.08;
    d=1.0;
    N1=100;
    N2=400;
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
    f_freq << left << setw(25) << setprecision(15);
    for (int i=0;i<=N1;i++){
        f_freq << a+i*dx1 << endl;
    }
    for (int i=1;i<=N2;i++){
        f_freq << b+i*dx2 << endl;
    }
    for (int i=1;i<=N3;i++){
        f_freq << c+i*dx3 << endl;
    }
}
void logarithm(double a,double b, int N, double lambda)
{
    ofstream f_freq("freqency");
    f_freq << left << setw(25) << setprecision(15);
    for (int i=0;i<=N/2;i++){
        f_freq << a*pow(lambda,i) << endl;
    }
    f_freq << 0 << endl;
    for (int i=N/2;i>=0;i--){
        f_freq << b*pow(lambda,i) << endl;
    }
}
void userdefined1(double a,double b,double c,double d, int N1, int N2, int N3, double lambda)
{
    double dx1,dx3;
    ofstream f_freq("freqency");
    f_freq << left << setw(25) << setprecision(15);
    dx1=(b-a)/N1;
    dx3=(d-c)/N3;
    f_freq << left << setw(25) << setprecision(15);
    for (int i=0;i<=N1;i++){
        f_freq << a+i*dx1 << endl;
    }
    for (int i=1;i<=N2/2;i++){
        f_freq << b*pow(lambda,i) << endl;
    }
    f_freq << 0 << endl;
    for (int i=N2/2;i>=0;i--){
        f_freq << c*pow(lambda,i) << endl;
    }
    for (int i=1;i<=N3;i++){
        f_freq << c+i*dx3 << endl;
    }
}
void userdefined2(double b, int N, double lambda)
{
    ofstream f_freq("freqency");
    f_freq << left << setw(25) << setprecision(15);
    f_freq << 0 << endl;
    for (int i=N;i>=0;i--){
        f_freq << b*pow(lambda,i) << endl;
    }
}
