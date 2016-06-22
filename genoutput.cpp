#include<iostream>
#include<iomanip>
#include<fstream>
#include"setup.h"
using namespace std;
void date_time(void);
void genoutput(void)
{
	cout << "  #############################################################################" << endl;
	cout << "  " << " U         = " << left << setw(15) << U          << " :Cloumb U" << endl;
	cout << "  " << " Ed        = " << left << setw(15) << Ed         << " :Level" << endl;
	cout << "  " << " num_kept  = " << left << setw(15) << num_kept   << " :Number of kept states" << endl;
	cout << "  " << " Beta      = " << left << setw(15) << Beta       << " :Temperature Beta" << endl;
	cout << "  " << " Beta_bar  = " << left << setw(15) << Beta_bar   << " :Beta_bar" << endl;
	cout << "  " << " N_T       = " << left << setw(15) << N          << " :Dot corresponds to temperature" << endl;
	cout << "  " << " N_max     = " << left << setw(15) << N_max      << " :Total number of dots" << endl;
	cout << "  " << " Lambda    = " << left << setw(15) << Lambda     << " :Logarithmic descritization parameter" <<  endl;
	cout << "  " << " alpha     = " << left << setw(15) << alpha      << " :Width in log-Gaussian" <<  endl;
	cout << "  " << " omega0    = " << left << setw(15) << omega0     << " :Width in Gaussian" <<  endl;
	cout << "  " << " smear     = " << left << setw(15) << smear << " :smear or unsmear" <<  endl;
	cout << "  #############################################################################" << endl;
	//char ch;
	//ifstream setup("../../setup.cpp");
	//while(setup.get(ch))
	//	cout << ch;
	//setup.close();
	//cout << "End of setup.h =========================================================>" << endl;
	//cout << "  Chain parameters in chain.dat ==========================================>" << endl;
	//ifstream chain("chain.dat");
	//while(chain.get(ch))
	//	cout << ch;
	//chain.close();
	//cout << "  End of chain.dat =======================================================>" << endl;
}
