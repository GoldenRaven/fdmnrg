#include<fstream>
using namespace std;
int main()
{
	bool findpeak;
	double x_error;
   	double U;// in unit of D
	double Ed;//in unite of D
	double temperature;//k_B*T/D
	double Lambda;
	double alpha;
	int dim_imp=4;
	int dim_dot=4;
	int num_kept;
	double Beta_bar=0.6;//unit same as Beta!
	bool smear=false;
	bool unsmear=true;
	double omega0=temperature;
	bool Q=false;
	bool Q_Sz=false;
	bool N_up_N_down=true;

	ifstream f_findpeak("findpeak");
	ifstream f_x_error("x_error");
	ifstream f_U("U");
	ifstream f_Ed("Ed");
	ifstream f_temperature("temperature");
	ifstream f_Lambda("Lambda");
	ifstream f_alpha("alpha");
	ifstream f_num_kept("num_kept");
	f_findpeak >> findpeak;
	f_x_error >> x_error;
	f_U >> U;
	f_Ed >> Ed;
	f_temperature >> temperature;
	f_Lambda >> Lambda;
	f_alpha >> alpha;
	f_num_kept >> num_kept;

	ofstream f_input("input");
	f_input << findpeak    << endl;
	f_input << x_error     << endl;
	f_input << U           << endl;
	f_input << Ed          << endl;
	f_input << temperature << endl;
	f_input << Lambda      << endl;
	f_input << alpha       << endl;
	f_input << num_kept    << endl;
	f_input << smear       << endl;
	f_input << unsmear     << endl;
	f_input << omega0      << endl;
	f_input << dim_imp     << endl;
	f_input << dim_dot     << endl;
	f_input << Beta_bar    << endl;
	f_input << Q           << endl;
	f_input << Q_Sz        << endl;
	f_input << N_up_N_down << endl;
}
