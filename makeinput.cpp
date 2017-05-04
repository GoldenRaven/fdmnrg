#include<fstream>
using namespace std;
int main()
{
   	double U=0.5;// in unit of D
	double Ed_up=-0.25;//in unite of D
	double Ed_down=-0.25;//in unite of D
	double temperature=1.0e-7;//k_B*T/D
	double Lambda=1.8;
	double alpha=0.588;
	int dim_imp=4;
	int dim_dot=4;
	int num_kept=512;
	double Beta_bar=0.6;//unit same as Beta!
	bool smear=false;
	bool unsmear=true;
	double omega0=temperature;
	bool Q=false;
	bool Q_Sz=false;
	bool N_up_N_down=true;
	bool occupation=true;

	ifstream f_U("U");
	ifstream f_Ed_up("Ed_up");
	ifstream f_Ed_down("Ed_down");
	ifstream f_temperature("temperature");
	ifstream f_Lambda("Lambda");
	ifstream f_num_kept("num_kept");
	ifstream f_alpha("alpha");
	ifstream f_occupation("occupation");
	f_U >> U;
	f_Ed_up >> Ed_up;
	f_Ed_down >> Ed_down;
	f_temperature >> temperature;
	f_Lambda >> Lambda;
	f_alpha >> alpha;
	f_num_kept >> num_kept;
	f_occupation >> occupation;

	ofstream f_input_total("input_total");
	f_input_total << U           << endl;
	f_input_total << Ed_up       << endl;
	f_input_total << Ed_down     << endl;
	f_input_total << temperature << endl;
	f_input_total << Lambda      << endl;
	f_input_total << alpha       << endl;
	f_input_total << num_kept    << endl;
	f_input_total << smear       << endl;
	f_input_total << unsmear     << endl;
	f_input_total << omega0      << endl;
	f_input_total << dim_imp     << endl;
	f_input_total << dim_dot     << endl;
	f_input_total << Beta_bar    << endl;
	f_input_total << Q           << endl;
	f_input_total << Q_Sz        << endl;
	f_input_total << N_up_N_down << endl;
	f_input_total << occupation  << endl;

	ofstream f_input_band("input_band");
	f_input_band << temperature << endl;
	f_input_band << Lambda      << endl;
	f_input_band << num_kept    << endl;
	f_input_band << dim_imp     << endl;
	f_input_band << dim_dot     << endl;
	f_input_band << Beta_bar    << endl;
	f_input_band << Q           << endl;
	f_input_band << Q_Sz        << endl;
	f_input_band << N_up_N_down << endl;
}
