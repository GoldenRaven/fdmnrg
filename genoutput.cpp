#include<iostream>
#include<iomanip>
#include<fstream>
#include"setup.h"
using namespace std;
void date_time(void);
void genoutput(void)
{
    cout << "  ##########################################################################################################################" << endl;
    cout << "  " << " job_id      = " << left << setw(30) << job_id     << " :job_id, 1 for entropy, 2 for DOS"     << endl;
    cout << "  " << " U           = " << left << setw(30) << U          << " :Cloumb U"                             << endl;
    cout << "  " << " Ed_up       = " << left << setw(30) << Ed_up      << " :Level up"                             << endl;
    cout << "  " << " Ed_down     = " << left << setw(30) << Ed_down    << " :Level down"                           << endl;
    cout << "  " << " num_kept    = " << left << setw(30) << num_kept   << " :Number of kept states"                << endl;
    cout << "  " << " Beta        = " << left << setw(30) << Beta       << " :Temperature Beta"                     << endl;
    cout << "  " << " Temperature = " << left << setw(30) << temperature<< " :Temperature Beta"                     << endl;
    cout << "  " << " Beta_bar    = " << left << setw(30) << Beta_bar   << " :Beta_bar"                             << endl;
    cout << "  " << " N_T         = " << left << setw(30) << N          << " :Dot corresponds to temperature"       << endl;
    cout << "  " << " N_max       = " << left << setw(30) << N_max      << " :Total number of dots"                 << endl;
    cout << "  " << " Lambda      = " << left << setw(30) << Lambda     << " :Logarithmic descritization parameter" << endl;
    cout << "  " << " smooth      = " << left << setw(30) << smooth     << " :smooth method"                        << endl;
    cout << "  " << " alpha       = " << left << setw(30) << alpha      << " :Width in log-Gaussian"                << endl;
    cout << "  " << " omega0      = " << left << setw(30) << omega0     << " :Width in Gaussian"                    << endl;
    cout << "  " << " Omega       = " << left << setw(30) << Omega      << " :cutoff in P_smooth"                   << endl;
    cout << "  " << " smear       = " << left << setw(30) << smear      << " :smear or unsmear"                     << endl;
    cout << "  " << " unsmear     = " << left << setw(30) << unsmear    << " :smear or unsmear"                     << endl;
    cout << "  " << " occupation  = " << left << setw(30) << occupation << " :calulate impurity occupation or not"  << endl;
    cout << "  " << " imp_dos     = " << left << setw(30) << imp_dos    << " :calulate imp DOS or not"              << endl;
    cout << "  " << " stm_dos     = " << left << setw(30) << stm_dos    << " :calulate stm DOS or not"              << endl;
    cout << "  " << " tc          = " << left << setw(30) << tc         << " :coupling tc"                          << endl;
    cout << "  " << " td          = " << left << setw(30) << td         << " :coupling td"                          << endl;
    cout << "  " << " N_max       = " << left << setw(30) << N_max      << " :max chain length"                     << endl;
    cout << "  " << " eig         = " << left << setw(30) << eig        << " :output energy flow or note"           << endl;
    if (Q){
        cout << "   Quantum number Charge(Q) is used in sorting." << endl;
    }else if(Q_Sz){
        cout << "   Quantum number Charge(Q) and total spin(Sz) is used in sorting." << endl;
    }else if(N_up_N_down){
        cout << "   Quantum number up number(N_up) and down number(N_down) is used in sorting." << endl;
    }
    cout << "  ##########################################################################################################################" << endl;
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
