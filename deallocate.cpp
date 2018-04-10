#include"setup.h"

void deallocate(void)
{
    for (int n=0;n<N_max;n++){
        for (int i=0;i<num_basis[n];i++){
            //delete [] vect[n][i];
            delete [] c_up_eigen[n][i];
            delete [] c_down_eigen[n][i];
        }
        //delete [] vect[n];
        //delete [] eigen[n];
        delete [] basis_ordered[n];
        delete [] c_up_eigen[n];
        delete [] c_down_eigen[n];
    }
    for (int n=0;n<N_max;n++){
        if (n==0){
            for (int i=0;i<dim_dot;i++){
                delete [] basis_kj[n][i];
            }
            delete [] basis_kj[n];
        }else{
            for (int i=0;i<num_eigen_kept[n-1];i++){
                delete [] basis_kj[n][i];
            }
            delete [] basis_kj[n];
        }
    }
    //delete [] vect;
    //delete [] eigen;
    delete [] basis_ordered;
    delete [] basis_kj;
    delete [] c_up_eigen;
    delete [] c_down_eigen;
    delete [] num_basis;
}
