#include <iostream>
#include "integral_read/integrals.h"
#include <H5Cpp.h>
#include <vector>
#include <string>
#include <chrono>
#include <Eigen/Dense>
#include <omp.h>
#include <tuple>


Eigen::MatrixXf Generate_new_overlap(Eigen::MatrixXf overlap_matrix, hsize_t dim1, hsize_t dim2) {

    //float occ_number = num_elec();

    //std::cout << "Occupation_number" << occ_number << std::endl;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> overlap_solver(overlap_matrix);

    Eigen::VectorXf eigenvalues = overlap_solver.eigenvalues().real();
    Eigen::MatrixXf eigen_vectors = overlap_solver.eigenvectors().real();

    Eigen::MatrixXf diagonal_mat(dim1,dim2);
    for (int i = 0; i < dim1; i++) {
        diagonal_mat(i,i) = 1.0/std::sqrt(eigenvalues[i]);
    }
    Eigen::MatrixXf new_overlap_matrix = eigen_vectors*diagonal_mat*eigen_vectors.transpose();

    return new_overlap_matrix;

}

Eigen::MatrixXf Generate_C0(Eigen::MatrixXf f, Eigen::MatrixXf new_overlap_matrix) {

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> fock_solver(f);

    Eigen::VectorXf eigenvalues_fock = fock_solver.eigenvalues().real();

    //f_T = L^T f_D L
    Eigen::MatrixXf eigen_vectors_fock = fock_solver.eigenvectors().real();

    //C0 = S^{-1/2}*L

    Eigen::MatrixXf C0 = new_overlap_matrix * eigen_vectors_fock;

    return C0;

}

Eigen::MatrixXf Generate_density_matrix(Eigen::MatrixXf new_overlap_matrix, float occ_number, Eigen::MatrixXf fock_matrix_0, hsize_t dim1, hsize_t dim2) {


    //Eigen::MatrixXf H_kin = readHDF5_h1e_kin(dim1,dim2);
    //Eigen::MatrixXf H_nuc = readHDF5_h1e_nuc(dim1,dim2);

   // Eigen::MatrixXf H_0 = H_kin + H_nuc;


    //F_T = S^{-1/2}.T * F * S^{-1/2}
    Eigen::MatrixXf f_T = new_overlap_matrix.transpose() * fock_matrix_0 * new_overlap_matrix; 

    Eigen::MatrixXf C0 = Generate_C0(f_T,new_overlap_matrix);
    
    
    Eigen::MatrixXf init_den(dim1,dim2);
    //std::cout << init_den << std::endl;
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < dim1; j++) {
        for (int k = 0; k < dim2; k++) {
            float pp = 0.0;
            for (int i = 0; i < occ_number; i++) {
                
                pp += C0(j,i)*C0(k,i);
                //sanity check
                //std::cout << pp << std::endl;
            }
            //std::cout << pp << std::endl;
            init_den(j,k) = pp;
        }
    }
    //std::cout << init_den << std::endl;
    return init_den;

}

Eigen::MatrixXf Generate_fock_matrix(Eigen::MatrixXf den_mat, Eigen::MatrixXf H_core, float occ_number, hsize_t dim1, hsize_t dim2, hsize_t dim3, hsize_t dim4) {
    //Construct the Fock_using the density matrix

    std::vector<std::vector<std::vector<std::vector<float>>>> eri = readHDF5_eri(dim1, dim2, dim3, dim4);

    Eigen::MatrixXf fock_matrix_int(dim1, dim2);

    #pragma omp parallel for collapse(2) // Parallelize the outer two loops

    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim1; j++) {
            float F = 0;
            for (int k = 0; k < dim1; k++) {
                for (int l = 0; l < dim1; l++) {
                    F += den_mat(k,l)*(2*eri[i][j][k][l] - eri[i][k][j][l]);
                }
            }
            fock_matrix_int(i,j) = H_core(i,j) + F;
        }
    }
    return fock_matrix_int;
}

float Compute_energy(Eigen::MatrixXf new_den_mat,  Eigen::MatrixXf new_fock_matrix, Eigen::MatrixXf H_core, hsize_t dim1) {
    float E = 0;
        for (int i = 0; i < dim1; i++) {
            for (int j = 0; j < dim1; j++) {
                E += new_den_mat(i,j)*(H_core(i,j) + new_fock_matrix(i,j));
            }
        }
    return E;
}


std::tuple<Eigen::MatrixXf, Eigen::MatrixXf>  SCF_procedure() {
    hsize_t dim1, dim2, dim3, dim4;
    Eigen::MatrixXf overlap_matrix = readHDF5_h1e_ao(dim1, dim2);
    Eigen::MatrixXf H_kin = readHDF5_h1e_kin(dim1,dim2);
    Eigen::MatrixXf H_nuc = readHDF5_h1e_nuc(dim1,dim2);
    float E_nuc = nuc_rep();

    Eigen::MatrixXf H_0 = H_kin + H_nuc;
    float occ_number = num_elec();


    Eigen::MatrixXf new_overlap = Generate_new_overlap(overlap_matrix, dim1 , dim2);

    
    //initialise the density matrix
    //Eigen::MatrixXf fock_old = H_0;
    Eigen::MatrixXf den_mat = Generate_density_matrix(new_overlap, occ_number , H_0, dim1, dim2);
    Eigen::MatrixXf fock_matrix = Generate_fock_matrix(den_mat, H_0, occ_number, dim1, dim2, dim3, dim4);

    float E_old = Compute_energy(den_mat, fock_matrix, H_0, dim1);
    float E0 = E_old+E_nuc;

    float tshld = 1;

    while (tshld > 1e-7) {
        Eigen::MatrixXf new_den_mat = Generate_density_matrix(new_overlap, occ_number , fock_matrix, dim1, dim2);
        Eigen::MatrixXf fock_matrix_new = Generate_fock_matrix(new_den_mat, H_0, occ_number, dim1, dim2, dim3, dim4);


        float E1 = Compute_energy(new_den_mat, fock_matrix_new, H_0, dim1);
        float E0 = E1+E_nuc;
        std::cout << "E_0\n" << E0 << std::endl;
        fock_matrix = fock_matrix_new;

        float E_diff = std::abs(E1-E_old);

        tshld = std::pow(E_diff, 2);

        E_old = E1;
    }

    Eigen::MatrixXf C0 = Generate_C0(fock_matrix, new_overlap);



    return std::make_tuple(fock_matrix, C0);
}

int main () {
    
    auto start = std::chrono::high_resolution_clock::now();
    Eigen::MatrixXf Fock_matrix, C0;

    std::tie(Fock_matrix, C0) = SCF_procedure();

    Eigen::MatrixXf fock_mo = C0.transpose() * Fock_matrix * C0;

    std::cout << Fock_matrix << std::endl;
    std::cout << std::endl;

    std::cout << fock_mo << std::endl;


    
auto end = std::chrono::high_resolution_clock::now();
    // Calculate the elapsed time

std::chrono::duration<double> elapsed = end - start;

// Print the elapsed time in seconds
std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
return 0;

}




