// integrals.h
#ifndef INTEGRALS_H
#define INTEGRALS_H

#include <vector>
#include <H5Cpp.h>

using namespace H5;

// Function prototypes
std::vector<std::vector<std::vector<std::vector<float>>>> readHDF5_eri(hsize_t& dim1, hsize_t& dim2, hsize_t& dim3, hsize_t& dim4);
std::vector<std::vector<float>> readHDF5_h1e_ao(hsize_t& dim1, hsize_t& dim2);
std::vector<std::vector<float>> readHDF5_h1e_kin(hsize_t& dim1, hsize_t& dim2);
std::vector<std::vector<float>> readHDF5_h1e_nuc(hsize_t& dim1, hsize_t& dim2);

#endif // INTEGRALS_H
