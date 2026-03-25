// main.cpp
#include "integrals.h"
#include <iostream>

int main() {
    hsize_t dim1, dim2, dim3, dim4;

    // Example usage of readHDF5_eri
    auto tensor_eri = readHDF5_eri(dim1, dim2, dim3, dim4);
    std::cout << "ERI tensor dimensions: " << dim1 << "x" << dim2 << "x" << dim3 << "x" << dim4 << std::endl;

    // Example usage of readHDF5_h1e_ao
    auto tensor_h1e_ao = readHDF5_h1e_ao(dim1, dim2);
    std::cout << "AO overlap tensor dimensions: " << dim1 << "x" << dim2 << std::endl;

    // Example usage of readHDF5_h1e_kin
    auto tensor_h1e_kin = readHDF5_h1e_kin(dim1, dim2);
    std::cout << "Kinetic energy tensor dimensions: " << dim1 << "x" << dim2 << std::endl;

    // Example usage of readHDF5_h1e_nuc
    auto tensor_h1e_nuc = readHDF5_h1e_nuc(dim1, dim2);
    std::cout << "Nuclear energy tensor dimensions: " << dim1 << "x" << dim2 << std::endl;

    return 0;
}
