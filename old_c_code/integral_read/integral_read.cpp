#include <iostream>
#include "integrals.h"
#include <vector>
#include <string>
#include <H5Cpp.h>
#include <chrono>


using namespace H5;

// Function to read HDF5 file and output the tensor
std::vector<std::vector<std::vector<std::vector<float>>>> readHDF5_eri(hsize_t& dim1, hsize_t& dim2, hsize_t& dim3, hsize_t& dim4) {
    try {
        // Open the HDF5 file
        H5File file("integrals.H5", H5F_ACC_RDONLY);

        // Open the dataset (assume it is named "two_electron_integrals")
        DataSet dataset = file.openDataSet("two_electron_integrals");

        // Get the dataspace of the dataset
        DataSpace dataspace = dataset.getSpace();

        // Get the dimensions of the dataset
        hsize_t dims[4];
        dataspace.getSimpleExtentDims(dims, nullptr);

        dim1 = dims[0];
        dim2 = dims[1];
        dim3 = dims[2];
        dim4 = dims[3];

        int dim = dim1 * dim2 * dim3 * dim4;
        std::vector<float> tensor_flattened(dim);
        std::vector<std::vector<std::vector<std::vector<float>>>> tensor_eri(dim1,
            std::vector<std::vector<std::vector<float>>>(dim2,
            std::vector<std::vector<float>>(dim3,
            std::vector<float>(dim4))));

        // Read the dataset into the tensor
        dataset.read(tensor_flattened.data(), PredType::NATIVE_FLOAT);

        hsize_t index = 0;
        for (hsize_t i1 = 0; i1 < dim1; i1++) {
            for (hsize_t i2 = 0; i2 < dim2; i2++) {
                for (hsize_t i3 = 0; i3 < dim3; i3++) {
                    for (hsize_t i4 = 0; i4 < dim4; i4++) {
                        tensor_eri[i1][i2][i3][i4] = tensor_flattened[index++];
                    }
                }
            }
        }

    

        // Return the tensor
        return tensor_eri;

    } catch (FileIException& error) {
        std::cerr << "File I/O error" << std::endl;
        throw;
    } catch (DataSetIException& error) {
        std::cerr << "Dataset I/O error" << std::endl;
        throw;
    } catch (DataSpaceIException& error) {
        std::cerr << "DataSpace I/O error" << std::endl;
        throw;
    } catch (Exception& error) {
        std::cerr << "HDF5 error" << std::endl;
        throw;
    }
}

// Function to read HDF5 file and output the tensor
std::vector<std::vector<float>> readHDF5_h1e_ao(hsize_t& dim1, hsize_t& dim2) {
    try {
        // Open the HDF5 file
        H5File file("integrals.H5", H5F_ACC_RDONLY);

        // Open the dataset (assume it is named "two_electron_integrals")
        DataSet dataset_h1_ao = file.openDataSet("overlap_integrals");

        // Get the dataspace of the dataset
        DataSpace dataspace_h1_ao = dataset_h1_ao.getSpace();

        // Get the dimensions of the dataset
        hsize_t dims[2];
        dataspace_h1_ao.getSimpleExtentDims(dims, nullptr);

        dim1 = dims[0];
        dim2 = dims[1];

        int dim = dim1 * dim2;
        std::vector<float> tensor_flattened(dim);
        std::vector<std::vector<float>> tensor(dim1,std::vector<float> (dim2));

        // Read the dataset into the tensor
        dataset_h1_ao.read(tensor_flattened.data(), PredType::NATIVE_FLOAT);
    
        for (int i1 = 0; i1 < dim1; i1++) {
            for (int i2 = 0; i2 < dim2; i2++) {
                int index = dim1*i1 + i2;
                tensor[i1][i2] = tensor_flattened[index];
            }
        }
        


        // Return the tensor
        return tensor;

    } catch (FileIException& error) {
        std::cerr << "File I/O error" << std::endl;
        throw;
    } catch (DataSetIException& error) {
        std::cerr << "Dataset I/O error" << std::endl;
        throw;
    } catch (DataSpaceIException& error) {
        std::cerr << "DataSpace I/O error" << std::endl;
        throw;
    } catch (Exception& error) {
        std::cerr << "HDF5 error" << std::endl;
        throw;
    }
}

// Function to read HDF5 file and output the kinetic energy integral tensor
std::vector<std::vector<float>> readHDF5_h1e_kin(hsize_t& dim1, hsize_t& dim2) {
    try {
        // Open the HDF5 file
        H5File file("integrals.H5", H5F_ACC_RDONLY);

        // Open the dataset (assume it is named "two_electron_integrals")
        DataSet dataset_h1_kin_ao = file.openDataSet("kinetic_energy_integrals");

        // Get the dataspace of the dataset
        DataSpace dataspace_h1_kin_ao = dataset_h1_kin_ao.getSpace();

        // Get the dimensions of the dataset
        hsize_t dims[2];
        dataspace_h1_kin_ao.getSimpleExtentDims(dims, nullptr);

        dim1 = dims[0];
        dim2 = dims[1];

        int dim = dim1 * dim2;
        std::vector<float> tensor_flattened(dim);
        std::vector<std::vector<float>> tensor(dim1,std::vector<float> (dim2));

        // Read the dataset into the tensor
        dataset_h1_kin_ao.read(tensor_flattened.data(), PredType::NATIVE_FLOAT);
    
        for (int i1 = 0; i1 < dim1; i1++) {
            for (int i2 = 0; i2 < dim2; i2++) {
                int index = dim1*i1 + i2;
                tensor[i1][i2] = tensor_flattened[index];
            }
        }
        


        // Return the tensor
        return tensor;

    } catch (FileIException& error) {
        std::cerr << "File I/O error" << std::endl;
        throw;
    } catch (DataSetIException& error) {
        std::cerr << "Dataset I/O error" << std::endl;
        throw;
    } catch (DataSpaceIException& error) {
        std::cerr << "DataSpace I/O error" << std::endl;
        throw;
    } catch (Exception& error) {
        std::cerr << "HDF5 error" << std::endl;
        throw;
    }
}


// Function to read HDF5 file and output the nuclear energy integral tensor
std::vector<std::vector<float>> readHDF5_h1e_nuc(hsize_t& dim1, hsize_t& dim2) {
    try {
        // Open the HDF5 file
        H5File file("integrals.H5", H5F_ACC_RDONLY);

        // Open the dataset (assume it is named "two_electron_integrals")
        DataSet dataset_h1_nuc = file.openDataSet("nuclear_integrals");

        // Get the dataspace of the dataset
        DataSpace dataspace_h1_nuc = dataset_h1_nuc.getSpace();

        // Get the dimensions of the dataset
        hsize_t dims[2];
        dataspace_h1_nuc.getSimpleExtentDims(dims, nullptr);

        dim1 = dims[0];
        dim2 = dims[1];

        int dim = dim1 * dim2;
        std::vector<float> tensor_flattened(dim);
        std::vector<std::vector<float>> tensor(dim1,std::vector<float> (dim2));

        // Read the dataset into the tensor
        dataset_h1_nuc.read(tensor_flattened.data(), PredType::NATIVE_FLOAT);
    
        for (int i1 = 0; i1 < dim1; i1++) {
            for (int i2 = 0; i2 < dim2; i2++) {
                int index = dim1*i1 + i2;
                tensor[i1][i2] = tensor_flattened[index];
            }
        }
        


        // Return the tensor
        return tensor;

    } catch (FileIException& error) {
        std::cerr << "File I/O error" << std::endl;
        throw;
    } catch (DataSetIException& error) {
        std::cerr << "Dataset I/O error" << std::endl;
        throw;
    } catch (DataSpaceIException& error) {
        std::cerr << "DataSpace I/O error" << std::endl;
        throw;
    } catch (Exception& error) {
        std::cerr << "HDF5 error" << std::endl;
        throw;
    }
}