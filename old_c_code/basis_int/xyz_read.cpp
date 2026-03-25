#include <iostream>
#include <string>
#include <vector>
#include <H5Cpp.h>
 

//#include "cint.h.in"

using namespace H5;

/*
template<typename T>
std::vector<T> read_hdf5_dataset(const std::string& filename, const std::string& dataset_name) {
    H5File file(filename, H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(dataset_name);
    
    DataSpace dataspace = dataset.getSpace();
    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, NULL);
    
    std::vector<T> data_basis(dims[0]);
    dataset.read(data_basis.data(), PredType::NATIVE_DOUBLE);
    return data_basis;
}
*/

int main() {

    //std::vector<int> bas = read_hdf5_dataset<int>("basis_info.H5", "bas");
    H5File file("basis_info.H5", H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet("bas");
    
    DataSpace dataspace = dataset.getSpace();
    
    // Get the dimensions of the dataset
    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims, NULL);
    
    // Create a vector to hold the data (make sure the type matches your HDF5 dataset)
    std::vector<int> bas(dims[0] * dims[1]); // Changed to double if the data type in HDF5 is double
    
    // Read the dataset into the vector
    dataset.read(bas.data(), PredType::NATIVE_INT);

    // Output the basis array elements
    std::cout << "Basis (bas) array elements:" << std::endl;
    std::cout << bas[4] << std::endl;
    std::cout << "Basis (bas) array elements:" << std::endl;
        for (size_t i = 0; i < dims[0]; i++) { // Loop over rows
            for (size_t j = 0; j < dims[1]; j++) { // Loop over columns
                std::cout << "bas[" << i << "][" << j << "]: " << bas[i * dims[1] + j] << std::endl;
            }
        }
// Open the dataset named "env"
    DataSet dataset_env = file.openDataSet("env");
    
    // Get the dataspace of the dataset
    DataSpace dataspace_env = dataset_env.getSpace();
    
    // Get the dimensions of the dataset
    hsize_t dims_env[1]; // env is 1D
    dataspace_env.getSimpleExtentDims(dims_env, NULL);
    
    // Create a vector to hold the data (64-bit floating points)
    std::vector<double> env(dims_env[0]); // Allocate based on the size
    
    // Read the dataset into the vector
    dataset_env.read(env.data(), PredType::NATIVE_DOUBLE); // Use NATIVE_DOUBLE for 64-bit

    // Output the env array elements
    std::cout << "Env (env) array elements:" << std::endl;
    for (size_t i = 0; i < dims_env[0]; i++) {
        std::cout << "env[" << i << "]: " << env[i] << std::endl;
    }
    //for (size_t i = 0; i < bas.size(); i++) {
    //    std::cout << "bas[" << i << "]: " << bas[i] << std::endl;
    //}

    DataSet dataset_coords = file.openDataSet("coords");
    
    DataSpace dataspace_coords = dataset_coords.getSpace();

    hsize_t dims_coords[2];
    
    // Get the dimensions of the dataset
    dataspace_coords.getSimpleExtentDims(dims_coords, NULL);
    
    // Create a vector to hold the data (make sure the type matches your HDF5 dataset)
    std::vector<double> coords(dims_coords[0] * dims_coords[1]); // Changed to double if the data type in HDF5 is double
    
    // Read the dataset into the vector
    dataset_coords.read(coords.data(), PredType::NATIVE_DOUBLE);

    // Output the basis array elements
    std::cout << "Basis (bas) array elements:" << std::endl;
    std::cout << "Basis (bas) array elements:" << std::endl;
        for (size_t i = 0; i < dims_coords[0]; i++) { // Loop over rows
            for (size_t j = 0; j < dims_coords[1]; j++) { // Loop over columns
                std::cout << "bas[" << i << "][" << j << "]: " << coords[i * dims[1] + j] << std::endl;
            }
        }

    return 0; // Return statement for main function
}