// Function to convert tensor to matrix
std::vector<std::vector<float>> tensorToMatrix(const std::vector<float>& tensor, hsize_t dim1, hsize_t dim2, hsize_t dim3, hsize_t dim4) {
    std::vector<std::vector<float>> matrix(dim1, std::vector<float>(dim1));

    for (int i1 = 0; i1 < dim1; i1++) {
        for (int i2 = 0; i2 < dim1; i2++) {
            int i3 = i1;
            int i4 = i2;
            size_t index = (i1 * dim2 * dim3 * dim4) + (i2 * dim3 * dim4) + (i3 * dim4) + i4;
            matrix[i1][i2] = tensor[index];
        }
    }

    return matrix;
