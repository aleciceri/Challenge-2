#include "algebra.hpp"
using namespace algebra;

int main(){
    std::string filename="Matrix.mtx";
    Matrix<double,StorageOrder::Row> M(filename);

    return 0;
}
