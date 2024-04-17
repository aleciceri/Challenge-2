#include "algebra.hpp"
using namespace algebra;

int main(){
    std::string filename="Matrix.mtx";
    Matrix<double,StorageOrder::Row> M(filename);
    std::cout<<M.is_compressed()<<std::endl;
    M.compress();
    M.print();

    return 0;
}
