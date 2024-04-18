#include "algebra.hpp"
using namespace algebra;

int main(){
    std::string filename="Matrix.mtx";
    Matrix<std::complex<double>,StorageOrder::Row> M(filename);
    std::cout<<M.is_compressed()<<std::endl;
    M.compress();
    //M.print();
    M.uncompress();
    //M.print();
    std::cout<<M(1,1)<<" "<<M(1,2)<<std::endl;
    std::complex<double> a=M(1,2);
    std::cout<<a<<std::endl;
    M(1,3)={2.5,0};
    M.print();

    std::vector<std::complex<double>> vec={{1.0,0},{1.0,0},{1.0,0},{1.0,0},{1.0,0}};
    for (size_t i = 0; i < vec.size(); i++)
    {
        std::cout<<vec[i]<<" ";
    }
    std::cout<<std::endl;
    std::vector<std::complex<double>> result=M*vec;
    for (size_t i = 0; i < result.size(); i++)
    {
        std::cout<<result[i]<<" ";
    }
    std::cout<<std::endl;

    //M.resize(5,6);
    //M.print();
    // to be fixed: operator () never calls const one

    return 0;
}
