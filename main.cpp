#include "algebra.hpp"
#include "chrono.hpp"
using namespace algebra;

int main(){
    std::string filename="lnsp_131.mtx";
    Matrix<double,StorageOrder::Row> M(filename);
    Timings::Chrono chronometer;

    std::vector<double> vec(131);
    chronometer.start();
    std::vector<double> result=M*vec;
    chronometer.stop();
    double time=chronometer.wallTime();
    std::cout<<"Time for uncompressed: "<<time<<" microseconds"<<std::endl;

    M.compress();

    chronometer.start();
    std::vector<double> result1=M*vec;
    chronometer.stop();
    double time1=chronometer.wallTime();
    std::cout<<"Time for compressed: "<<time1<<" microseconds"<<std::endl;

    std::cout<<"One norm: "<<M.norm<Norm::One>()<<std::endl;
    std::cout<<"Frobenius norm: "<<M.norm<Norm::Frobenius>()<<std::endl;
    std::cout<<"Inifnity norm: "<<M.norm<Norm::Infinity>()<<std::endl;

    Matrix<double,StorageOrder::Row> M2(131,1);
    std::vector<double> result2=M*M2;
    for (size_t i = 0; i < 131; i++)
    {
        std::cout<<result[i]<<" ";
    }


    return 0;
}
