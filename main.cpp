#include "algebra.hpp"
#include "chrono.hpp"
using namespace algebra;

int main(){
    std::string filename="Matrix.mtx";
    Matrix<double,StorageOrder::Column> M(filename);
    Timings::Chrono chronometer;
    //M.compress();
    //M.print();
    //M.uncompress();
    //M.print();
    /*std::cout<<M(1,1)<<" "<<M(1,2)<<std::endl;
    std::complex<double> a=M(1,2);
    std::cout<<a<<std::endl;
    M(1,3)={2.5,0};
    M.print();*/

    std::vector<double> vec={1.0,1.0,1.0,1.0,1.0};
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
    //M.resize(5,6);
    //M.print();
    // to be fixed: operator () never calls const one

    std::cout<<"One norm: "<<M.norm<Norm::One><<std::endl;
    std::cout<<"One norm: "<<M.norm<Norm::One><<std::endl;
    std::cout<<"One norm: "<<M.norm<Norm::One><<std::endl;

    return 0;
}
