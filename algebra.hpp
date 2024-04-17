#ifndef ALGEBRA_HPP
#define ALGEBRA_HPP

#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <complex>
#include <fstream>
#include <sstream>
#include <map>

// definition of the StorageOrder enumerated class that will be useful for the template class Matrix 
enum class StorageOrder : unsigned int
{
  Row = 0,
  Column = 1
};

// definition of the comparator for the map of the Matrix class, so that it will be ordered differently based on the StorageOrder
// this is the default one
template <StorageOrder S>
struct MyTypeComparator {
    bool operator()(const std::array<std::size_t,2>& lhs, const std::array<std::size_t,2>& rhs) const {
            // if the row is equal, the ordering is based on the column, otherwise it is based on the row
            if(lhs[0]==rhs[0])
                return lhs[1] < rhs[1];
            else return lhs[0] < rhs[0];
        }
};

// this is the specialization for the column type matrix, since there are only two possible cases,
// it chooses the default if it is Row-wise, while it chooses the following if it is column-wise
template <>
struct MyTypeComparator<StorageOrder::Column>{
    bool operator()(const std::array<std::size_t,2>& lhs, const std::array<std::size_t,2>& rhs) const {
            // if the column is equal, the ordering is based on the row, otherwise it is based on the column
            if(lhs[1]==rhs[1])
                return lhs[0] < rhs[0];
            else return lhs[1] < rhs[1];
    }
};


namespace algebra{
    // definition of the template class Matrix
    template<typename T,StorageOrder S> class Matrix{
        public:
            // construntor that takes number of rows and number of columns
            Matrix(std::size_t n_rows,std::size_t n_cols):nrows(n_rows),ncols(n_cols){};
            // constructor that takes the filename of the Matrix Market format file
            Matrix(std::string filename){
                std::ifstream file(filename);
                if (!file.is_open()) {
                    std::cerr << "Error opening file: " << filename << std::endl;
                }

                std::string line;
                // Skip comment lines
                while (getline(file, line) && line[0] == '%');

                std::stringstream ss(line);
                std::size_t numRows, numCols, numNonZeros;
                if (!(ss >> numRows >> numCols >> numNonZeros)) {
                    std::cerr << "Invalid matrix market file format" << std::endl;
                }
                nrows=numRows;
                ncols=numCols;

                std::size_t row, col;
                T value;
                while (file >> row >> col >> value) {
                    std::array<std::size_t,2> pos={row-1,col-1};
                    elements[pos] = value; // Matrix Market is 1-indexed
                }

                file.close();

                std::cout << "Matrix read successfully from file: " << filename << std::endl;
            };
            // method taking the filename of the the Matrix Market format file
            void read_mtx(std::string filename){
                std::ifstream file(filename);
                if (!file.is_open()) {
                    std::cerr << "Error opening file: " << filename << std::endl;
                }

                std::string line;
                // Skip comment lines
                while (getline(file, line) && line[0] == '%');

                std::stringstream ss(line);
                std::size_t numRows, numCols, numNonZeros;
                if (!(ss >> numRows >> numCols >> numNonZeros)) {
                    std::cerr << "Invalid matrix market file format" << std::endl;
                }
                nrows=numRows;
                ncols=numCols;

                size_t row, col;
                T value;
                while (file >> row >> col >> value) {
                    std::array<std::size_t,2> pos={row-1,col-1};
                    elements[pos] = value; // Matrix Market is 1-indexed
                }

                file.close();

                std::cout << "Matrix read successfully from file: " << filename << std::endl;
                compressed=false;
                return;
            };
            // non-const call operator returns the reference to the value in the map
            T& operator()(std::size_t row, std::size_t col){
                // if values are out of bounds, it raises an error
                if(row>=nrows || col>=ncols)
                    std::cerr<<"Values out of bounds"<<std::endl;
                // if the matrix is compressed it checks if the value is stored, if not it raises an error, otherwise return the error
                if(compressed){
                    std::size_t index = find_elements(row,col);
                    if(index=values.size())
                        std::cerr<<"Cannot add elements if the matrix is compressed"<<std::endl;
                    return values[index];
                }
                // if the matrix is not compressed, it checks if element is stored, if not it stores a 0 in the called positions
                std::array<std::size_t,2> pos={row,col};
                if(!elements.find(std::array(row,col))){
                    elements[pos]=0;
                    }
                // returns the value
                return elements[pos];
            };
            // const call operator, return the value
            T operator()(std::size_t row, std::size_t col) const{
                // if values are out of bounds, it raises an error
                if(row>=nrows || col>=ncols)
                    std::cerr<<"Values out of bounds"<<std::endl;
                // if the matrix is not compressed, it checks if the element is stored, if so returns the value, 0 otherwise
                std::array<std::size_t,2> pos={row,col};
                if(!compressed){
                    if(!elements.find(pos))
                        return 0;
                    else
                        return elements[pos];
                }
                // if the matrix is compressed, it checks if the element is stored with another method, if so returns the value, 0 otherwise
                std::size_t index = find_elements(row,col);
                // if the element is not stored, it returns values.size() as default value for not stored indexes
                if(index==values.size())
                    return 0;
                else 
                    return values[index];
            };
            // method for checking if the matrix is compressed or not
            inline bool is_compressed(){return compressed;};
            // compress and uncompress methods to change the Matrix storage in the memory
            void compress(){
                if(compressed){
                    std::cout<<"Matrix already compressed"<<std::endl;
                    return;}
                if constexpr(S==StorageOrder::Row){
                    first_indexes.resize(nrows+1);
                    second_indexes.resize(elements.size());
                    values.resize(elements.size());
                    first_indexes[0]=0;
                    for(std::size_t i=0;i<nrows;i++){
                        std::array<std::size_t,2> pos={i,0};
                        auto lower=elements.lower_bound(pos);
                        pos={i,ncols-1};
                        auto upper=elements.upper_bound(pos);
                        std::size_t count=0;
                        for(auto k=lower;k!=upper;k++){
                            second_indexes.emplace_back(k->first[1]);
                            values.emplace_back(k->second);
                            count++;
                        }
                        first_indexes[i+1]=first_indexes[i]+count+1;
                    }
                }
                else{
                    first_indexes.resize(ncols+1);
                    second_indexes.resize(elements.size());
                    values.resize(elements.size());
                    first_indexes[0]=0;
                    for(std::size_t i=0;i<nrows;i++){
                        std::array<std::size_t,2> pos={0,i};
                        auto lower=elements.lower_bound(pos);
                        pos={nrows-1,i};
                        auto upper=elements.upper_bound(pos);
                        std::size_t count=0;
                        for(auto k=lower;k!=upper;k++){
                            second_indexes.emplace_back(k->first[0]);
                            values.emplace_back(k->second);
                            count++;
                        }
                        first_indexes[i+1]=first_indexes[i]+count+1;
                    }
                }
                elements.clear();
                compressed=true;
                std::cout<<"Matrix succesfully compressed"<<std::endl;
                return;
            };
            void uncompress(){
                if(!compressed){
                    std::cout<<"Matrix already uncompressed"<<std::endl;
                    return;
                }
                if constexpr(S=StorageOrder::Row){
                    for(std::size_t i=0;i<nrows;i++){
                        for(std::size_t j=first_indexes[i];j<first_indexes[i+1];j++){
                            std::array<std::size_t,2> pos={i,second_indexes[j]};
                            elements[pos]=values[j];
                        }
                    }
                }
                else{
                    for(std::size_t i=0;i<ncols;i++){
                        for(std::size_t j=first_indexes[i];j<first_indexes[i+1];j++){
                            std::array<std::size_t,2> pos={second_indexes[j],i};
                            elements[pos]=values[j];
                        }
                    }
                }
                compressed=false;
                first_indexes.clear();
                second_indexes.clear();
                values.clear();
                return;
            };
            // resize method for resizing the matrix
            void resize(std::size_t new_nrows,std::size_t new_ncols){
                if(new_nrows*new_ncols<nrows*ncols){
                    std::cerr<<"New dimensions are not feasible, too small dimensions"<<std::endl;
                }
                std::map<std::array<std::size_t,2>,T,MyTypeComparator<S>> elements_new;
                if constexpr(S==StorageOrder::Row){
                    for(auto i=elements.begin();i!=elements.end();i++){
                        std::size_t new_row=static_cast<std::size_t>(i->first[0]*i->first[1]/new_ncols);
                        std::size_t new_col=static_cast<std::size_t>(i->first[0]*i->first[1]%new_ncols);
                        std::array<std::size_t,2> pos={new_row,new_col};
                        elements_new[pos]=i->second;
                    }
                }
                else{
                    for(auto i=elements.begin();i!=elements.end();i++){
                        std::size_t new_col=static_cast<std::size_t>(i->first[0]*i->first[1]/new_nrows);
                        std::size_t new_row=static_cast<std::size_t>(i->first[0]*i->first[1]%new_nrows);
                        std::array<std::size_t,2> pos={new_row,new_col};
                        elements_new[pos]=i->second;
                    }        
                }
                elements.clear();
                elements=elements_new;
                nrows=new_nrows;
                ncols=new_ncols;
                return;
            };
            // * friend operator between a Matrix and a vector
            friend std::vector<T> operator*(Matrix<T,S>& M,const std::vector<T>& vec){
                if(vec.size()!=M.ncols)
                    std::cerr<<"Matrix-vector multiplication unfeasible, dimensions do not match"<<std::endl;
                std::vector<T> result(M.nrows);
                if constexpr(S==StorageOrder::Row){
                    if(M.compressed){
                        for (std::size_t i = 0; i < M.nrows; i++)
                        {
                            T element=0;
                            for (std::size_t j = M.first_indexes[i]; i < M.first_indexes[i+1]; j++)
                            {
                                element+=M.values[j]*vec[M.second_indexes[j]];
                            }
                            result[i]=element;
                        }
                    }
                    else{
                        for (std::size_t i = 0; i < M.nrows; i++)
                        {
                            T element=0;
                            std::array<std::size_t,2> pos={i,0};
                            auto lower=M.elements.lower_bound(pos);
                            pos={i,M.ncols-1};
                            auto upper=M.elements.upper_bound(pos);
                            for (auto j = lower; i != upper; j++)
                            {
                                element+=j->second*vec[j->first[1]];
                            }
                            result[i]=element;
                        }
                    } 
                }
                else{
                    for(std::size_t i=0;i<M.ncols;i++){
                        std::array<std::size_t,2> pos={0,i};
                        auto lower=M.elements.lower_bound(pos);
                        pos={M.nrows-1,i};
                        auto upper=M.elements.upper_bound(pos);
                        for(auto iter=lower;iter!=upper;iter++){
                            result[iter->first[0]]+=vec[iter->first[1]]*iter->second;
                        }
                    }

                }
                return result;
            };
            // method to print the matrix in mtx format
            void print(){
                if(!compressed){
                    std::cout<<nrows<<" "<<ncols<<" "<<elements.size()<<std::endl;
                    for(const auto& iter : elements){
                        std::cout<<iter.first[0]+1<<" "<<iter.first[1]+1 <<" "<<iter.second<<std::endl;
                    }
                }
                else{
                    std::cout<<nrows<<" "<<ncols<<" "<<values.size()<<std::endl;
                    if constexpr(S==StorageOrder::Row)
                    {
                        for (std::size_t i = 0; i < nrows; i++){
                            for (std::size_t j = first_indexes[i]; j < first_indexes[i+1]; j++)
                            {
                                std::cout<<i+1<<" "<<second_indexes[j]+1<<" "<<values[j]<<std::endl;
                            }
                        }
                    }
                    else{
                        for (std::size_t i = 0; i < ncols; i++){
                            for (std::size_t j = first_indexes[i]; j < first_indexes[i+1]; j++)
                            {
                                std::cout<<second_indexes[j]+1<<" "<<i+1<<" "<<values[j]<<std::endl;
                            }
                        }
                    }
                }
            };
        private:
            // number of rows and columns
            std::size_t nrows;
            std::size_t ncols;
            // map that stores the elements
            std::map<std::array<std::size_t,2>,T,MyTypeComparator<S>> elements;
            // boolean for keeping which storage is being used, false by default
            bool compressed=false;
            // vectors for the compressed 
            std::vector<std::size_t> first_indexes;
            std::vector<std::size_t> second_indexes;
            std::vector<T> values;
            // method for checking if an element is stored or not, used in the const call operator 
            std::size_t find_elements(std::size_t row,std::size_t col){
                // since I control the bounds before the call of this function, I don't need to do that again
                if(S==StorageOrder::Row){
                    // if I have zero 'non-zero' elements on the given row, for sure the element is not stored
                    if(first_indexes[row]==first_indexes[row+1])
                        return values.size(); 
                    // otherwise I look for the second index in the right span of indexes given by the first vector
                    for(std::size_t i=first_indexes[row];i<first_indexes[row+1];i++){
                        if(second_indexes[i]==col)
                            return i;
                    }
                }
                else{
                    // if I have zero 'non-zero' elements on the given column, for sure the element is not stored
                    if(first_indexes[col]==first_indexes[col+1])
                        return values.size(); 
                    // otherwise I look for the second index in the right span of indexes given by the first vector
                    for(std::size_t i=first_indexes[col];i<first_indexes[col+1];i++){
                        if(second_indexes[i]==row)
                            return i;
                }
                }
                return values.size();
            };
    };

};

/*
template<typename T> std::vector<T>& operator+=(std::vector<T>& lhs,const std::vector<T>& rhs);

template<typename T> std::vector<T>& operator*(const T lhs,const std::vector<T>& rhs);
*/

#endif