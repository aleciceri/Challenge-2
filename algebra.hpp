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

enum class Norm : unsigned int
{
  One = 0,
  Infinity = 1,
  Frobenius = 2
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
            Matrix(std::string filename);
            // method taking the filename of the the Matrix Market format file
            void read_mtx(std::string filename);
            // non-const call operator returns the reference to the value in the map
            T& operator()(std::size_t row, std::size_t col);
            // const call operator, return the value
            T operator()(std::size_t row, std::size_t col) const;
            // method for checking if the matrix is compressed or not
            bool is_compressed(){return compressed;};
            // compress and uncompress methods to change the Matrix storage in the memory
            void compress();
            void uncompress();
            // resize method for resizing the matrix
            void resize(std::size_t new_nrows,std::size_t new_ncols);
            // * friend operator between a Matrix and a vector
            friend std::vector<T> operator*(Matrix<T,S>& M,const std::vector<T>& vec){
                // control => check for right dimensions
                if(vec.size()!=M.ncols)
                    std::cerr<<"Matrix-vector multiplication unfeasible, dimensions do not match"<<std::endl;
                // definition of result vector that will be returned, with known dimension
                std::vector<T> result(M.nrows);
                // differentiation between the StorageOrder
                if constexpr(S==StorageOrder::Row){
                    // differentiation between compressed or not matrix, to know where to look for the data
                    if(M.compressed){
                        // loop over the rows
                        for (std::size_t i = 0; i < M.nrows; i++)
                        {
                            // value for the i-th element of the result vector
                            T element=0;
                            // loop over non-zero values of row i, for the row-vector vector multiplication
                            for (std::size_t j = M.first_indexes[i]; j < M.first_indexes[i+1]; j++)
                            {
                                // for every non-zero value, I sum the right increment
                                element+=M.values[j]*vec[M.second_indexes[j]];
                            }
                            // assignment of value
                            result[i]=element;
                        }
                    }
                    else{
                        // loop over rows
                        for (std::size_t i = 0; i < M.nrows; i++)
                        {
                            // value for the i-th element of the result vector
                            T element=0;
                            // defintion of lower and upper, iterators for firt element of the row and first iterator after the last element of the row
                            std::array<std::size_t,2> pos={i,0};
                            auto lower=M.elements.lower_bound(pos);
                            pos={i,M.ncols-1};
                            auto upper=M.elements.upper_bound(pos);
                            // loop over non-zero values of i-th row, with iterators, for the row vector multiplication
                            for (auto j = lower; j != upper; j++)
                            {
                                // for every non-zero-value, I sum the increment
                                element+=j->second*vec[j->first[1]];
                            }
                            // assignment of the value
                            result[i]=element;
                        }
                    } 
                }
                else{
                    if(!M.compressed){
                        // loop over columns
                        for(std::size_t i=0;i<M.ncols;i++){
                            // defintion of lower and upper, iterators for firt element of the column and first iterator after the last element of the column
                            std::array<std::size_t,2> pos={0,i};
                            auto lower=M.elements.lower_bound(pos);
                            pos={M.nrows-1,i};
                            auto upper=M.elements.upper_bound(pos);
                            // loop over non-zero elements of the i-th column
                            // I sum the i-th column vector times the i-th element of the vector element-wise, but still following the given formula
                            for(auto iter=lower;iter!=upper;iter++){
                                // summing each increment of the i-th column vector
                                result[iter->first[0]]+=vec[iter->first[1]]*iter->second;
                            }
                        }
                    }
                    else{
                        // loop over columns
                        for(std::size_t i=0;i<M.ncols;i++){
                            // loop over non-zero elements of column i
                            // I sum the i-th column vector times the i-th element of the vector element-wise, but still following the given formula
                            for (std::size_t j = M.first_indexes[i]; j < M.first_indexes[i+1]; j++)
                            {
                                // summing each increment of the i-th column vector
                                result[M.second_indexes[j]]+=M.values[j]*vec[i];
                            }
                        }
                    }
                }
                // return the result vector
                return result;
            };
            // method to print the matrix in mtx format
            void print();
            // template method for the norm, the default one will be the One Norm, then I'll specialize for the other 2 cases
            template<Norm N> 
            double norm<N>();
        private:
            // number of rows and columns
            std::size_t nrows;
            std::size_t ncols;
            // map that stores the elements
            std::map<std::array<std::size_t,2>,T,MyTypeComparator<S>> elements;
            // boolean for keeping which storage is being used, false by default
            bool compressed=false;
            // vectors for the compressed matrix
            std::vector<std::size_t> first_indexes;
            std::vector<std::size_t> second_indexes;
            std::vector<T> values;
            // method for checking if an element is stored or not in case of compressed matrix, used in the const call operator 
            std::size_t find_elements(std::size_t row,std::size_t col);
    };

};

// Following there are all the missing definitions of the method of the Matrix class
template<typename T, StorageOrder S>
algebra::Matrix<T,S>::Matrix(std::string filename){
                // definition of file stream
                std::ifstream file(filename);
                // if it is not possible to open the file, it creates an error
                if (!file.is_open()) {
                    std::cerr << "Error opening file: " << filename << std::endl;
                }

                std::string line;
                // Skip comment lines
                while (getline(file, line) && line[0] == '%');

                // string stream for reding the linea
                std::stringstream ss(line);
                std::size_t numRows, numCols, numNonZeros;
                // if the first line is not of the type "nrows ncols nonzeros", it is an error in the matrix format
                if (!(ss >> numRows >> numCols >> numNonZeros)) {
                    std::cerr << "Invalid matrix market file format" << std::endl;
                }
                // assignment of number of columns and rows
                nrows=numRows;
                ncols=numCols;

                std::size_t row, col;
                T value;
                // knowing the format, while I have the triple asignment means that I'm still reading the matrix
                // in the while condition I read the 3 elements
                while (file >> row >> col >> value) {
                    // I create an array for the map key, then I assign the value to the key
                    std::array<std::size_t,2> pos={row-1,col-1};
                    elements[pos] = value;
                }
                // before ending I have to close the file
                file.close();
                // Message for user
                std::cout << "Matrix read successfully from file: " << filename << std::endl;
            }

template<typename T, StorageOrder S>
void algebra::Matrix<T,S>::read_mtx(std::string filename){
                // I'm overwriting the matrix, so I have to erase the previous data
                if(compressed){
                    first_indexes.clear();
                    second_indexes.clear();
                    values.clear();
                }
                else{
                    elements.clear();
                }
                // definition of file stream
                std::ifstream file(filename);
                // if it is not possible to open the file, it creates an error
                if (!file.is_open()) {
                    std::cerr << "Error opening file: " << filename << std::endl;
                }

                std::string line;
                // Skip comment lines
                while (getline(file, line) && line[0] == '%');

                // string stream for reding the linea
                std::stringstream ss(line);
                std::size_t numRows, numCols, numNonZeros;
                // if the first line is not of the type "nrows ncols nonzeros", it is an error in the matrix format
                if (!(ss >> numRows >> numCols >> numNonZeros)) {
                    std::cerr << "Invalid matrix market file format" << std::endl;
                }
                // assignment of number of columns and rows
                nrows=numRows;
                ncols=numCols;

                std::size_t row, col;
                T value;
                // knowing the format, while I have the triple asignment means that I'm still reading the matrix
                // in the while condition I read the 3 elements
                while (file >> row >> col >> value) {
                    // I create an array for the map key, then I assign the value to the key
                    std::array<std::size_t,2> pos={row-1,col-1};
                    elements[pos] = value;
                }
                // before ending I have to close the file
                file.close();
                // Message for user
                std::cout << "Matrix read successfully from file: " << filename << std::endl;
                // I'm overwriting the matrix and putting it not compressed by default
                compressed=false;
                return;
            }

template<typename T, StorageOrder S>
T& algebra::Matrix<T,S>::operator()(std::size_t row, std::size_t col){
                // if values are out of bounds, it raises an error
                if(row>=nrows || col>=ncols)
                    std::cerr<<"Values out of bounds"<<std::endl;
                // if the matrix is compressed it checks if the value is stored, if not it raises an error, otherwise return the error
                if(compressed){
                    std::size_t index = find_elements(row,col);
                    if(index==values.size())
                        std::cerr<<"Cannot add elements if the matrix is compressed"<<std::endl;
                    return values[index];
                }
                // if the matrix is not compressed, it checks if element is stored, if not it stores a 0 in the called positions
                std::array<std::size_t,2> pos={row,col};
                if(elements.find(pos)==elements.end()){
                    elements[pos]=0;
                }
                // returns the value
                return elements[pos];
            }

template<typename T, StorageOrder S>
T algebra::Matrix<T,S>::operator()(std::size_t row, std::size_t col) const{
                // if values are out of bounds, it raises an error
                if(row>=nrows || col>=ncols)
                    std::cerr<<"Values out of bounds"<<std::endl;
                // if the matrix is not compressed, it checks if the element is stored, if so returns the value, 0 otherwise
                std::array<std::size_t,2> pos={row,col};
                if(!compressed){
                    if(elements.find(pos)==elements.end())
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
            }

template<typename T, StorageOrder S>
void algebra::Matrix<T,S>::compress(){
                // control => if the matrix is already compressed it tells the user and returns immediately
                if(compressed){
                    std::cout<<"Matrix already compressed"<<std::endl;
                    return;
                }
                // I have to differ the cases of Row-wise and Column-wise
                if constexpr(S==StorageOrder::Row){
                    // I know a priori the sizes needed, so I can already give them the needed space in memory
                    first_indexes.resize(nrows+1);
                    // I reserve and not resize the second_indexes and values so that I can use emplace back
                    second_indexes.reserve(elements.size());
                    values.reserve(elements.size());
                    // the index of the first-non zero of the first row is for sure 0, if there are
                    first_indexes[0]=0;
                    // loop over the rows
                    for(std::size_t i=0;i<nrows;i++){
                        // definition of upper and lower, the iterators for the first elements of the row i and the first iterator after the last element of row i
                        // since the map is given the ordering specialized on the StorageOrder, this works as written the comment above
                        std::array<std::size_t,2> pos={i,0};
                        auto lower=elements.lower_bound(pos);
                        pos={i,ncols-1};
                        auto upper=elements.upper_bound(pos);
                        // count for the non zero elements of the row
                        std::size_t count=0;
                        // for loop betweeen iterators
                        for(auto k=lower;k!=upper;k++){
                            // I store the column of the element and its value in the right vector
                            second_indexes.emplace_back(k->first[1]);
                            values.emplace_back(k->second);
                            // increase the count of nonzero elements in the line
                            count++;
                        }
                        // definition of next value of first_indexes according to the CSR
                        first_indexes[i+1]=first_indexes[i]+count;
                    }
                }
                else{
                    // I know a priori the sizes needed, so I can already give them the needed space in memory
                    first_indexes.resize(ncols+1);
                    // I reserve and not resize the second_indexes and values so that I can use emplace back
                    second_indexes.reserve(elements.size());
                    values.reserve(elements.size());
                    first_indexes[0]=0;
                    // loop over the columns
                    for(std::size_t i=0;i<ncols;i++){
                        // definition of upper and lower, the iterators for the first elements of the column i and the first iterator after the last element of column i
                        // since the map is given the ordering specialized on the StorageOrder, this works as written the comment above
                        std::array<std::size_t,2> pos={0,i};
                        auto lower=elements.lower_bound(pos);
                        pos={nrows-1,i};
                        auto upper=elements.upper_bound(pos);
                        // count for the non zero elements of the column
                        std::size_t count=0;
                        // loop between iterators
                        for(auto k=lower;k!=upper;k++){
                            // I store the row of the element and its value in the right vector
                            second_indexes.emplace_back(k->first[0]);
                            values.emplace_back(k->second);
                            // increase the count of nonzero elements in the line
                            count++;
                        }
                        // definition of next value of first_indexes according to the CSR
                        first_indexes[i+1]=first_indexes[i]+count;
                    }
                }
                // once I have taken all the data, I can clear elements
                elements.clear();
                // I change the value of compressed
                compressed=true;
                // message for user
                std::cout<<"Matrix succesfully compressed"<<std::endl;
                return;
            }

template<typename T, StorageOrder S>
void algebra::Matrix<T,S>::uncompress(){
                // control => if matrix is already uncompressed, It tells the user and stops
                if(!compressed){
                    std::cout<<"Matrix already uncompressed"<<std::endl;
                    return;
                }
                // differentiation of the function in case of Row-wise or Column-wise ordering
                if constexpr(S==StorageOrder::Row){
                    // loop over the rows
                    for(std::size_t i=0;i<nrows;i++){
                        // loop over non-zero elements of the row
                        for(std::size_t j=first_indexes[i];j<first_indexes[i+1];j++){
                            // definition of the key of the map, so row nd column of the element
                            std::array<std::size_t,2> pos={i,second_indexes[j]};
                            // assignment of the value
                            elements[pos]=values[j];
                        }
                    }
                }
                else{
                    // loop over the columns
                    for(std::size_t i=0;i<ncols;i++){
                        // loop over non-zero elements of the column
                        for(std::size_t j=first_indexes[i];j<first_indexes[i+1];j++){
                            // definition of the key of the map, so row and column of the element
                            std::array<std::size_t,2> pos={second_indexes[j],i};
                            // assignment of the value
                            elements[pos]=values[j];
                        }
                    }
                }
                // I change the compressed value since the matrix is now uncompressed
                compressed=false;
                // now I can free the memory space of the compressed matrix since I'm not using it
                first_indexes.clear();
                second_indexes.clear();
                values.clear();
                // message for the user
                std::cout<<"Matrix succesfuly uncompressed"<<std::endl;
                return;
            }

template<typename T, StorageOrder S>
void algebra::Matrix<T,S>::resize(std::size_t new_nrows, std::size_t new_ncols){
                // control => if the new size has less elements, the resize cannot be done, if there are more, the new elements are 0 by default
                if(new_nrows*new_ncols<nrows*ncols){
                    std::cerr<<"New dimensions are not feasible, too small dimensions"<<std::endl;
                }
                // this method can be called only if the matrix is uncompressed
                if(compressed==true){
                    std::cerr<<"Cannot resize if compressed"<<std::endl;
                }
                // new map for the new position of values in the resized matrix
                std::map<std::array<std::size_t,2>,T,MyTypeComparator<S>> elements_new;
                // differentiation between the two StorageOrder
                if constexpr(S==StorageOrder::Row){
                    // loop over the elements of the matrix
                    for(auto i=elements.begin();i!=elements.end();i++){
                        // computation of new row and new column index of the element
                        std::size_t new_row=static_cast<std::size_t>((i->first[0]*ncols + i->first[1])/new_ncols);
                        std::size_t new_col=static_cast<std::size_t>((i->first[0]*ncols + i->first[1])%new_ncols);
                        // definition of the key for the map
                        std::array<std::size_t,2> pos={new_row,new_col};
                        // assignment of the value in the new map
                        elements_new[pos]=i->second;
                    }
                }
                else{
                    // loop over the elements of the matrix
                    for(auto i=elements.begin();i!=elements.end();i++){
                        // computation of new row and new column index of the element
                        std::size_t new_col=static_cast<std::size_t>((i->first[0] + i->first[1]*nrows)/new_nrows);
                        std::size_t new_row=static_cast<std::size_t>((i->first[0] + i->first[1]*nrows)%new_nrows);
                        // definition of the key for the map
                        std::array<std::size_t,2> pos={new_row,new_col};
                        // assignment of the value in the new map
                        elements_new[pos]=i->second;
                    }        
                }
                // I can free the space with the old data and assign the new data with the new positions
                elements.clear();
                elements=elements_new;
                // assignment of new number of columns and new number of rows
                nrows=new_nrows;
                ncols=new_ncols;
                // message for user
                std::cout<<"Matrix succesfully resized"<<std::endl;
                return;
            }

template<typename T, StorageOrder S>
void algebra::Matrix<T,S>::print(){
                // differentiation over compressed or not
                if(!compressed){
                    // in uncompressed case, no differentiation for the StorageOrder of the matrix
                    // the ordering of the matrix already chooses the right order and there is not difference in the implementation
                    // first row of Matrix Market Format
                    std::cout<<nrows<<" "<<ncols<<" "<<elements.size()<<std::endl;
                    // loop over the elements of the map
                    for(const auto& iter : elements){
                        // printing "row column value"
                        std::cout<<iter.first[0]+1<<" "<<iter.first[1]+1 <<" "<<iter.second<<std::endl;
                    }
                }
                else{
                    // first row of Matrix Market Format
                    std::cout<<nrows<<" "<<ncols<<" "<<values.size()<<std::endl;
                    // differentiation for the StorageOrder
                    if constexpr(S==StorageOrder::Row)
                    {
                        // loop over rows
                        for (std::size_t i = 0; i < nrows; i++){
                            // loop over non-zero elements of the ith row
                            for (std::size_t j = first_indexes[i]; j < first_indexes[i+1]; j++)
                            {
                                // printing "row column value"
                                std::cout<<i+1<<" "<<second_indexes[j]+1<<" "<<values[j]<<std::endl;
                            }
                        }
                    }
                    else{
                        // loop over the columns
                        for (std::size_t i = 0; i < ncols; i++){
                            // loop over non-zero elements of ith column
                            for (std::size_t j = first_indexes[i]; j < first_indexes[i+1]; j++)
                            {
                                // printing "row column value"
                                std::cout<<second_indexes[j]+1<<" "<<i+1<<" "<<values[j]<<std::endl;
                            }
                        }
                    }
                }
            }

template<typename T, StorageOrder S>
std::size_t algebra::Matrix<T,S>::find_elements(std::size_t row, std::size_t col){
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
            }

template <typename T, StorageOrder S>
template <Norm N>
double algebra::Matrix<T,S>::norm<N>()
{
    // definition of norm, initialized as 0
    double norm=0;
    // differentiation with respect to the type of norm
    if constexpr (N==Norm::One){
        // differentiation with respect to the StorageOrder
        if constexpr (S==StorageOrder::Column){
            // differentiation with respect to the compressed or not
            if(compressed){
                // temporary variable
                double temp=0;
                // loop over the columns
                for(std::size_t i=0;i<ncols;i++){
                    // initialization of temporary for each column
                    temp=0;
                    // computation of sum of absolute values of elements of i-th column
                    for(std::size_t j=first_indexes[i];j<first_indexes[i+1];j++){
                        temp+=std::abs(values[j]);
                    }
                    // since I look for the max, if temp>norm, I store temp in norm and at the end I'll have the max
                    if(temp>norm)
                        norm=temp;
                }
            }
            else{
                // temporary value
                double temp=0
                // loop over the columns
                for(std::size_t i=0;i<ncols;i++){
                    // initialization of temporary for each column
                    temp=0;
                    // looking for iterator of first element of i-th column and first iterator after the last element of i-th column
                    std::array<std::size_t,2> pos={0,i};
                    auto lower=elements.lower_bound(pos);
                    pos={i,ncols-1};
                    auto upper=elements.upper_bound(pos);
                    // computation of sum of absolute values of elements of i-th column
                    for(auto iter=lower;iter!=upper;iter++){
                        temp+=std::abs(iter->second);
                    }
                    // since I look for the max, if temp>norm, I store temp in norm and at the end I'll have the max
                    if(temp>norm)
                        norm=temp;
                }
            }
        }
        else{
            // vector of temporaries, one for each column
            std::vector<double> temp(ncols);
            // differentiation with respect to compressed or not
            if(compressed){
                // loop over non zero values
                for(std::size_t j=0;j<values.size();j++){
                    // I add the absolute value of the element if the column-corresponent element of temp
                    temp[second_indexes[j]]+=std::abs(values[j]);
                }
            }
            else{
                // loop over non-zero elements
                for(auto iter : elements){
                    // I add the absolute value of the element if the column-corresponent element of temp
                    temp[iter->first[1]]+=std::abs(iter->second);
                }
            }
            // the norm is the maximum element in the vector of temporaries
            norm=std::max(temp);
        }
    }
    else if constexpr (N==Norm::Infinity){
        // differentiation with respect to compressed or not
        if(compressed){
            // loop over non-zero elements
            for(std::size_t i=0;i<values.size();i++){
                // if absolute value of element is greater than norm, I store it
                if(std::abs(values[i])>norm)
                    norm=std::abs(values[i]);
            }
        }
        else{
            // loop over non-zero elements
            for(auto iter : elements){
                // if absolute value of element is greater than norm, I store it
                if(std::abs(elements->second)>norm)
                    norm=std::abs(elements->second);
            }
        }
    }
    else{
        // differentiation with respect to compressed or not
        if(compressed){
            // loop over non-zero elements
            for(std::size_t i=0;i<values.size();i++){
                // for every non-zero element, I add the square of its absolute value
                norm+=square(std::abs(values[i]));
            }
        }
        else{
            // loop over non-zero elements
            for(auto iter : elements){
                // for every non-zero element, I add the square of its absolute value
                norm+=square(std::abs(elements->second));
            }
        }
        // store the square root of the sum of absolute values of all non zero elements
        norm=std::sqrt(norm);        
    }
    // return the norm
    return norm;
}

// function for square of a number, to allow me not to use std::pow or to call uselessy 2 times std::abs
double square(double num){
    return num*num;
}

#endif