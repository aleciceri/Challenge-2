#include "algebra.hpp"

using namespace algebra;

// defintion of non const call operator
template<typename T,StorageOrder S>
T& Matrix<T,S>::operator()(std::size_t row, std::size_t col){
    // if values are out of bounds, it raises an error
    if(row>=nrows || col>=ncols)
        std::cerr<<'Values out of bounds'<<std::endl;
    // if the matrix is compressed it checks if the value is stored, if not it raises an error, otherwise return the error
    if(compressed){
        index = find_elements(row,col);
        if(index=values.size())
            std::cerr<<'Cannot add elements if the matrix is compressed'<<std::endl;
        return values[index];
    }
    // if the matrix is not compressed, it checks if element is stored, if not it stores a 0 in the called positions
    if(!elements.find(std::array(row,col)))
        elements[std::array(row,col)]=0;
    // returns the value
    return elements[std::array(row,col)];
}

template<typename T,StorageOrder S>
T Matrix<T,S>::operator()(std::size_t row, std::size_t col) const{
    // if values are out of bounds, it raises an error
    if(row>=nrows || col>=ncols)
        std::cerr<<'Values out of bounds'<<std::endl;
    // if the matrix is not compressed, it checks if the element is stored, if so returns the value, 0 otherwise
    if(!compressed){
        if(!elements.find(std::array(row,col)))
            return 0;
        else
            return elements[std::array(row,col)];
    }
    // if the matrix is compressed, it checks if the element is stored with another method, if so returns the value, 0 otherwise
    std::size_t index = find_elements(row,col);
    // if the element is not stored, it returns values.size() as default value for not stored indexes
    if(index==values.size())
        return 0;
    else 
        return values[index];
}

// definition of the helper method to find the index of a stored element or if the element is not stored in the compressed case
template<typename T>
std::size_t Matrix<T,StorageOrder::Row>::find_elements(std::size_t row,std::size_t col){
    if(first_indexes[row]==first_indexes[row+1])
        return values.size(); 
    for(std::size_t i=first_indexes[row];i<first_indexes[row+1];i++){
        if(second_index[i]==col)
            return i;
    }
    return values.size();
}

template<typename T>
std::size_t Matrix<T,StorageOrder::Column>::find_elements(std::size_t row,std::size_t col){
    if(first_indexes[col]==first_indexes[col+1])
        return values.size(); 
    for(std::size_t i=first_indexes[col];i<first_indexes[col+1];i++){
        if(second_index[i]==row)
            return i;
    }
    return values.size();
}












