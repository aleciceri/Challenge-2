#ifndef JACOBIANFACTORY_HPP
#define JACOBIANFACTORY_HPP

#include<vector>
#include<array>

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
            // non-const call operator returns the reference to the value in the map
            T& operator()(std::size_t row, std::size_t col);
            // const call operator, return the value
            T operator()(std::size_t row, std::size_t col) const;
            // method for checking if the matrix is compressed or not
            inline bool is_compressed(){return compressed;};
            // compress and uncompress methods to change the Matrix storage in the memory
            void compress();
            void uncompress();
            // resize method for resizing the matrix
            void resize(std::size_t new_nrows,std::size_t new_ncols);
            // * friend operator between a Matrix and a vector
            friend std::vector<T> operator*(Matrix<T,S>& M,const std::vector<T>& vec);
        private:
            // number of rows and columns
            std::size_t nrows;
            std::size_t ncols;
            // map that stores the elements
            std::map<std::array<std::size_t,2>,T,MyTypeComparator<S>> elements;
            // boolean for keeping which storage is being used, false by default
            bool compressed=false;
            // vectors for the compressed 
            std::vector<std::size_t> first_indexes(nrows+1);
            std::vector<std::size_t> second_indexes;
            std::vector<T> values;
            // method for checking if an element is stored or not, used in the const call operator 
            std::size_t find_elements(std::size_t row,std::size_t col);
    };

};


template<typename T> std::vector<T>& operator+=(std::vector<T>& lhs,const std::vector<T>& rhs);

template<typename T> std::vector<T>& operator*(const T lhs,const std::vector<T>& rhs);

#endif