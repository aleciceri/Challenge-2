# Challenge-2
Second challenge for PACS course.

Implementation of Matrix template class.

In this repository there are these files:
- Makefile
- main.cpp
- Insp_131.mtx, Matrix Market file for the given test matrix
- Insp_131.mtx:Zone.Identifier, file downloaded with the Matrix Market matrix file
- Matrix.mtx, Matrix Market file, for writing smaller matrices to test the code in an easier way
- algebra.hpp, header-only file for the definition of the requested class

I will now describe the implementation of the Matrix class in algebra.hpp file.

Firstly, there is the defintion of the StorageOrder class needed for the template class Matrix, which requires a type for the elements and a StorageOrder class between Row and Column.

It is also present the defintion of a template ordering struct MyTypeComparator, which will be given to the map in the Matrix class so that the map will order the keys Row-wise or Column-wise depending on the given StorageOrder.

This struct has a default implementation and a specialized implementation for StorageOrder::Column. Since there are only 2 cases I know that in the StorageOrder case it will use the default implementation, while in the StorageOrder::Column it will use the specialized one, with no need for further specializations.

The private members of the class are:
- nrows: number of rows of the matrix
- ncols: number of columns of the matrix
- elements: map for the uncompressed storage of the non-zero values
- compressed: a boolean for storing if the matrix is now compressed or not
- first_indexes, second_indexes and values that are vectors for the storage of the matrix in the compressed case, according to CSR or CSC depending on the StorageOrder given
- find_elements: a function for finding if, in the compressed case, an element with given row and column is stored or not, used in the const call operator

Now I will describe the public methods.

The class, defined inside the namespace algebra as requested, has 2 constructors:
- the first takes only the number of rows and columns
- the second takes a string of a file and reads it as a Matrix Market file

It has also a method read_mtx that takes the filename and reads the file as a Matrix Market file that work in the very same way, with the detail that it clears the previous stored data since it overwrites the matrix. Also it changes the compressed boolean since by default it overwrites the Matrix and stores the new one in a uncompressed way.

I decided to make also a constructor even if only the method was requested because I think that it can be useful to have both, saving useless definition with the first constructor just to overwrite thw matrix the line after and still having the possibility with the method to read a Matrix Market file and overwrite the stored matrix.

