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

I decided to make also a constructor even if only the method was requested because I think that it can be useful to have both, saving useless definition with the first constructor just to overwrite thw matrix the line after and still having the possibility with the method to read a Matrix Market file and overwrite a stored matrix.

It follows the non-const call operator(), that returns a reference that can allow to change or, in case of uncompressed matrix, even add a new element. In fact if the matrix is compressed it only returns the reference to already existing stored elements, so if you call with non stored locations, it will give you a returned code and give you an error. Otherwise, if the matrix is compressed, it will create a new element in the map with the given location, set it 0 by default and then return you the reference.

In both compressed and uncompressed cases, if you call with already stored locations, it will return the already existing reference.

The next method is is_compressed(), which simply returns the private member compressed boolean.

In the const call operator(), on the other hand, it is made the same control to check if a location has been stored or not, but just to see if it is the case to return 0 or to find and return the right value. This difference is given by the fact that here you only need to return the value, so in the case of non stored value you don't need to create a new element which reference is to returned.

Then there are the methods compress() and uncompress(), which aim is to change the way to store the matrix from CSR (or CSC) to COOmap in uncompress() and the other way in compress(). In both method at the end the not used container is cleared to save us from using useless memory.

It follows the method resize, that takes two size_t inputs as new wanted dimension of the matrix and tries to change not only the dimension of the matrix, but also the locations accordingly. Firstly it check if the new dimension are fine, so that the total nomber of element is equal or greater, if not it raises an error because we don't shrink matrices, we only expand them by letting default 0 values in the new locations.

The new locations are computed based on the StorageOrder of the matrix object, filling the rows (or columns) with the zeros and non-zeros elements: for example if we increase the number of columns by 3 in a row ordered Matrix, in the first line those 3 elements will be the first 3 of the second line, the new second line starts from the fourth element of the current second line and so on. The same but opposite reasoning is made if we decrease the number of columns. It also holds for column ordered matrices, just flipping the rows and the columns in the process.

Then there is the friend call operator* between the Matrix and a vector. It firt checks if the dimensions are coherent, if not it raises an error. If the dimensions are correct, it computes the product following the given formulas (notice that in the case of column ordere Matrices, it add element-wise the j-th column multiplied by the j-th element of the vector, like in the formula, even if it is less readable than the row ordering case that is easier to be understood in the code).

I decided to do a print() method to print the Matrix in the Matrix Market format, which is very useful to have a visual representation of the result, still being given an easy readable format for sparse matrices.

