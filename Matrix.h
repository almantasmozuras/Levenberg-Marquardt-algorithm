#ifndef H_MATRIX_H
#define H_MATRIX_H
#include<iostream>
#include<cmath>

std::string MATRIX_ERR_INVERT_NONSQUARE = "Cannot invert non-square matrix";
std::string MATRIX_ERR_INVERT_DETERMINANT = "Determinant is zero; cannot invert";
std::string MATRIX_ERR_MULTIPLY_DIMENSION = "Cannot mutliply matrices; height_1 must equal width_2";
std::string MATRIX_ERR_COPY_DIMENSION = "Cannot copy matrix as dimensions are not the same";
std::string MATRIX_ERR_COPYTRANSPOSE_DIMENSION = "Cannot copy matrix transpose as dimensions are not the same";
std::string MATRIX_ERR_ADD_DIMENSION = "Cannot add matrices as the dimensions are not the same";
std::string MATRIX_ERR_DETERMINANT_NONSQUARE = "Cannot compute the determinant of a non-square matrix";

class Matrix {
private:
	bool _canDelete = true;
	double** _matrix;
	unsigned int _width;
	unsigned int _height;

	void deleteMatrix(double** matrix, unsigned int height) {
		//delete the given 2D array
		for (unsigned int i = 0; i < height; i++) {
			delete[] matrix[i];
		}
		delete[] matrix;
	}
	void deleteMatrix() {
		//delete _matrix
		if (!_canDelete) {
			return;
		}
		deleteMatrix(_matrix, _height);
	}
	void printMatrix(double** matrix, unsigned int width, unsigned int height) {
		for (unsigned int i = 0; i < height; i++) {
			for (unsigned int j = 0; j < width; j++) {
				std::cout << matrix[i][j] << "\t";
			}
			std::cout << std::endl;
		}
	}
	double getDeterminant(double** matrix, unsigned int size) {
		//recursive method to get the determinant of nxn matrix
		if (size == 1) {
			return matrix[0][0];
		}
		if (size == 2) {
			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
		}
		int sign = 1;
		double sum = 0;
		for (unsigned int i = 0; i < size; i++) {
			double** submatrix = getSubmatrix(matrix, i, 0, size, size);
			sum += (double) pow(-1,i)*matrix[0][i]*getDeterminant(submatrix,size-1);
			deleteMatrix(submatrix,size-1);
		}
		return sum;
	}
	
public:
	Matrix(unsigned int width, unsigned int height) :_width(width), _height(height) {
		_matrix = createMatrix(_width, _height);
	}
	Matrix():_width(1),_height(1) {
		_matrix = createMatrix(_width, _height);
	}
	Matrix(double** matrix, unsigned int width, unsigned int height) :_width(width), _height(height),_matrix(matrix) {
	}
	~Matrix() {
		if (_canDelete) {
			deleteMatrix();
		}
		_canDelete = false;
	}
	void printMatrix() {
		printMatrix(_matrix, _width, _height);
	}
	void printMatrix(std::string str) {
		std::cout << "======================\n" << str.c_str() << std::endl;
		std::cout<<"width: "<<_width<<", height: "<<_height<<"\n======================" << std::endl;
		printMatrix(_matrix, _width, _height);
	}
	void set(unsigned int x, unsigned int y, double val) {
		_matrix[y][x] = val;
	}
	void setMatrix(double** newMatrix, unsigned int newWidth, unsigned int newHeight) {
		deleteMatrix();
		_matrix = newMatrix;
		_width = newWidth;
		_height = newHeight;
	}
	void transpose() {
		//convert mxn to nxm matrix
		double** matrix2 = createMatrix(_height, _width);
		for (unsigned int i = 0; i < _width; i++) {
			for (unsigned int j = 0; j < _height; j++) {
				matrix2[i][j] = _matrix[j][i];
			}
		}
		deleteMatrix();
		_matrix = matrix2;
		unsigned int temp = _width;
		_width = _height;
		_height = temp;
	}
	void randomise(int max) {
		//
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] = (double)(rand() % max);
			}
		}
	}
	void fromMatrixMultiplication(Matrix *a, Matrix *b) {
		//set this->_matrix to the matrix created from multilplying axb
		if (a->_width != b->_height) {
			std::cout << MATRIX_ERR_MULTIPLY_DIMENSION.c_str() << std::endl;
			throw (std::string) MATRIX_ERR_MULTIPLY_DIMENSION;
		}
		double** newMatrix = _matrix;
		if (_width !=b->_width && _height !=a->_height) {
			newMatrix = createMatrix(b->_width, a->_height);
			_width = b->_width;
			_height = a->_height;
		}
		for (unsigned int i = 0; i < a->_height; i++) {
			for (unsigned int j = 0; j < b->_width; j++) {
				double sum = 0;
				for (unsigned int k = 0; k < a->_width; k++) {
					sum += a->_matrix[i][k] * b->_matrix[k][j];
				}
				newMatrix[i][j] = sum;
			}
		}
	}
	void copyMatrix(double** matrix, unsigned int width, unsigned int height) {
		//copies a 2D array to this->_matrix
		if (_width != width || _height != height) {
			std::cout << MATRIX_ERR_COPY_DIMENSION.c_str() << std::endl;
			throw (std::string) MATRIX_ERR_COPY_DIMENSION;
		}
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] = matrix[i][j];
			}
		}
	}
	void copyMatrix(Matrix *matrix) {
		copyMatrix(matrix->_matrix, matrix->_width, matrix->_height);
	}
	void copyMatrixTranspose(double** matrix, unsigned int width, unsigned int height) {
		//copies a 2D array to this->_matrix
		if (_width != height || _height != width) {
			std::cout << MATRIX_ERR_COPY_DIMENSION.c_str() << std::endl;
			throw (std::string) MATRIX_ERR_COPY_DIMENSION;
		}
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] = matrix[j][i];
			}
		}
	}
	void copyMatrixTranspose(Matrix *matrix) {
		copyMatrixTranspose(matrix->_matrix, matrix->_width, matrix->_height);
	}
	void diagonalise() {
		//set all elements except the leading diagonal to zero
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				if (i != j) {
					_matrix[i][j] = 0;
				}
			}
		}
	}
	void invert() {
		//invert this matrix
		// to get inverse, get cofactor matrix, divide by determinant and transpose

		double determinant = getDeterminant();
		if (_width != _height) {
			std::cout << MATRIX_ERR_INVERT_NONSQUARE.c_str() << std::endl;
			throw (std::string) MATRIX_ERR_INVERT_NONSQUARE;
		}
		if (determinant == 0) {
			std::cout << MATRIX_ERR_INVERT_DETERMINANT.c_str() << std::endl;
			throw (std::string) MATRIX_ERR_INVERT_DETERMINANT;
		}
		double** inverse = createMatrix(_width, _height);
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				//swapping i and j for array indexes as need to take transpose
				inverse[j][i] = (double) pow(-1,i+j) / determinant * getDeterminant(getSubmatrix(_matrix, j, i, _width, _height),_width-1);
			}
		}
		deleteMatrix();
		_matrix = inverse;
	}
	void addScalar(double val) {
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] += val;
			}
		}
	}
	void addMatrix(Matrix *other) {
		if (_height != other->_height || _width != other->_width) {
			std::cout << MATRIX_ERR_ADD_DIMENSION.c_str() << std::endl;
			throw (std::string) MATRIX_ERR_ADD_DIMENSION;
		}
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] += other->_matrix[i][j];
			}
		}
	}
	void addMatrix(double** matrix, unsigned int width, unsigned int height) {
		if (_height != height || _width != width) {
			std::cout << MATRIX_ERR_ADD_DIMENSION.c_str() << std::endl;
			throw (std::string) MATRIX_ERR_ADD_DIMENSION;
		}
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] += matrix[i][j];
			}
		}
	}
	void subtractMatrix(double** matrix, unsigned int width, unsigned int height) {
		if (_height != height || _width != width) {
			std::cout << MATRIX_ERR_ADD_DIMENSION.c_str() << std::endl;
			throw (std::string) MATRIX_ERR_ADD_DIMENSION;
		}
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] -= matrix[i][j];
			}
		}
	}
	void multiplyScalar(double val) {
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] *= val;
			}
		}
	}
	
	double** getMatrix() {
		return _matrix;
	}
	double** multiplyMatrix(Matrix *other) {
		//multiply this matrix by other (this matrix is on LHS)
		if (_width != other->_height) {
			std::cout << MATRIX_ERR_MULTIPLY_DIMENSION.c_str() << std::endl;
			throw (std::string) MATRIX_ERR_MULTIPLY_DIMENSION;
		}
		double** newMatrix = createMatrix(other->_width, _height);
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < other->_width; j++) {
				double sum = 0;
				for (unsigned int k = 0; k < _width; k++) {
					sum += _matrix[i][k] * other->_matrix[k][j];
				}
				newMatrix[i][j] = sum;
			}
		}
		return newMatrix;
	}
	double** getSubmatrix(double** matrix, unsigned int x, unsigned int y, unsigned int currentWidth, unsigned int currentHeight) {
		//return an (height-1)x(width-1) matrix
		//double
		double** submatrix = createMatrix(currentWidth - 1, currentHeight - 1);

		int yIndex = 0;
		for (unsigned int i = 0; i < currentHeight; i++) {
			if (i != y) {
				int xIndex = 0;
				for (unsigned int j = 0; j < currentWidth; j++) {
					if (j != x) {
						submatrix[yIndex][xIndex] = matrix[i][j];
						xIndex++;
					}
				}
				yIndex++;
			}
		}
		return submatrix;
	}
	double** createMatrix(unsigned int width, unsigned int height) {
		double** matrix = new double*[height];
		for (unsigned int i = 0; i < height; i++) {
			matrix[i] = new double[width];
			for (unsigned int j = 0; j < width; j++) {
				matrix[i][j] = 0;
			}
		}
		return matrix;
	}
	double get(unsigned int x, unsigned int y) {
		return _matrix[y][x];
	}
	double getDeterminant() {
		//return the determinant for this matrix; can only be done for square matrix
		if (_width != _height) {
			std::cout << MATRIX_ERR_DETERMINANT_NONSQUARE.c_str() << std::endl;
			throw (std::string) MATRIX_ERR_DETERMINANT_NONSQUARE;
		}
		return getDeterminant(_matrix,_width);
	}
	double getSum() {
		double sum = 0; 
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				sum += _matrix[i][j];
			}
		}
		return sum;
	}
	double getSquareSum() {
		double sum = 0;
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				sum += _matrix[i][j]*_matrix[i][j];
			}
		}
		return sum;
	}

	void operator+=(double val) {
		addScalar(val);
	}
	void operator+=(Matrix *other) {
		addMatrix(other);
	}
	void operator-=(double val) {
		addScalar(-val);
	}
	void operator/=(double val) {
		multiplyScalar(1 / val);
	}
	void operator*=(double val) {
		multiplyScalar(val);
	}
	void operator*=(Matrix *other) {
		unsigned int newWidth = 0; 
		unsigned int newHeight = 0;
		double** newMatrix = multiplyMatrix(other);
		deleteMatrix();
		_matrix = newMatrix;
		_width = other->_width;
	}
	double** operator*(Matrix *other) {
		return multiplyMatrix(other);
	}
};
#endif