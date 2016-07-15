#ifndef H_MATRIX_H
#define H_MATRIX_H
#include<iostream>
#include<cmath>

class Matrix {
private:
	bool _canDelete = true;
	float** _matrix;
	unsigned int _width;
	unsigned int _height;

	void deleteMatrix(float** matrix, unsigned int height) {
		for (unsigned int i = 0; i < height; i++) {
			delete[] matrix[i];
		}
		delete[] matrix;
	}
	void deleteMatrix() {
		if (!_canDelete) {
			return;
		}
		for (unsigned int i = 0; i < _height; i++) {
			delete[] _matrix[i];
		}
		delete[] _matrix;
	}
	void addScalar(float val) {
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] += val;
			}
		}
	}
	void addMatrix(Matrix &other) {
		if (_height != other._height || _width != other._width) {
			throw (std::string) "Cannot add matrices; not same dimensions";
		}
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] += other._matrix[i][j];
			}
		}
	}
	void multiplyScalar(float val) {
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] *= val;
			}
		}
	}
	float** multiplyMatrix(Matrix &other) {
		if (_width != other._height) {
			throw (std::string) "Cannot mutliply matrices; height_1 must equal width_2";
		}
		float** newMatrix = createMatrix(other._width, _height);
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < other._width; j++) {
				float sum = 0;
				for (unsigned int k = 0; k < _width; k++) {
					sum += _matrix[i][k] * other._matrix[k][j];
				}
				newMatrix[i][j] = sum;
			}
		}
		return newMatrix;
	}
	float getDeterminant(float** matrix, unsigned int size) {
		//std::cout << "size " << size << std::endl;
		if (size == 1) {
			return matrix[0][0];
		}
		if (size == 2) {
			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
		}
		int sign = 1;
		float sum = 0;
		for (unsigned int i = 0; i < size; i++) {
			float** submatrix = getSubmatrix(matrix, i, 0, size, size);
			sum += pow(-1,i)*matrix[0][i]*getDeterminant(submatrix,size-1);
			deleteMatrix(submatrix,size-1);
		}
		return sum;
	}
public:
	Matrix(unsigned int width, unsigned int height) :_width(width), _height(height) {
		_matrix = createMatrix(_width, _height);
	}
	Matrix(float** matrix, unsigned int width, unsigned int height) :_width(width), _height(height),_matrix(matrix) {
	}
	~Matrix() {
		deleteMatrix();
		_canDelete = false;
	}

	void printMatrix() {
		printMatrix(_matrix, _width, _height);
	}
	void printMatrix(float** matrix, unsigned int width, unsigned int height) {
		for (unsigned int i = 0; i < height; i++) {
			for (unsigned int j = 0; j < width; j++) {
				std::cout << matrix[i][j] << "\t";
			}
			std::cout << std::endl;
		}
	}
	void set(unsigned int x, unsigned int y, float val) {
		_matrix[y][x] = val;
	}
	void setMatrix(float** newMatrix, unsigned int newWidth, unsigned int newHeight) {
		deleteMatrix();
		_matrix = newMatrix;
		_width = newWidth;
		_height = newHeight;
	}
	void transpose() {
		float** matrix2 = createMatrix(_height, _width);
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
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] = (float)(rand() % max);
			}
		}
	}
	void fromMatrixMultiplication(Matrix &a, Matrix &b) {
		if (a._width != b._height) {
			throw (std::string) "Cannot mutliply matrices; height_1 must equal width_2";
		}
		float** newMatrix = _matrix;
		if (_width !=b._width && _height !=a._height) {
			newMatrix = createMatrix(b._width, a._height);
			_width = b._width;
			_height = a._height;
		}
		for (unsigned int i = 0; i < a._height; i++) {
			for (unsigned int j = 0; j < b._width; j++) {
				float sum = 0;
				for (unsigned int k = 0; k < a._width; k++) {
					sum += a._matrix[i][k] * b._matrix[k][j];
				}
				newMatrix[i][j] = sum;
			}
		}
	}
	void copyMatrix(float** matrix, unsigned int width, unsigned int height) {
		//copies a 2D array to this->_matrix
		if (_width != width || _height != height) {
			throw (std::string) "Cannot copy matrix as dimensions are not the same";
		}
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] = matrix[i][j];
			}
		}
	}
	void copyMatrix(Matrix &matrix) {
		copyMatrix(matrix._matrix, matrix._width, matrix._height);
	}
	void copyMatrixTranspose(float** matrix, unsigned int width, unsigned int height) {
		//copies a 2D array to this->_matrix
		if (_width != height || _height != width) {
			throw (std::string) "Cannot copy matrix transpose as dimensions are not the same";
		}
		for (unsigned int i = 0; i < _height; i++) {
			for (unsigned int j = 0; j < _width; j++) {
				_matrix[i][j] = matrix[j][i];
			}
		}
	}
	void copyMatrixTranspose(Matrix &matrix) {
		copyMatrixTranspose(matrix._matrix, matrix._width, matrix._height);
	}
	void diagonalise() {
		for (int i = 0; i < _height; i++) {
			for (int j = 0; j < _width; j++) {
				if (i != j) {
					_matrix[i][j] = 0;
				}
			}
		}
	}


	float** getSubmatrix(float** matrix, unsigned int x, unsigned int y, unsigned int currentWidth, unsigned int currentHeight) {
		//return an (height-1)x(width-1) matrix
		//float
		float** submatrix = createMatrix(currentWidth - 1, currentHeight - 1);

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
	float** createMatrix(unsigned int width, unsigned int height) {
		float** matrix = new float*[height];
		for (unsigned int i = 0; i < height; i++) {
			matrix[i] = new float[width];
			for (unsigned int j = 0; j < width; j++) {
				matrix[i][j] = 0;
			}
		}
		return matrix;
	}
	float get(unsigned int x, unsigned int y) {
		return _matrix[y][x];
	}
	float getDeterminant() {
		if (_width != _height) {
			throw (std::string) "Cannot take determinant of non-square matrix";
		}
		return getDeterminant(_matrix,_width);
	}

	void operator+=(float val) {
		addScalar(val);
	}
	void operator+=(Matrix &other) {
		addMatrix(other);
	}
	void operator-=(float val) {
		addScalar(-val);
	}
	void operator/=(float val) {
		multiplyScalar(1 / val);
	}
	void operator*=(float val) {
		multiplyScalar(val);
	}
	void operator*=(Matrix &other) {
		unsigned int newWidth = 0; 
		unsigned int newHeight = 0;
		float** newMatrix = multiplyMatrix(other);
		deleteMatrix();
		_matrix = newMatrix;
		_width = other._width;
	}
	float** operator*(Matrix &other) {
		return multiplyMatrix(other);//, other._width, other._height);
	}
};
#endif
