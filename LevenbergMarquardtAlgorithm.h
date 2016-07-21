#ifndef H_LMA_H
#define H_LMA_H
#include "MultiVariableFunction.h"
#include "Matrix.h"
#include<vector>
class LevenbergMarquardtAlgorithm {
protected:
	//_foo ->foo member func of class
	//foo_ ->foo should not be deleted (created on heap outside class or is returned to user at some point)
	bool print = true;
	MultiVariableFunction* _func;
	unsigned int _nParams = 0;//no of params, passed in constructor
	unsigned int _nPoints = 0;//number of data points, passed in run() args
	unsigned int _currentIterations = 0;
	bool _canDelete = true;
	bool _initialised = false;
	float _lambda = 0.95f;//damping parameter
	float* _currentParams_;
	float* _initGuess_;//initial parameter guess; passed in constructor
	float* _partialDerivatives;//array of partial derivatives
	float** _data_;//passed in constructor
	float** _funcValues;
	float** _jacobian;
	Matrix *matrixJ, *matrixJT, *matrixDiagJTJ, *matrixYF,*matrixF, *matrixLHS, *matrixRHS, *matrixDelta;
	std::vector<Matrix*> matrixArr = { matrixJ,matrixJT,matrixDiagJTJ,matrixYF,matrixF,matrixLHS,matrixRHS,matrixDelta };
	void deleteAll() {
		if (!_canDelete) {
			return;
		}

		printStatus("DeleteAll");
		delete[] _partialDerivatives;
		delete[] _funcValues;
		delete2DArray(_jacobian, _nPoints);
		for (Matrix* a : matrixArr) {
			delete a;
		}

	}
	void delete2DArray(float** arr, unsigned int height) {
		printStatus("Delete2DArray");
		for (unsigned int i = 0; i < height; i++) {
			delete[] arr[i];
		}
		delete[] arr;
	}
	void initialise() {
		printStatus("initialise");
		_initialised = true;
		_partialDerivatives = new float[_nParams];
		_currentParams_ = new float[_nParams];
		_jacobian = new float*[_nPoints];
		_funcValues = new float*[1];
		_funcValues[0] = new float[_nPoints];
		for (unsigned int i = 0; i < _nPoints; i++) {
			_jacobian[i] = new float[_nParams];
		}
		for (unsigned int i = 0; i < _nParams; i++) {
			_currentParams_[i] = _initGuess_[i];
			std::cout << _currentParams_[i] << std::endl;
		}
	}
	void adjustParameters(float* delta) {
		printStatus("adjustParameters");
		for (unsigned int i = 0; i < _nParams; i++) {
			_currentParams_[i] += delta[i];
		}
	}
	void formJacobian() {
		printStatus("formJacobian");
		for (unsigned int i = 0; i < _nPoints; i++) {
			_funcValues[0][i] = _data_[0][i] - _func->getValue(_currentParams_, _data_[0][i]);
			_partialDerivatives = _func->getJacobian(_currentParams_, _partialDerivatives, _data_[0][i]);
			for (unsigned int j = 0; j < _nParams; j++) {
				_jacobian[i][j] = _partialDerivatives[j];
			}
		}
		std::cout << "forming jacobian done" << std::endl;
	}
	void runAlgorithm() {	
		printStatus("runAlgorithm");
		matrixJT = new Matrix(_nPoints, _nParams);		//
		matrixJ = new Matrix(_nParams, _nPoints);		//
		matrixDiagJTJ = new Matrix(_nParams, _nParams);	//formed from lambda*diag(JT*J)
		matrixYF = new Matrix(1, _nParams);
		matrixF = new Matrix(1, _nParams);
		matrixLHS = new Matrix(_nParams, _nParams);		//formed from (JT*J - lambda*diag(JT*J)).inverse
		matrixRHS = new Matrix(1, _nParams);
		matrixDelta = new Matrix(1, _nParams);
		
		_currentIterations = 0;
		while (true) {
			formJacobian();
			matrixLHS = getLHS();
			matrixLHS->printMatrix("LHS");
			//matrixRHS = getRHS();
			//matrixDelta->fromMatrixMultiplication(matrixLHS, matrixRHS);
			//adjustParameters(matrixDelta->getMatrix()[0]);
			_currentIterations++;
			break;
		}
	}
	Matrix* getLHS() {
		//return the left hand side of the LMA equation (JT * J  +  lambda*diag(JT*J))
		printStatus("getLHS");
		matrixJ->copyMatrix(_jacobian,_nParams,_nPoints);
		matrixJT->copyMatrixTranspose(matrixJ);
		matrixLHS->fromMatrixMultiplication(matrixJT, matrixJ);
		matrixDiagJTJ->copyMatrix(matrixLHS);
		matrixDiagJTJ->diagonalise();
		*matrixDiagJTJ *= _lambda;
		*matrixLHS += matrixDiagJTJ;
		matrixLHS->invert();
		return matrixLHS;
	}
	Matrix* getRHS() {
		//_y_ and _funcValues are float[1][_nPoints] (ie height 1), therefore form from matrix transpose
		printStatus("getRHS");
		matrixYF->copyMatrixTranspose(_funcValues, 1, _nPoints);
		matrixRHS->fromMatrixMultiplication(matrixJT, matrixYF);
		return matrixRHS;
	}
	void printStatus(std::string str) {
		if (print) {
			std::cout << str.c_str() << std::endl;
		}
	}
public:
	LevenbergMarquardtAlgorithm(MultiVariableFunction* func):_func(func),_nParams(func->getNParams()) {
	}
	LevenbergMarquardtAlgorithm() {

	}
	~LevenbergMarquardtAlgorithm() {
		printStatus("destructor");
		if (!_canDelete) {
			return;
		}
		deleteAll();
		_canDelete = false;
	}
	float* run(float** data, float* initGuess, unsigned int nPoints) {
		printStatus("run");
		_data_ = data;
		_initGuess_ = initGuess;
		_nPoints = nPoints;
		if (!_initialised) {
			initialise();
		}
		runAlgorithm();
		return _currentParams_;
	}
};
#endif