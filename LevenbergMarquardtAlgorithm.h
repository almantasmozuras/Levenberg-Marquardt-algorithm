#ifndef H_L MA_H
#define H_LMA_H
#include "MultiVariableFunction.h"
#include "Matrix.h"
class LevenbergMarquardtAlgorithm {
protected:
	//_foo ->foo member func of class
	//foo_ ->foo should not be deleted (created on heap outside class or is returned to user at some point)
	MultiVariableFunction _func;
	unsigned int _nParams;//no of params, passed in constructor
	unsigned int _nPoints;//number of data points, passed in run() args
	unsigned int _currentIterations = 0;
	bool _canDelete = true;
	float _lambda = 0.95f;//damping parameter
	float* _currentParams_;
	float* _delta;//how much each parameter should be changed by
	float* _initGuess_;//initial parameter guess; passed in constructor
	float* _partialDerivatives;//array of 
	float* _funcValues;
	float** _data_;//passed in constructor
	float** _jacobian;
	void deleteAll() {
		if (!_canDelete) {
			return;
		}
		delete[] _delta;
		delete[] _partialDerivatives;
		delete[] _funcValues;
		delete2DArray(_jacobian, _nParams);
	}
	void delete2DArray(float** arr, unsigned int height) {
		for (unsigned int i = 0; i < height; i++) {
			delete[] arr[i];
		}
		delete arr;
	}
	void initialise() {
		_partialDerivatives = new float[_nParams];
		_currentParams_ = new float[_nParams];
		_jacobian = new float*[_nParams];
		_funcValues = new float[_nPoints];
		_delta = new float[_nParams];
		for (unsigned int i = 0; i < _nParams; i++) {
			_jacobian[i] = new float[_nPoints];
			_partialDerivatives[i] = 0;
			_currentParams_[i] = _initGuess_[i];
		}
	}
	void adjustParameters() {
		for (unsigned int i = 0; i < _nParams; i++) {
			_currentParams_[i] += _delta[i];
		}
	}
	void formJacobian() {
		for (unsigned int i = 0; i < _nPoints; i++) {
			_funcValues[i] = _func.getValue(_currentParams_, _data_[0][i]);
			_partialDerivatives = _func.getJacobian(_currentParams_, _partialDerivatives, _data_[0][i]);
			for (unsigned int j = 0; j < _nParams; j++) {
				_jacobian[j][i] = _partialDerivatives[j];
			}
		}
	}
	void runAlgorithm() {
		Matrix diagJJT(_nParams, _nParams);	//formed from lambda*diag(J*JT)
		Matrix LHS(_nParams, _nParams);		//formed from (JJT - lambda*diag(JJT)).inverse
		Matrix J(_nPoints, _nParams);		//
		Matrix JT(_nParams, _nPoints);		//
		Matrix JJT(_nParams, _nParams);
		Matrix F(1, _nParams);				//
		_currentIterations = 0;
		while (true) {

			formJacobian();

			_currentIterations++;
		}
	}
	void performIteration() {

	}
	Matrix& getLHS(Matrix &LHS, Matrix &JJT, Matrix &diagJJT) {
		LHS.copyMatrix(JJT);
		LHS += diagJJT;
		LHS.invert();
		return LHS;
	}
	Matrix& getDiagJJT(Matrix& diagJJT, Matrix& J, Matrix& JT) {
		diagJJT.fromMatrixMultiplication(JT, J);
		diagJJT.diagonalise();
		diagJJT *= _lambda;
		return diagJJT;
	}
	Matrix& getRHS() {

	}
public:
	LevenbergMarquardtAlgorithm(MultiVariableFunction &func):_func(func),_nParams(func.getNParams()) {
		initialise();
	}
	~LevenbergMarquardtAlgorithm() {
		deleteAll();
		_canDelete = false;
	}
	float* run(float** data, float* initGuess, unsigned int nPoints) {
		_data_ = data;
		_initGuess_ = initGuess;
		_nPoints = nPoints;
		runAlgorithm();
		return _currentParams_;
	}

};

#endif
