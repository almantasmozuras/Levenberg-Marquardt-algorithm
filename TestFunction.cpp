#include "MultiVariableFunction.h"
#include "LevenbergMarquardtAlgorithm.h"
#include "Matrix.h"
#include <cmath>
#include <iostream>

double** getData(MultiVariableFunction* f, double* trueParams, unsigned int len);

class TestFunction :public MultiVariableFunction {
public:
	TestFunction():MultiVariableFunction(3) {

	}
	//Function of the form a^2*x^2 +b*x + a*c
	double getValue(double* params, double x) {
		double a = params[0];
		double b = params[1];
		double c = params[2];
		return a*x*x + b*x + c;
	}
	double* getJacobian(double* params, double* returnedArray, double x) {
		double a = params[0];
		double b = params[1];
		double c = params[2];
		returnedArray[0] = x*x;
		returnedArray[1] = x;
		returnedArray[2] = 1;
		return returnedArray;
	}
};

double getRand(double lower, double upper) {
	double randNum = (rand() % 1000) / 1000.f;
	double d = upper - lower;
	return lower + d*randNum;
}

double** getData(MultiVariableFunction* f, double* trueParams, unsigned int len) {
	double** data = new double*[2];
	data[0] = new double[len];
	data[1] = new double[len];
	double start = 0;
	double end = 10;

	for (unsigned int i = 0; i < len; i++) {
		data[0][i] = start + (end - start)*i / 100;
		data[1][i] = f->getValue(trueParams, data[0][i])*0.95;
	}
	return data;
}

int main() {
	double trueParams[3] = { 3,4,5};
	double initGuess[3] = { 1,2,1};
	unsigned int len = 200;
	MultiVariableFunction* func = new TestFunction();
	LevenbergMarquardtAlgorithm lma(func);
	double** data = getData(func, trueParams, len);
	lma.run(data, initGuess, len);
}

/*

*/