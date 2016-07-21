#ifndef H_TESTFUNCTION_H
#define H_TESTFUNCTION_H
#include "MultiVariableFunction.h"
#include "LevenbergMarquardtAlgorithm.h"

class TestFunction :public MultiVariableFunction {
public:
	TestFunction():MultiVariableFunction(3) {

	}
	//Function of the form a^2*x^2 +b*x + a*c
	float getValue(float* params, float &x) {
		float a = params[0];
		float b = params[1];
		float c = params[2];
		return a*a*x*x + b*x + a*c;
	}
	float* getJacobian(float* params, float* returnedArray, float &x) {
		float a = params[0];
		float b = params[1];
		float c = params[2];
		returnedArray[0] = 2*a*x*x + c;
		returnedArray[1] = x;
		returnedArray[2] = a;
		//std::cout << returnedArray[0] << "; " << returnedArray[1] << "; " << returnedArray[2] << std::endl;
		return returnedArray;
	}
};

#endif