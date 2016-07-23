#include "MultiVariableFunction.h"
#include "LevenbergMarquardtAlgorithm.h"
#include "Matrix.h"
#include <cmath>
#include <iostream>

float** getData(MultiVariableFunction* f, float* trueParams, unsigned int len);

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

float getRand(float lower, float upper) {
	float randNum = (rand() % 1000) / 1000.f;
	float d = upper - lower;
	return lower + d*randNum;
}

float** getData(MultiVariableFunction* f, float* trueParams, unsigned int len) {
	float** data = new float*[2];
	data[0] = new float[len];
	data[1] = new float[len];
	float start = -10;
	float end = 10;

	for (int i = 0; i < 100; i++) {
		data[0][i] = start + (end - start)*i / 100.f;
		data[1][i] = f->getValue(trueParams, data[0][i]) *getRand(0.97f, 1.03f);
	}
	return data;
}

int main() {
	float trueParams[3] = { 3,4,5 };
	float initGuess[3] = { 5,4,3 };
	unsigned int len = 100;
	MultiVariableFunction* func = new TestFunction();
	LevenbergMarquardtAlgorithm lma(func);
	float** data = getData(func, trueParams, len);
	lma.run(data, initGuess, len);
	//std::cout << "finished in main" << std::endl;
}