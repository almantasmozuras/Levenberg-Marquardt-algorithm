#include "stdafx.h"

#include "MultiVariableFunction.h"
#include "LevenbergMarquardtAlgorithm.h"
#include "Matrix.h"
#include <cmath>
#include <iostream>

void getData(MultiVariableFunction* f, double** data, double* trueParams, unsigned int nPoints, double noise);
double getRand(double lower, double upper);

class TestFunction :public MultiVariableFunction {
public:
	TestFunction() :MultiVariableFunction(3) {

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

void getData(MultiVariableFunction* f, double** data, double* trueParams, unsigned int nPoints, double noise) {
	//double** data = new double*[2];
	//data[0] = new double[nPoints];
	//data[1] = new double[nPoints];
	double start = 0;
	double end = 10;

	for (unsigned int i = 0; i < nPoints; i++) {
		data[0][i] = start + (end - start)*i / 100;
		data[1][i] = f->getValue(trueParams, data[0][i]) + getRand(-noise, noise);//*0.95;;//
	}
	//return data;
}

int main()
{
	double noise = 10;//0.1;
	std::cout << "noise =	" << noise << std::endl;
	double trueParams[3] = { 3,4,5 };//{ 10,8,9 };

	double initGuess[3] = { 1,2,1 };//{ 0,0,0 }; 
	unsigned int nPoints = 200;
	MultiVariableFunction* func = new TestFunction();

	LevenbergMarquardtAlgorithm lma(func);
	double** data = new double*[2];
	data[0] = new double[nPoints];
	data[1] = new double[nPoints];
	getData(func, data, trueParams, nPoints, noise);
	std::cout <<"x" << "               " << "y" << std::endl;
	for (unsigned int j = 0; j < nPoints; j++)
	{
		std::cout << data[0][j]<<"               "<< data[1][j] << std::endl;	
	}



	lma.run(data, initGuess, nPoints);

	std::cout << "true params =  " << trueParams[0]<< "   " << trueParams[1] << "   "<< trueParams[2] << std::endl;

	for (unsigned int i = 0; i < 2; i++) {
		delete[] data[i];
	}
	delete[] data;
	
	return 0;
}
