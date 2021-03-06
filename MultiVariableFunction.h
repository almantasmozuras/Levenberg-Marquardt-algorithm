#ifndef H_MULTIVARIABLEFUNCTION_H
#define H_MULTIVARIABLEFUNCTION_H

//class that the user inherits when creating their own function for use in LMA, where f(x,a,b,c...). getValue returns the
//value of the function f given the current value of the parameters a, b, c... at the position x
//getJacobian returns the Jacobian (array of the partial derivates of f with respect to a, b, c, etc. when evaluated
//with the current parameters
class MultiVariableFunction {
protected:
	unsigned int _nParams;
public:
	MultiVariableFunction(unsigned int nParams) :_nParams(nParams) {

	}
	unsigned int getNParams() {
		return _nParams;
	}
	virtual double getValue(double* params, double x) = 0;
	virtual double* getJacobian(double* params, double* returnedArray, double x) = 0;
};

#endif