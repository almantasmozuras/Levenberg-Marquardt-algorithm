#ifndef H_MULTIVARIABLEFUNCTION_H
#define H_MULTIVARIABLEFUNCTION_H

//class that the user inherits when creating their own function for use in LMA, where f(a,b,c...). getValue returns the
//value of the function f given the current value of the parameters a, b, c...
//getJacobian returns the Jacobian (array of the partial derivates of f with respect to a, b, c, etc. when evaluated
//with the current parameters
class MultiVariableFunction {
private:
	int _nParams;
public:
	MultiVariableFunction(int nParams) :_nParams(nParams) {

	}
	int getNParams() {
		return _nParams;
	}
	virtual float getValue(float* params) {
		return 0;
	}
	virtual float* getJacobian(float* params) {
		return 0;
	}
};

#endif