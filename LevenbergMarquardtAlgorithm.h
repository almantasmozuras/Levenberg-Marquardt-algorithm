#ifndef H_LMA_H
#define H_LMA_H
#include "MultiVariableFunction.h"
class LevenbergMarquardtAlgorithm {
private:
	MultiVariableFunction _func;
public:
	LevenbergMarquardtAlgorithm(MultiVariableFunction func):_func(func) {

	}
};

#endif
