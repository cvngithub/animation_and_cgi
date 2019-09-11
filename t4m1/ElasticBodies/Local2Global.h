#ifndef __LOCAL2GLOBAL_H__
#define __LOCAL2GLOBAL_H__

#include "../MathDefs.h"
#include <Eigen/Core>
#include <iostream>

// Transform global Jacobian based on local one
void Local2Global(const std::vector<int> &i_points, const MatrixXs &hessE_local, MatrixXs &hessE_global);

#endif
