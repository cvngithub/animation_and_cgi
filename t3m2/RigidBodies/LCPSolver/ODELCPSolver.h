#ifndef __ODE_LCP_SOLVER_H__
#define __ODE_LCP_SOLVER_H__

#include "FOSSSim/MathUtilities.h"
#include "lcp.h"

namespace lcputils
{

// TODO: IMPORTANT NOTE IF USING SOMETHING OTHER THAN DOUBLES
//  ODE assumes you are using floats or doubles, and also makes some assumptions
//  about the number of bits used to store these values. Be careful if moving
//  to long doubles or floats with this code base.
void solveLCPwithODE( const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& lambda );

}

#endif
