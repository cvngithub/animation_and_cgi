#include "ImplicitEuler.h"

bool ImplicitEuler::stepScene( TwoDScene& scene, scalar dt )
{
  VectorXs& x = scene.getX();
  VectorXs& v = scene.getV();
  const VectorXs& m = scene.getM();
  assert(x.size() == v.size());
  assert(x.size() == m.size());

  // Implement implicit euler here for extra credit!

  
  
  /////////////
  // Some examples of working with vectors, matrices, etc.
  
  // How to get the force Jacobian from two d scene
  int ndof = x.size();
  assert( ndof%2 == 0 );
  // Note that the system's state is passed to two d scene as a change from the last timestep's solution
  VectorXs dx = dt*v;
  VectorXs temp1 = VectorXs::Zero(ndof);
  MatrixXs temp = MatrixXs::Zero(ndof,ndof);
  VectorXs dv = temp1;
  MatrixXs A = temp;
  VectorXs B = temp1;
  MatrixXs Mass_m = temp;

  scalar num = std::numeric_limits<scalar>::infinity();

  for (int iter = 0; iter < ndof; ++iter)
  {
    Mass_m(iter, iter) = m[iter];
  }

  while (num >= 1.0e-9)
  {
    VectorXs Force = temp1;
    scene.accumulateGradU(Force, dx, dv);
    B = - Mass_m * dv - Force * dt;

    MatrixXs Force_pos_diff = temp;
    MatrixXs Force_vel_diff = temp;
    scene.accumulateddUdxdx(Force_pos_diff, dx, dv);
    scene.accumulateddUdxdv(Force_vel_diff, dx, dv);
    A = Mass_m - dt * dt * Force_pos_diff - Force_vel_diff *dt;


    for (int iter = 0; iter < scene.getNumParticles(); ++iter)
    {
      if (scene.isFixed(iter)) 
      {
        B.segment<2>(2 * iter).setZero();
        A.row(2 * iter).setZero();
        A.row(2 * iter + 1).setZero();
        A.col(2 * iter).setZero();
        A.col(2 * iter + 1).setZero();
        A(2 * iter, 2 * iter) = A(2 * iter + 1, 2 * iter + 1) = 1;
      }
    }

    VectorXs sln = A.fullPivLu().solve(B);
    dv += sln;
    dx = (v + dv) * dt;
    num = sln.norm();
  }

  v = v + dv;
  x = x + dt*v;
  
  return true;
}
