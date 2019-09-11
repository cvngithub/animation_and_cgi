#include "LinearizedImplicitEuler.h"

bool LinearizedImplicitEuler::stepScene( TwoDScene& scene, scalar dt )
{
  VectorXs& x = scene.getX();
  VectorXs& v = scene.getV();
  const VectorXs& m = scene.getM();
  assert(x.size() == v.size());
  assert(x.size() == m.size());

  // Implement implicit euler here!
  int ndof = x.size();
  assert( ndof%2 == 0 );
  MatrixXs temp = MatrixXs::Zero(ndof,ndof);
  VectorXs temp1 = VectorXs::Zero(ndof);
  MatrixXs Mass_m = temp;
  VectorXs dv = temp1;

  for (int iter = 0; iter < ndof; ++iter)
  {
    Mass_m(iter, iter) = m[iter];
  }

  
  MatrixXs Force_pos_diff = temp;
  MatrixXs Force_vel_diff = temp;

  scene.accumulateddUdxdx(Force_pos_diff, v*dt, dv);
  scene.accumulateddUdxdv(Force_vel_diff, v*dt, dv);

  MatrixXs A = Mass_m - dt * dt * Force_pos_diff - Force_vel_diff *dt;

  VectorXs Force = temp1;
  scene.accumulateGradU(Force, v*dt, dv);
  VectorXs B = - Force * dt;

  for (int iter = 0; iter < scene.getNumParticles(); ++iter)
    {
      if (scene.isFixed(iter)) {
        B.segment<2>(2 * iter).setZero();
        A.row(2 * iter).setZero();
        A.row(2 * iter + 1).setZero();
        A.col(2 * iter).setZero();
        A.col(2 * iter + 1).setZero();
        A(2 * iter, 2 * iter) = A(2 * iter + 1, 2 * iter + 1) = 1;
      }
    }


  VectorXs sln = A.fullPivLu().solve(B);
  v = v + sln;
  x = x + dt*v;
  
  return true;
}
