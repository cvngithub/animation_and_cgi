#include "SymplecticEuler.h"

bool SymplecticEuler::stepScene( TwoDScene& scene, scalar dt )
{
  /* Add milestone 2 code here.      */
  VectorXs& x = scene.getX();
  VectorXs& v = scene.getV();
  VectorXs& m = scene.getM();
  VectorXs gradU = VectorXs::Zero(x.size());
  // get all the forces
  scene.accumulateGradU( gradU );

  VectorXs fixed = VectorXs::Zero(m.size());
  for (int i = 0; i < m.size()/2; i++)
  {
    fixed[2*i+0] = !scene.isFixed(i);
    fixed[2*i+1] = !scene.isFixed(i);
  }

  MatrixXs inverse_m = m.asDiagonal().inverse();
  v += fixed.cwiseProduct(dt * ( -inverse_m * gradU));
  x += fixed.cwiseProduct(dt * v);

  return true;
}






