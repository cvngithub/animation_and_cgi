#include "GravitationalForce.h"

void GravitationalForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );

  // Compute the force Jacboian here!
  Vector2s pos_1 = x.segment<2>(2 * m_particles.first);
  Vector2s pos_2 = x.segment<2>(2 * m_particles.second);

  Vector2s n_vector = pos_2 - pos_1;
  Vector2s n_hat_vector = n_vector / n_vector.norm();

  Matrix2s Identity = Matrix2s::Identity();
  Matrix2s K_matrix = - (m_G * m(2* m_particles.first) * m(2* m_particles.second) / pow(n_vector.norm(),3)) * (Identity - 3 *(n_hat_vector * n_hat_vector.transpose()));

  hessE.block<2,2>(2 * m_particles.first, 2 * m_particles.first) +=  K_matrix;
  hessE.block<2,2>(2 * m_particles.first, 2 * m_particles.second) -=  K_matrix;
  hessE.block<2,2>(2 * m_particles.second, 2 * m_particles.first) -=  K_matrix;
  hessE.block<2,2>(2 * m_particles.second, 2 * m_particles.second) +=  K_matrix;

}

void GravitationalForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );
  // Nothing to do.
}
