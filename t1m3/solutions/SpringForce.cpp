#include "SpringForce.h"

void SpringForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Implement force Jacobian here!
  
  // Contribution from elastic component

  // Contribution from damping
  Vector2s pos_1 = x.segment<2>(2 * m_endpoints.first);
  Vector2s pos_2 = x.segment<2>(2 * m_endpoints.second);

  Vector2s vel_1 = v.segment<2>(2 * m_endpoints.first);
  Vector2s vel_2 = v.segment<2>(2 * m_endpoints.second);

  

  Vector2s n_vector = pos_2 - pos_1;
  Vector2s v_vector = vel_2 - vel_1;
  Vector2s n_hat_vector = n_vector / n_vector.norm();

  Matrix2s Identity = Matrix2s::Identity();
  
  Matrix2s K_matrix = -m_k * ((n_hat_vector * n_hat_vector.transpose()) + ((n_vector.norm() - m_l0) / n_vector.norm())*(Identity - (n_hat_vector * n_hat_vector.transpose()))); 

  Matrix2s K_damp_matrix = (-m_b / n_vector.norm()) * (n_hat_vector.dot(v_vector)* Identity + n_hat_vector * (v_vector.transpose())) * ( Identity - (n_hat_vector * n_hat_vector.transpose()));

  hessE.block<2,2>(2 * m_endpoints.first, 2 * m_endpoints.first) += K_matrix + K_damp_matrix;
  hessE.block<2,2>(2 * m_endpoints.first, 2 * m_endpoints.second) -= K_matrix + K_damp_matrix;
  hessE.block<2,2>(2 * m_endpoints.second, 2 * m_endpoints.first) -= K_matrix + K_damp_matrix;
  hessE.block<2,2>(2 * m_endpoints.second, 2 * m_endpoints.second) += K_matrix + K_damp_matrix;

}

void SpringForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Implement force Jacobian here!

  // Contribution from damping
  Vector2s pos_1 = x.segment<2>(2 * m_endpoints.first);
  Vector2s pos_2 = x.segment<2>(2 * m_endpoints.second);
  Vector2s n_vector = pos_2 - pos_1;
  Vector2s n_hat_vector = n_vector / n_vector.norm();

  Matrix2s Beta_matrix  = m_b * (n_hat_vector * n_hat_vector.transpose());

  hessE.block<2,2>(2 * m_endpoints.first, 2 * m_endpoints.first) -= Beta_matrix;
  hessE.block<2,2>(2 * m_endpoints.first, 2 * m_endpoints.second) += Beta_matrix;
  hessE.block<2,2>(2 * m_endpoints.second, 2 * m_endpoints.first) += Beta_matrix;
  hessE.block<2,2>(2 * m_endpoints.second, 2 * m_endpoints.second) -= Beta_matrix;
}
