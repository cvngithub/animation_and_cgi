#include "GravitationalForce.h"

void GravitationalForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%2 == 0 );
  assert( m_particles.first >= 0 );  assert( m_particles.first < x.size()/2 );
  assert( m_particles.second >= 0 ); assert( m_particles.second < x.size()/2 );

  // Add milestone 2 code here.
  const Vector2s& ri = x.segment<2>(2*m_particles.first);
  const Vector2s& rj = x.segment<2>(2*m_particles.second);
  const scalar& mi = m[2*m_particles.first];  
  const scalar& mj = m[2*m_particles.second]; 
  const scalar l = (ri - rj).norm();

  E += -m_G * mi * mj / l;

}

void GravitationalForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  assert( m_particles.first >= 0 );  assert( m_particles.first < x.size()/2 );
  assert( m_particles.second >= 0 ); assert( m_particles.second < x.size()/2 );

  // Add milestone 2 code here.
  const Vector2s& ri = x.segment<2>(2*m_particles.first);
  const Vector2s& rj = x.segment<2>(2*m_particles.second);
  const scalar& mi = m[2*m_particles.first];  
  const scalar& mj = m[2*m_particles.second]; 
  const scalar l = (ri - rj).norm();
  const Vector2s n = (ri - rj) / l;

  scalar l2 = l * l;

  Vector2s toadd = m_G * mi * mj * n / l2;

  gradE.segment<2>(2*m_particles.first)  += toadd;
  gradE.segment<2>(2*m_particles.second) += -toadd;
  

}
