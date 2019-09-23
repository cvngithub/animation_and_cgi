#include "SpringForce.h"
#include <iostream>

void SpringForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Add milestone 2 code here.
  const Vector2s& ri = x.segment<2>(2*m_endpoints.first);
  const Vector2s& rj = x.segment<2>(2*m_endpoints.second);
  const scalar l = (ri - rj).norm();
  
  E += 0.5 * m_k * (l - m_l0) * (l - m_l0);
}

void SpringForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Add milestone 2 code here.
  // spring elastic force
  const Vector2s& ri = x.segment<2>(2*m_endpoints.first);
  const Vector2s& rj = x.segment<2>(2*m_endpoints.second);
  const scalar l = (ri - rj).norm();
  const Vector2s n = (ri - rj) / l;
  Vector2s toadd = +m_k * ( l - m_l0 ) * n;

  // spring damping force
  const Vector2s& vi = v.segment<2>(2*m_endpoints.first);
  const Vector2s& vj = v.segment<2>(2*m_endpoints.second);
  toadd += m_b * n.dot(vi - vj) * n; 

  // add to the accumulator
  gradE.segment<2>(2*m_endpoints.first)  += toadd;
  gradE.segment<2>(2*m_endpoints.second) += -toadd;

}
