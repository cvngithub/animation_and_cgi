#include "RigidBodySpringForce.h"

scalar RigidBodySpringForce::computePotentialEnergy( const std::vector<RigidBody>& rbs )
{
  assert( m_rb0 >= -1 ); assert( m_rb0 < (int) rbs.size() );
  assert( m_rb1 >= -1 ); assert( m_rb1 < (int) rbs.size() );
  
    scalar p = 0;
    scalar l = (computeFirstEndpoint(rbs) - computeSecondEndpoint(rbs)).norm();
    p += 0.5 * m_k * (l - m_l0) * (l - m_l0);
    return p;
}

void RigidBodySpringForce::computeForceAndTorque( std::vector<RigidBody>& rbs )
{
  assert( m_rb0 >= -1 ); assert( m_rb0 < (int) rbs.size() );
  assert( m_rb1 >= -1 ); assert( m_rb1 < (int) rbs.size() );
  
  // COMPLETE THIS CODE
  // for all rigid bodies i rbs[i].getForce()  += ... some force you compute ...
  //                        rbs[i].getTorque() += ... some torque you compute ...
  Vector2s l = computeFirstEndpoint(rbs) - computeSecondEndpoint(rbs);
  if (l.norm() == 0 && m_l0 == 0) return;
  if (m_rb0 != -1) {
      Vector2s f1 = -m_k * (l.norm() - m_l0) * l / l.norm();
      Vector2s x = computeFirstEndpoint(rbs) - rbs[m_rb0].getX();
      rbs[m_rb0].getForce() += f1;
      rbs[m_rb0].getTorque() += x(0) * f1(1) - x(1) * f1(0);
  }
  if (m_rb1 != -1) {
      Vector2s f2 = m_k * (l.norm() - m_l0) * l / l.norm();
      Vector2s x = computeSecondEndpoint(rbs) - rbs[m_rb1].getX();
      rbs[m_rb1].getForce() += f2;
      rbs[m_rb1].getTorque() += x(0) * f2(1) - x(1) * f2(0);
  }
}
