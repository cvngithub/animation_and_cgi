#include "RigidBody.h"

Vector2s RigidBody::computeTotalMomentum() const
{
  // COMPLETE THIS CODE
  return m_V * m_M;
}

scalar RigidBody::computeCenterOfMassAngularMomentum() const
{
  // COMPLETE THIS CODE
  Vector2s momentum = computeTotalMomentum();
  return m_X(0) * momentum(1) - m_X(1) * momentum(0);
}

scalar RigidBody::computeSpinAngularMomentum() const
{
  // COMPLETE THIS CODE
  return m_omega * m_I;
}

scalar RigidBody::computeCenterOfMassKineticEnergy() const
{
  // COMPLETE THIS CODE
  return 0.5*m_M*m_V.squaredNorm();
}

scalar RigidBody::computeSpinKineticEnergy() const
{
  // COMPLETE THIS CODE
  return 0.5*m_I*m_omega*m_omega;
}

scalar RigidBody::computeTotalMass( const VectorXs& masses ) const
{
  // COMPLETE THIS CODE
  scalar totalMass = 0;
  for(int i=0;i<masses.size();i++) {
      totalMass += masses(i);
  }
  return totalMass;}

Vector2s RigidBody::computeCenterOfMass( const VectorXs& vertices, const VectorXs& masses ) const
{
  // COMPLETE THIS CODE
  Vector2s com;
  com.setZero();
  for(int i=0;i<masses.size();i++) {
      com+=masses(i)*vertices.segment<2>(2*i);
  }
  com /= m_M;
  return com;
}

scalar RigidBody::computeMomentOfInertia( const VectorXs& vertices, const VectorXs& masses ) const
{
  assert( vertices.size()%2 == 0 );
  assert( 2*masses.size() == vertices.size() );
  
  // COMPLETE THIS CODE
  scalar I = 0;
  for(int i=0;i<masses.size();i++) {
      I += masses(i)*(m_X-vertices.segment<2>(2*i)).squaredNorm();
  }
  return I;}
