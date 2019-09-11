#include "RigidBodyGravityForce.h"

scalar RigidBodyGravityForce::computePotentialEnergy( const std::vector<RigidBody>& rbs )
{
  // COMPLETE THIS CODE
  scalar p = 0;
    for (std::vector<RigidBody>::size_type i = 0; i < rbs.size(); i++) {
        p += -rbs[i].getM() * m_g.dot(rbs[i].getX());
    }
    return p;
}

void RigidBodyGravityForce::computeForceAndTorque( std::vector<RigidBody>& rbs )
{
  // COMPLETE THIS CODE
  // for all rigid bodies i rbs[i].getForce()  += ... some force you compute ...
  //                        rbs[i].getTorque() += ... some torque you compute ...
  for (std::vector<RigidBody>::size_type i = 0; i < rbs.size(); i++) {
        Vector2s f = rbs[i].getM() * m_g;
        rbs[i].getForce() += f;
    }
}
