#include "RigidBodySymplecticEuler.h"

bool RigidBodySymplecticEuler::stepScene( RigidBodyScene& scene, scalar dt )
{
  // Clear any previously computed forces in the rigid bodies
  std::vector<RigidBody>& rbs = scene.getRigidBodies();
  for( std::vector<RigidBody>::size_type i = 0; i < rbs.size(); ++i ) rbs[i].getForce().setZero();
  for( std::vector<RigidBody>::size_type i = 0; i < rbs.size(); ++i ) rbs[i].getTorque() = 0.0;
  
  // Add each force's contribution to each rigid body using previous step's state
  std::vector<RigidBodyForce*>& frcs = scene.getForces();
  for( std::vector<RigidBodyForce*>::size_type i = 0; i < frcs.size(); ++i ) frcs[i]->computeForceAndTorque(rbs);
  
  // Useful method:
  //   RigidBody::getX()
  //   RigidBody::getV()
  //   RigidBody::getTheta()
  //   RigidBody::getOmega()
  //   RigidBody::getForce()
  //   RigidBody::getTorque()
  
  // For each rigid body
  for( std::vector<RigidBody>::size_type i = 0; i < rbs.size(); ++i )
  {
    // COMPLETE THIS CODE
    Vector2s f = rbs[i].getForce();
    scalar m = rbs[i].getM();
    rbs[i].getV() += dt * f / m;
    rbs[i].getX() += rbs[i].getV() * dt;
    if (rbs[i].getI() != 0)
        rbs[i].getOmega() += dt * rbs[i].getTorque() / rbs[i].getI();
    rbs[i].getTheta() += rbs[i].getOmega() * dt;
  }
  
  for( std::vector<RigidBodyForce*>::size_type i = 0; i < rbs.size(); ++i ) rbs[i].updateDerivedQuantities();
  
  return true;
}
