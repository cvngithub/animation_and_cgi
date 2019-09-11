#include "RigidBodyWindForce.h"

void RigidBodyWindForce::computeForceAndTorque( std::vector<RigidBody>& rbs )
{
  // COMPLETE THIS CODE
  // for all rigid bodies i rbs[i].getForce()  += ... some force you compute ...
  //                        rbs[i].getTorque() += ... some torque you compute ...
  for (std::vector<RigidBody>::size_type i = 0; i < rbs.size(); i++) {
      for (int j = 0; j < rbs[i].getNumEdges(); j++) {
          Vector2s e = rbs[i].computeWorldSpaceEdge(j);
          Vector2s v = rbs[i].getWorldSpaceVertex(j);
          Vector2s nhat;
          nhat(0) = -e(1);
          nhat(1) = e(0);
          nhat /= nhat.norm();
          for (int k = 0; k < m_num_quadrature_points; k++) {
              Vector2s l = e / (2 * m_num_quadrature_points);
              Vector2s samplePoint = v + l * (2 * k + 1);
              Vector2s f = m_beta * e.norm() / m_num_quadrature_points *
              (m_wind - rbs[i].computeWorldSpaceVelocity(samplePoint)).dot(nhat) * nhat;
              rbs[i].getForce() += f;
              Vector2s x = samplePoint - rbs[i].getX();
              rbs[i].getTorque() += x(0) * f(1) - x(1) * f(0);
          }
      }
  }
}
