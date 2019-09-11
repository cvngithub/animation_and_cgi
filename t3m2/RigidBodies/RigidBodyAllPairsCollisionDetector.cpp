#include "RigidBodyAllPairsCollisionDetector.h"
#include <iostream>

using namespace std;

void RigidBodyAllPairsCollisionDetector::detectCollisions( const std::vector<RigidBody>& rbs, std::set<RigidBodyCollision>& collisions )
{
  // Compute all vertex-edge collisions between all pairs of rigid bodies
  int N = rbs.size();  // number of rigid bodies in the system
  for (int i=0; i<N; i++) {  // for each rigid body
      scalar r_vertex = rbs[i].getRadius();
      int Nvi = rbs[i].getNumVertices();  // same as number of edges
      for (int v=0; v<Nvi; v++) {  // for each vertex in rigid body i0
          Vector2s x1 = rbs[i].getWorldSpaceVertex(v);
          for (int j=0; j<N; j++) {  // for every other rigid body i1
              if (i==j) continue;  // ignore same body
              int Nvj = rbs[j].getNumVertices();  // same as number of edges
              scalar r_edge = rbs[j].getRadius();
              for (int e=0; e<Nvj; e++) {  // for each edge in rigid body i1
                  //Vector2s Edge = rbs[j].computeWorldSpaceEdge(e);
                  Vector2s x2 = rbs[j].getWorldSpaceVertex(e);
                  Vector2s x3 = rbs[j].getWorldSpaceVertex((e+1)%Nvj);
                  scalar alpha = (x1 -x2).dot(x3-x2)/(x3-x2).squaredNorm();
                  alpha = min(max(alpha, 0.0), 1.0);
                  Vector2s x_alpha = x2 + alpha*(x3-x2);
                  Vector2s n = x_alpha - x1;  // let the contact normal be from i to j
                  if ((n.squaredNorm()+0e-6) < ((r_vertex+r_edge)*(r_vertex+r_edge))) {
                      Vector2s nhat = n.normalized(); //
                      Vector2s Xc2 = x_alpha - nhat*r_edge; // point of contact
                      Vector2s Xc = x1 + nhat*r_vertex; // point of contact
                      //cout << "contact error: " << (Xc2-Xc).norm() << endl;
                      //cout << "point of contact " << Xc << endl;
                      Vector2s ri = Xc - rbs[i].getX(); // vector from center of mass to point of contact (world space)
                      //Vector2s rj = Xc2 - rbs[j].getX();// vector from center of mass to point of contact (world space)
                      Vector2s rj = Xc - rbs[j].getX();// vector from center of mass to point of contact (world space)
                      //cout << "collision detected: ";
                      //cout << "i=" << i << " j=" << j << " ";
                      //cout << "ri " << ri[0] << " " << ri[1] << " ";
                      //cout << "rj " << rj[0] << " " << rj[1] << " ";
                      //cout << nhat[0] << " " << nhat[1] << endl;
                      RigidBodyCollision rbc(i,j,ri,rj,nhat);
                      collisions.insert(rbc);  //add to set collisions
                  }
              }
          }
      }
  }
}
