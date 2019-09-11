#include "RigidBodyLCPCollisionResolver.h"
#include <iostream>

using namespace std;

void removeRows(MatrixXs& matrix, int r, int num) {
    int n = matrix.rows();
    int m = matrix.cols();
    // matrix.block(r,c,p,q)  means starting at (r,c), block of size (p,q)
    if (r+num < n)
        matrix.block(r,0,n-num-r,m) = matrix.block(r+num,0,n-num-r,m);
    matrix.conservativeResize(n-num,m);
}

void RigidBodyLCPCollisionResolver::resolveCollisions( std::vector<RigidBody>& rbs, const std::set<RigidBodyCollision>& rbcs )
{
  int Nrb = rbs.size();  // number of rigid bodies in the system
  int F = 0;  // count how many rigid bodies are fixed
  for (int rb=0; rb<Nrb; rb++) {
      if (rbs[rb].isFixed()) {
          F++;
      }
  }
  int Nfree = Nrb-F;  // number of non-fixed rigid bodies
  MatrixXs M = MatrixXs::Zero(3*Nfree,3*Nfree);  // M is sized according to non-fixed rigid bodies
  MatrixXs Gradr_i(2,3);
  MatrixXs Gradr_j(2,3);
  VectorXs qdot(3*Nfree);  // qdot is a vector of length 3n
  int rb_free = 0;  // index for free rigid bodies
  for (int rb=0; rb<Nrb; rb++) {
      if (rbs[rb].isFixed()) {
              continue;
      }
      else {
          M(3*rb_free,3*rb_free) = rbs[rb].getM();  // get mass of this rigid body
          M(3*rb_free+1,3*rb_free+1) = rbs[rb].getM();  // get mass of this rigid body
          M(3*rb_free+2,3*rb_free+2) = rbs[rb].getI();  // get moment of inertia for this rigid body
          Vector2s v = rbs[rb].getV();  // this is the velocity of the center of mass
          qdot(3*rb_free) = v[0];
          qdot(3*rb_free+1) = v[1];
          qdot(3*rb_free+2) = rbs[rb].getOmega();  // omega of rigid body
          rb_free++;
      }
  }
    
    
  MatrixXs Minv = M.inverse();
  //cout << "Minv: " << endl << Minv << endl;
  
  // Compute N matrix, A, b
  int Ncoll = rbcs.size();  // number of collisions
  //cout << "Ncoll: " << Ncoll << endl;
  if (Ncoll > 0) {
      int k=0;  // collision iterator
      int i;
      int j;
      Vector2s ri;
      Vector2s rj;
      Vector2s nhat;
      MatrixXs N(3*Nrb,Ncoll);
      VectorXs cn(3*Nrb);
      MatrixXs Gamma(2,3*Nrb);
      for (set<RigidBodyCollision>::iterator it = rbcs.begin(); it != rbcs.end(); ++it) {  // for each collision i, j
          const RigidBodyCollision rbc = *it; // Note the "*" here. rbc = set{i,j,r0,r1,nhat}
          i = rbc.i0;  // index to first rigid body
          j = rbc.i1;  // index to second rigid body
          ri = rbc.r0;  // vector from X[i] to ri point
          rj = rbc.r1;  // vector from X[j] to rj point
          nhat = rbc.nhat;  //normal from i to j
          //cout << "rcbs retrieved: ";
          //cout << "i=" << i << " j=" << j << " ";
          //cout << "ri " << ri[0] << " " << ri[1] << " ";
          //cout << "rj " << rj[0] << " " << rj[1] << " ";
          //cout << nhat[0] << " " << nhat[1] << endl;
          scalar theta_i = rbs[i].getTheta();
          scalar theta_j = rbs[j].getTheta();
          //Gradr_i << 1.0, 0.0, -sin(theta_i)*ri[0]-cos(theta_i)*ri[1],
          //           0.0, 1.0, cos(theta_i)*ri[0]-sin(theta_i)*ri[1];  // ri is in body space
          //Gradr_j << 1.0, 0.0, -sin(theta_j)*rj[0]-cos(theta_j)*rj[1],
          //           0.0, 1.0, cos(theta_j)*rj[0]-sin(theta_j)*rj[1];  // rj is in body space
          Gradr_i << 1.0, 0.0, -ri[1],
                     0.0, 1.0, ri[0];  // ri is in world space
          Gradr_j << 1.0, 0.0, -rj[1],
                     0.0, 1.0, rj[0];  // rj is in world space
          Gamma = MatrixXs::Zero(2,3*Nrb);
          Gamma.block<2,3>(0,3*i) = Gradr_i;  // stuff in Gradr_i
          Gamma.block<2,3>(0,3*j) = -Gradr_j;   // stuff in Gradr_j
          //cout << "Gamma " << "i=" << i << " j=" << j << endl;
          //cout << Gamma << endl;
          //cout << "nhat:" << endl;
          //cout << nhat[0] << " " << nhat[1] << endl;
          cn = Gamma.transpose() * (-nhat);  // Size = 3*Nrb
          //cout << cn << endl;
          N.col(k) = cn;  // stuff cn into column k
          k++;  // increment iterator
      }
      // Downsize N
      //cout << "sizeN: " << int(N.size()) << endl;
      //cout << Nfree << " " << Nrb << endl;
      if (Nfree < Nrb) {
          for (int rb=0; rb<Nrb; rb++) {
              if (rbs[rb].isFixed()) {
                  removeRows(N, 3*rb, 3);  // remove rows 3*rb, 3*rb+1, 3*rb+2
              }
          }
      }
      //
      //cout << "N: " << endl << N << endl;
      
      
      // Compute A matrix, and b vector
      MatrixXs NT = N.transpose();
      MatrixXs A = NT * Minv * N;  // Compute matrix A  size should be Ncoll x Ncoll
      VectorXs b = NT * qdot;  // Compute vector b  size should be Ncoll
      //cout << "A: " << endl << A << endl;
      //cout << "b:" << endl << b << endl;
      //cout << N.transpose() << endl;
      //cout << qdot << endl;
      //
      VectorXs lambda(int(b.size()));  // size should be #collisions
      lcputils::solveLCPwithODE( A, b, lambda );  // Now optimize A*lambda + b --> therefore the size of lambda is num_coll
      //cout << "lambda: " << endl << lambda << endl;
      // TO complete this method: modify V and Omega of the rigid bodies to prevent penetration
      VectorXs dV = Minv * N * lambda;
      //cout << "dV: " << endl << dV << endl;
      //cout << "qdot: " << endl << qdot << endl;
      //
      rb_free = 0;  // index to dV for free rigid bodies
      for (int rb=0; rb<Nrb; rb++) {
          if (rbs[rb].isFixed()) {
              continue;
          }
          else {
              rbs[rb].getV() += dV.segment<2>(3*rb_free);    //Block of n elements, starting at i: vector.segment<n>(i);
              rbs[rb].getOmega() += dV(3*rb_free+2);
              rb_free++;
          }
      }
      //
  }
}
  // Example of using the lcp solver
  /*
  MatrixXs A(4,4);
  A << 0.828869, 0.337798, -0.28125, -0.21875,
       0.337798, 0.828869, -0.21875, -0.28125,
       -0.28125, -0.21875, 0.828869, 0.337798,
       -0.21875, -0.28125, 0.337798, 0.828869;

  VectorXs b(4);
  b << -1, -1, -1, -1;

  VectorXs lambda(4);
  lcputils::solveLCPwithODE( A, b, lambda );
  */

