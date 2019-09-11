#include "RigidBodyVelocityProjectionCollisionResolver.h"
#include <iostream>

using namespace std;

void rmRows(MatrixXs& matrix, int r, int num) {
    int n = matrix.rows();
    int m = matrix.cols();
    // matrix.block(r,c,p,q)  means starting at (r,c), block of size (p,q)
    if (r+num < n)
        matrix.block(r,0,n-num-r,m) = matrix.block(r+num,0,n-num-r,m);
    matrix.conservativeResize(n-num,m);
}

void RigidBodyVelocityProjectionCollisionResolver::resolveCollisions( std::vector<RigidBody>& rbs, const std::set<RigidBodyCollision>& rbcs )
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
    
  // Compute N matrix, A, b
  int Ncoll = rbcs.size();  // number of collisions
    
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
          scalar theta_i = rbs[i].getTheta();
          scalar theta_j = rbs[j].getTheta();
          Gradr_i << 1.0, 0.0, -ri[1],
                     0.0, 1.0, ri[0];  // ri is in world space
          Gradr_j << 1.0, 0.0, -rj[1],
                     0.0, 1.0, rj[0];  // rj is in world space
          Gamma = MatrixXs::Zero(2,3*Nrb);
          Gamma.block<2,3>(0,3*i) = Gradr_i;  // stuff in Gradr_i
          Gamma.block<2,3>(0,3*j) = -Gradr_j;   // stuff in Gradr_j
          cn = Gamma.transpose() * (-nhat);  // Size = 3*Nrb
          N.col(k) = cn;  // stuff cn into column k
          k++;  // increment iterator
      }
      if (Nfree < Nrb) {  // Downsize N
          for (int rb=0; rb<Nrb; rb++) {
              if (rbs[rb].isFixed()) {
                  rmRows(N, 3*rb, 3);  // remove rows 3*rb, 3*rb+1, 3*rb+2
              }
          }
      }
      
      VectorXs MinusMqdot = -M * qdot;
      //cout << "M: " << endl << M << endl;
      //cout << "-Mqdot: " << endl << MinusMqdot << endl;
      //cout << "N: " << endl << N << endl;
      //cout << qdot << endl;
      
      // Matrix in quadratic form of objective
      QuadProgPP::Matrix<scalar> G;
      // Equality constraints
      QuadProgPP::Matrix<scalar> CE;
      // Inequality constraints
      QuadProgPP::Matrix<scalar> CI;
  
      QuadProgPP::Vector<scalar> g0;
      QuadProgPP::Vector<scalar> ce0;
      QuadProgPP::Vector<scalar> ci0;
      QuadProgPP::Vector<scalar> x;

      // insert values from M matrix here  OBVIOUSLY this has to scale with size M
      int Mrows = M.rows();
      int Mcols = M.cols();
      G.resize(Mrows,Mcols);
      for(int i=0; i< Mrows; i++) for(int j=0; j< Mcols; j++) G[i][j]=0;
      for(int i=0; i< Mcols; i++) G[i][i] = M(i,i);

      // Insert values from -M*qdot here
      int Nmmqd = MinusMqdot.size();
      g0.resize(Nmmqd);
      for(int i=0; i< Nmmqd; i++) g0[i] = MinusMqdot[i];
      
      // No equality constraints, currently
      CE.resize(Mrows,0);  // Check this size !!!!!!!!!!
      ce0.resize(0);

      // Compute the number of inequality constraints
      int Nrows = N.rows();
      int Ncols = N.cols();
      CI.resize(Nrows,Ncols);
      for( int i = 0; i < Nrows; ++i ) for( int j = 0; j < Ncols; ++j ) CI[i][j] = N(i,j);

      // Constant term added to inequality constraints
      ci0.resize(Ncols);
      for(int i=0; i< Ncols; i++) ci0[i] = 0;

      // Solution
      x.resize(Nrows);
      solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
      
      rb_free = 0;  // index to dV for free rigid bodies
      for (int rb=0; rb<Nrb; rb++) {
          if (rbs[rb].isFixed()) {
              continue;
          }
          else {
              Vector2s v(x[3*rb_free],x[3*rb_free+1]);
              scalar omega = x[3*rb_free+2];
              rbs[rb].getV() = v;
              rbs[rb].getOmega() = omega;
              rb_free++;
          }
      }
  }
}

  /*
  // Example of using QuadProg++
  // For detailed documentation, please see FOSSSim/quadprog/QuadProg++.hh
 
  // Matrix in quadratic form of objective
  QuadProgPP::Matrix<scalar> G;
  // Equality constraints
  QuadProgPP::Matrix<scalar> CE;
  // Inequality constraints
  QuadProgPP::Matrix<scalar> CI;
  
  QuadProgPP::Vector<scalar> g0;
  QuadProgPP::Vector<scalar> ce0;
  QuadProgPP::Vector<scalar> ci0;
  QuadProgPP::Vector<scalar> x;

  // M = 16  0 0
  //      0 16 0
  //      0  0 289.75
  G.resize(3,3);
  for(int i=0; i< 3; i++) for(int j=0; j< 3; j++) G[i][j]=0;
  G[0][0] = 16; 
  G[1][1] = 16;
  G[2][2] = 289.75;

  // -M \dot q = -0
  //             139.52
  //            -0
  g0.resize(3);
  g0[0] = 0;
  g0[1] = 139.52;
  g0[2] = 0;

  // No equality constraints, currently
  CE.resize(3,0);
  ce0.resize(0);

  // Compute the number of inequality constraints
  CI.resize(3,24);

  MatrixXs tempN(3,24);
  tempN << 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           1,  1,  1,  1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          -5, -4, -3, -2, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  for( int i = 0; i < 3; ++i ) for( int j = 0; j < 24; ++j ) CI[i][j] = tempN(i,j);

  // Constant term added to inequality constraints
  ci0.resize(24);
  for(int i=0; i< 24; i++) ci0[i] = 0;

  // Solution
  x.resize(3);

  solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
  */

