#include "PenaltyForce.h"
#include "TwoDScene.h"
#define LOG(x) std::cout << x << std::endl

void PenaltyForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
    // Feel free to implement if you feel like doing so.
}

// BEGIN STUDENT CODE


// Adds the gradient of the penalty potential (-1 * force) for a pair of 
// particles to the total.
// Read the positions of the particles from the input variable x. Radii can
// be obtained from the member variable m_scene, the penalty force stiffness 
// from member variable m_k, and penalty force thickness from member variable
// m_thickness.
// Inputs:
//   x:    The positions of the particles in the scene. 
//   idx1: The index of the first particle, i.e. the position of this particle
//         is ( x[2*idx1], x[2*idx1+1] ).
//   idx2: The index of the second particle.
// Outputs:
//   gradE: The total gradient of penalty force. *ADD* the particle-particle
//          gradient to this total gradient.
void PenaltyForce::addParticleParticleGradEToTotal(const VectorXs &x, int idx1, int idx2, VectorXs &gradE)
{
    VectorXs x1 = x.segment<2>(2*idx1);
    VectorXs x2 = x.segment<2>(2*idx2);
    
    double r1 = m_scene.getRadius(idx1);
    double r2 = m_scene.getRadius(idx2);
    
    // your implementation here

    scalar distance = (x2 - x1).norm();
    VectorXs nhat = (x2 - x1) / distance;

//    Easy way, no matrices
//    VectorXs toAddGradE = (  distance <= r1 + r2 + m_thickness ? 1 : 0) * m_k * (distance -r1 -r2 -m_thickness) * nhat;
//    gradE.segment<2>(2*idx1) += -toAddGradE;
//    gradE.segment<2>(2*idx2) += toAddGradE;

    // dimensionality
    int dim = x1.size();
    // compute matrix: (grad n)T. It is a 2d x d matrix (here d = 2 is the spatial dim)
    MatrixXs grad_n(dim, 2*dim);
    grad_n.block<2,2>(0,0) = -Matrix2s::Identity();
    grad_n.block<2,2>(0,2) = Matrix2s::Identity();

    VectorXs localGradE(4);
    localGradE.setZero();
    if(distance <= r1 + r2 + m_thickness)
    {
        localGradE = (m_k * (distance -(r1 + r2 + m_thickness)) * grad_n.transpose() * nhat);
    }    

    gradE.segment<2>(2*idx1) += (!m_scene.isFixed(idx1)) * localGradE.segment<2>(0);
    gradE.segment<2>(2*idx2) += (!m_scene.isFixed(idx2)) * localGradE.segment<2>(dim);

}

// Adds the gradient of the penalty potential (-1 * force) for a particle-edge
// pair to the total.
// Read the positions of the particle and edge endpoints from the input
// variable x.
// Inputs:
//   x:    The positions of the particles in the scene.
//   vidx: The index of the particle.
//   eidx: The index of the edge, i.e. the indices of the particle making up the
//         endpoints of the edge are given by m_scene.getEdge(eidx).first and 
//         m_scene.getEdges(eidx).second.
// Outputs:
//   gradE: The total gradient of penalty force. *ADD* the particle-edge
//          gradient to this total gradient.
void PenaltyForce::addParticleEdgeGradEToTotal(const VectorXs &x, int vidx, int eidx, VectorXs &gradE)
{
    VectorXs x1 = x.segment<2>(2*vidx);
    VectorXs x2 = x.segment<2>(2*m_scene.getEdge(eidx).first);
    VectorXs x3 = x.segment<2>(2*m_scene.getEdge(eidx).second);
    
    if(m_scene.getEdge(eidx).first != vidx && m_scene.getEdge(eidx).second != vidx)
    {
        double r1 = m_scene.getRadius(vidx);
        double r2 = m_scene.getEdgeRadii()[eidx];

        // your implementation here
        // your implementation here
        // compute alpha


        scalar alpha = (x1 - x2).dot(x3 - x2) / (x3 - x2).squaredNorm();
        alpha = alpha > 1.0 ? 1.0 : (alpha < 0.0 ? 0.0 : alpha); // clamp to [0,1]

        VectorXs xalpha = x2 + alpha * (x3 - x2);

        scalar distance = (xalpha - x1).norm();
        VectorXs nhat = (xalpha - x1) / distance;

        // dimensionality
        int dim = x1.size();
        // compute matrix: (grad n)T. It is a 2d x d matrix (here d = 2 is the spatial dim)
        MatrixXs grad_n(dim, 3*dim);
        grad_n.block<2,2>(0,0)     = -Matrix2s::Identity();
        grad_n.block<2,2>(0,dim)   = (1 - alpha) * Matrix2s::Identity();
        grad_n.block<2,2>(0,2*dim) = alpha * Matrix2s::Identity();


        VectorXs localGradE(3*dim);
        localGradE.setZero();
        if(distance <= r1 + r2 + m_thickness)
        {
            localGradE = (m_k * (distance -(r1 + r2 + m_thickness)) * grad_n.transpose() * nhat);
        }    

        gradE.segment<2>(2*vidx)                            += (!m_scene.isFixed(vidx)) * localGradE.segment<2>(0);
        gradE.segment<2>(2*m_scene.getEdges()[eidx].first)  += (!m_scene.isFixed(m_scene.getEdges()[eidx].first)) * localGradE.segment<2>(dim);
        gradE.segment<2>(2*m_scene.getEdges()[eidx].second) += (!m_scene.isFixed(m_scene.getEdges()[eidx].second)) * localGradE.segment<2>(2*dim);
    }

}

// Adds the gradient of the penalty potential (-1 * force) for a particle-
// half-plane pair to the total.
// Read the positions of the particle from the input variable x.
// Inputs:
//   x:    The positions of the particles in the scene.
//   vidx: The index of the particle.
//   pidx: The index of the half-plane, i.e. the position and normal vectors
//         for the half-plane can be retrieved by calling
//         m_scene.getHalfplane(pidx).
// Outputs:
//   gradE: The total gradient of the penalty force. *ADD* the particle-
//          half-plane gradient to this total gradient.
void PenaltyForce::addParticleHalfplaneGradEToTotal(const VectorXs &x, int vidx, int pidx, VectorXs &gradE)
{
    VectorXs x1 = x.segment<2>(2*vidx);
    VectorXs nh = m_scene.getHalfplane(pidx).second;
    
    // your implementation here

    VectorXs toAddGradE(2);
    toAddGradE.setZero();
    VectorXs px = m_scene.getHalfplane(pidx).first;
    scalar r1 = m_scene.getRadius(vidx);
    
    nh.normalize();
    VectorXs n = (px - x1).dot(nh) * nh;
    scalar dist = n.norm();
    if ( dist <=  (r1 + m_thickness)  )// overlap
    {
        toAddGradE = m_k * (dist - r1 - m_thickness) * nh;
    }
    gradE.segment<2>(2*vidx) += (!m_scene.isFixed(vidx)) * toAddGradE;

}
