#include "SimpleCollisionHandler.h"
#include <iostream>
#include <set>

// BEGIN STUDENT CODE //


// Detects whether two particles are overlapping (including the radii of each)
// and approaching.
// If the two particles overlap and are approaching, returns true and sets 
// the vector n to be the vector between the first and second particle.
// Inputs:
//   scene: The scene data structure. The positions and radii of the particles
//          can be obtained from here.
//   idx1:  The index of the first particle. (Ie, the degrees of freedom
//          corresponding to this particle are entries 2*idx1 and 2*idx1+1 in
//          scene.getX().
//   idx2:  The index of the second particle.
// Outputs:
//   n: The vector between the two particles.
//   Returns true if the two particles overlap and are approaching.
bool SimpleCollisionHandler::detectParticleParticle(TwoDScene &scene, int idx1, int idx2, Vector2s &n)
{
    VectorXs x1 = scene.getX().segment<2>(2*idx1);
    VectorXs x2 = scene.getX().segment<2>(2*idx2);
    
    // your implementation here
    // find the closest vector
    VectorXs v1 = scene.getV().segment<2>(2*idx1);
    VectorXs v2 = scene.getV().segment<2>(2*idx2);

    const scalar& r1 = scene.getRadius(idx1);
    const scalar& r2 = scene.getRadius(idx2);
    scalar distance2 = (x1 - x2).squaredNorm();

    if (distance2 < (r1 + r2)*(r1 + r2))
    {
        if((x1 - x2).dot(v1 - v2) < 0 ) 
        {
            n = (x2 - x1);
            return true;
        }
        else
        {            
            return false;
        }   
    }
    else
    {
        return false;
    }
}

// Detects whether a particle and an edge are overlapping (including the radii 
// of both) and are approaching.
// If the two objects overlap and are approaching, returns true and sets the 
// vector n to be the shortest vector between the particle and the edge.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   eidx:  The index of the edge. (Ie, the indices of particle with index e are
//          scene.getEdges()[e].first and scene.getEdges()[e].second.)
// Outputs:
//   n: The shortest vector between the particle and the edge.
//   Returns true if the two objects overlap and are approaching.
bool SimpleCollisionHandler::detectParticleEdge(TwoDScene &scene, int vidx, int eidx, Vector2s &n)
{
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    VectorXs x2 = scene.getX().segment<2>(2*scene.getEdges()[eidx].first);
    VectorXs x3 = scene.getX().segment<2>(2*scene.getEdges()[eidx].second);
    
    // your implementation here
    // compute the closest vector between the two bodies
    scalar alpha = (x1 - x2).dot(x3 - x2) / (x3 - x2).squaredNorm();
    alpha = alpha > 1.0 ? 1.0 : (alpha < 0.0 ? 0.0 : alpha); // clamp to [0,1]
        
    VectorXs xalpha = x2 + alpha * (x3 - x2);
    n = xalpha - x1;

    scalar re = scene.getEdgeRadii()[eidx];
    scalar r1 = scene.getRadius(vidx); 

    if(n.squaredNorm() < (r1 + re)*(r1 + re)) // overlap
    {
        // compute valpha
        VectorXs v1 = scene.getV().segment<2>(2*vidx);
        VectorXs v2 = scene.getV().segment<2>(2*scene.getEdges()[eidx].first);
        VectorXs v3 = scene.getV().segment<2>(2*scene.getEdges()[eidx].second);
        VectorXs valpha = v2 + alpha * (v3 - v2);

        if ( (valpha - v1).dot(n) < 0 ) // approaching
        {
            return true;
        }
        else
        {
            return false;
        }
             
    }
    else
    {
        return false;
    }
    
}

// Detects whether a particle and a half-plane are overlapping (including the 
// radius of the particle) and are approaching.
// If the two objects overlap and are approaching, returns true and sets the 
// vector n to be the shortest vector between the particle and the half-plane.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   pidx:  The index of the halfplane. The vectors (px, py) and (nx, ny) can
//          be retrieved by calling scene.getHalfplane(pidx).
// Outputs:
//   n: The shortest vector between the particle and the half-plane.
//   Returns true if the two objects overlap and are approaching.
bool SimpleCollisionHandler::detectParticleHalfplane(TwoDScene &scene, int vidx, int pidx, Vector2s &n)
{
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    VectorXs px = scene.getHalfplane(pidx).first;
    VectorXs pn = scene.getHalfplane(pidx).second;
    
    // your implementation here
    scalar r = scene.getRadius(vidx);
    VectorXs v1 = scene.getV().segment<2>(2 * vidx); // velocity particle
    pn.normalize();
    n = (px - x1).dot(pn) * pn;

    if ( n.squaredNorm() <  r * r )// overlap
    {
        if( v1.dot(pn) < 0 ) // approaching
        {
            return true;
        }
        else
        {
            return false;
        }
        
    }
    else
    {
        return false;
    }    
}


// Responds to a collision detected between two particles by applying an impulse
// to the velocities of each one.
// You can get the COR of the simulation by calling getCOR().
// Inputs:
//   scene: The scene data structure.
//   idx1:  The index of the first particle.
//   idx2:  The index of the second particle.
//   n:     The vector between the first and second particle.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleParticle(TwoDScene &scene, int idx1, int idx2, const Vector2s &n)
{
    const VectorXs &M = scene.getM();
    VectorXs &v = scene.getV();
    
    // your implementation here
    scalar distance2 = n.squaredNorm();
    
    scalar m1 = M[2*idx1];
    scalar m2 = M[2*idx2];
    
    VectorXs v1 = v.segment<2>(2*idx1);
    VectorXs v2 = v.segment<2>(2*idx2);

    scalar cor_coeff = 0.5 * (1.0 + this->getCOR());
    VectorXs impulse_numerator =  cor_coeff * (2.0 * (v2 - v1).dot(n)) * n / (distance2);
    
    v.segment<2>(2*idx1) +=  (!scene.isFixed(idx1)) * impulse_numerator / (1.0 +   (!scene.isFixed(idx2)) * m1/m2);
    v.segment<2>(2*idx2) += -(!scene.isFixed(idx2)) * impulse_numerator / (1.0 +   (!scene.isFixed(idx1)) * m2/m1);

}

// Responds to a collision detected between a particle and an edge by applying
// an impulse to the velocities of each one.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   eidx:  The index of the edge.
//   n:     The shortest vector between the particle and the edge.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleEdge(TwoDScene &scene, int vidx, int eidx, const Vector2s &n)
{
    const VectorXs &M = scene.getM();
    
    int eidx1 = scene.getEdges()[eidx].first;
    int eidx2 = scene.getEdges()[eidx].second;
    
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    VectorXs x2 = scene.getX().segment<2>(2*eidx1);
    VectorXs x3 = scene.getX().segment<2>(2*eidx2);
    
    VectorXs v1 = scene.getV().segment<2>(2*vidx);
    VectorXs v2 = scene.getV().segment<2>(2*eidx1);
    VectorXs v3 = scene.getV().segment<2>(2*eidx2);
    
    // your implementation here
    //compute alpha
    scalar alpha = (x1 - x2).dot(x3 - x2) / (x3 - x2).squaredNorm();
    alpha = alpha > 1 ? 1 : (alpha < 0 ? 0 : alpha); // clamp to [0,1]

    VectorXs xalpha = x2 + alpha * (x3 - x2);
    // compute valpha
    VectorXs valpha = v2 + alpha * (v3 - v2);
    
    scalar cor_coeff = 0.5 * (1.0 + this->getCOR());
    VectorXs impulse_numerator = cor_coeff * 2.0 * (valpha - v1).dot(n) * n / (n.squaredNorm());
    scalar m1 = M[2*vidx];
    scalar me1 = M[2*eidx1];
    scalar me2 = M[2*eidx2];

    scene.getV().segment<2>(2*vidx)  +=   (!scene.isFixed(vidx))  * impulse_numerator / (1 + ((1 - alpha) * (1- alpha) * (!scene.isFixed(eidx1)) *  m1 / me1) + (alpha * alpha * (!scene.isFixed(eidx2)) * m1 / me2));
    scene.getV().segment<2>(2*eidx1) += - (!scene.isFixed(eidx1)) * (1 - alpha)*impulse_numerator / ( (!scene.isFixed(vidx)) * me1/m1 + (1-alpha)*(1-alpha) + alpha*alpha*(!scene.isFixed(eidx2))*me1/me2 );
    scene.getV().segment<2>(2*eidx2) += - (!scene.isFixed(eidx2)) * alpha*impulse_numerator/( (!scene.isFixed(vidx)) * me2/m1 + (1-alpha)*(1-alpha)* (!scene.isFixed(eidx1)) * me2 / me1 + alpha*alpha );

}


// Responds to a collision detected between a particle and a half-plane by 
// applying an impulse to the velocity of the particle.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   pidx:  The index of the half-plane.
//   n:     The shortest vector between the particle and the half-plane.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleHalfplane(TwoDScene &scene, int vidx, int pidx, const Vector2s &n)
{
    VectorXs nhat = n;
    
    // your implementation here
    nhat.normalize();
    VectorXs v = scene.getV().segment<2>(2*vidx);
    
    scalar cor_coeff = 0.5 * (1.0 + this->getCOR());
    scene.getV().segment<2>(2*vidx) += - cor_coeff * 2 * (v.dot(nhat)) * nhat / scene.getM()[2*vidx];
}
