#include "ExplicitEuler.h"

bool ExplicitEuler::stepScene( TwoDScene& scene, scalar dt )
{
    // Your code goes here!
    
    // Some tips on getting data from TwoDScene:
    // A vector containing all of the system's position DoFs. x0, y0, x1, y1, ...
    //VectorXs& x = scene.getX();
    // A vector containing all of the system's velocity DoFs. v0, v0, v1, v1, ...
    //VectorXs& v = scene.getV();
    // A vector containing the masses associated to each DoF. m0, m0, m1, m1, ...
    //const VectorXs& m = scene.getM();
    // Determine if the ith particle is fixed
    // if( scene.isFixed(i) )
    VectorXs& x = scene.getX();
    VectorXs& v = scene.getV();
    const VectorXs& m = scene.getM();
	VectorXs updatedX = x + (dt*v);
    
    VectorXs F(x.size());
    F.setZero(x.size());
    scene.accumulateGradU(F);
    for (int i=0; i<scene.getNumParticles(); i++) {
        if (scene.isFixed(i)) {
            F(2*i) = 0;
            F(2*i + 1) = 0;
        }
    }    
    VectorXs updatedV = v + (dt*m.cwiseInverse().cwiseProduct(-F));
    x = updatedX;
    v = updatedV;    
    return true;
}
