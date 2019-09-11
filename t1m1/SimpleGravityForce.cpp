#include "SimpleGravityForce.h"

void SimpleGravityForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size()%2 == 0 );
    
    // Your code goes here!
    
    VectorXs m_gravity_dash = m_gravity.replicate(x.size()/2, 1);
    E +=  -m_gravity_dash.transpose().cwiseProduct(m.transpose()) * x; 
}

void SimpleGravityForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size() == gradE.size() );
    assert( x.size()%2 == 0 );
    
    // Your code goes here!
    
    VectorXs m_gravity_dash = m_gravity.replicate(x.size()/2, 1);
    gradE += -m_gravity_dash.transpose().cwiseProduct(m.transpose());
}
