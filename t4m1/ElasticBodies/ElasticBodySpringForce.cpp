#include "ElasticBodySpringForce.h"
#include "Local2Global.h"
#include <assert.h>

/*
Python script to produce energy, its gradient and hessian

# from sympy.printing.cxxcode import cxxcode
from sympy import collect, expand, symbols, numbered_symbols, sqrt, cse, ccode
import numpy as np


def main():
    # Produce the code for the energy

    x = symbols('x(i:j)((0:2))', real=True)

    print(x)

    print()

    E, m_alpha, m_l0, norm, norm2 = symbols('w, m_alpha, m_l0, norm, norm2', real=True)

    norm2 = (x[0] - x[2]) ** 2 + (x[1] - x[3]) ** 2

    norm = sqrt(norm2)

    E = m_alpha / (2.0 * m_l0) * (norm - m_l0) ** 2

    print('scalar E =', ccode(E) + ';')

    print()

    # Produce the gradient

    partials = [E.diff(var) for var in x]

    variables_sigma = numbered_symbols('sigma_')
    replacements, reduced = cse(partials, symbols=variables_sigma)

    for key, val in replacements:
        print('scalar', key, '=', ccode(val) + ';')

    print()

    for i, r in enumerate(reduced):
        print('gradE_update({})'.format(i), '=', ccode(r) + ';')

    # Produce the hessian

    partials_second = [E.diff(ix).diff(iy) for ix in x for iy in x]

    print()

    variables_sigma = numbered_symbols('sigma_')
    replacements, reduced = cse(partials_second, symbols=variables_sigma)

    for key, val in replacements:
        print('scalar', key, '=', ccode(val) + ';')

    print()

    for i, r in enumerate(reduced):
        print('hessE_update({row},{column})'.format(row=int(i/4), column=i%4), '=', ccode(r) + ';')


main()
 */

void ElasticBodySpringForce::addEnergyToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, scalar &E) {
    assert(x.size() == v.size());
    assert(x.size() % 2 == 0);
    assert(m_idx1 >= 0);
    assert(m_idx1 < x.size() / 2);
    assert(m_idx2 >= 0);
    assert(m_idx2 < x.size() / 2);

    const Vector2s xi = x.segment<2>(2 * m_idx1);
    const Vector2s xj = x.segment<2>(2 * m_idx2);

    E += 0.5 * m_alpha * pow(-m_l0 + sqrt(pow(xi(0) - xj(0), 2) + pow(xi(1) - xj(1), 2)), 2) / m_l0;

}

void ElasticBodySpringForce::addGradEToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, VectorXs &gradE) {
    assert(x.size() == v.size());
    assert(x.size() == gradE.size());
    assert(x.size() % 2 == 0);
    assert(m_idx1 >= 0);
    assert(m_idx1 < x.size() / 2);
    assert(m_idx2 >= 0);
    assert(m_idx2 < x.size() / 2);

    const Vector2s xi = x.segment<2>(2 * m_idx1);
    const Vector2s xj = x.segment<2>(2 * m_idx2);

    VectorXs gradE_update(4);

    scalar sigma_0 = xi(0) - xj(0);
    scalar sigma_1 = xi(1) - xj(1);
    scalar sigma_2 = sqrt(sigma_0 * sigma_0 + sigma_1 * sigma_1);
    scalar sigma_3 = m_alpha * (-m_l0 + sigma_2) / (m_l0 * sigma_2);

    gradE_update(0) = sigma_0 * sigma_3;
    gradE_update(1) = sigma_1 * sigma_3;
    gradE_update(2) = -gradE_update(0);
    gradE_update(3) = -gradE_update(1);

    gradE.segment<2>(2 * m_idx1) += gradE_update.segment<2>(0);
    gradE.segment<2>(2 * m_idx2) += gradE_update.segment<2>(2);

}

void ElasticBodySpringForce::addHessXToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, MatrixXs &hessE) {
    assert(x.size() == v.size());
    assert(x.size() == m.size());
    assert(x.size() == hessE.rows());
    assert(x.size() == hessE.cols());
    assert(x.size() % 2 == 0);
    assert(m_idx1 >= 0);
    assert(m_idx1 < x.size() / 2);
    assert(m_idx2 >= 0);
    assert(m_idx2 < x.size() / 2);

    const Vector2s xi = x.segment<2>(2 * m_idx1);
    const Vector2s xj = x.segment<2>(2 * m_idx2);

    std::vector<int> i_points;
    i_points.push_back(m_idx1);
    i_points.push_back(m_idx2);

    MatrixXs hessE_update(4, 4);

    scalar sigma_0 = xi(0) - xj(0);
    scalar sigma_1 = pow(sigma_0, 2);
    scalar sigma_2 = xi(1) - xj(1);
    scalar sigma_3 = pow(sigma_2, 2);
    scalar sigma_4 = sigma_1 + sigma_3;
    scalar sigma_5 = 1.0*m_alpha/m_l0;
    scalar sigma_6 = sigma_5/sigma_4;
    scalar sigma_7 = sqrt(sigma_4);
    scalar sigma_8 = sigma_5*(-m_l0 + sigma_7);
    scalar sigma_9 = sigma_8/sigma_7;
    scalar sigma_10 = -xi(0) + xj(0);
    scalar sigma_11 = sigma_8/pow(sigma_4, 3.0/2.0);
    scalar sigma_12 = sigma_0*sigma_11;
    scalar sigma_13 = sigma_10*sigma_12 + sigma_9;
    scalar sigma_14 = sigma_0*sigma_6;
    scalar sigma_15 = sigma_14*sigma_2;
    scalar sigma_16 = -xi(1) + xj(1);
    scalar sigma_17 = sigma_12*sigma_16;
    scalar sigma_18 = -sigma_9;
    scalar sigma_19 = sigma_10*sigma_14 + sigma_18;
    scalar sigma_20 = sigma_14*sigma_16;
    scalar sigma_21 = sigma_12*sigma_2;
    scalar sigma_22 = sigma_11*sigma_2;
    scalar sigma_23 = sigma_10*sigma_22;
    scalar sigma_24 = sigma_16*sigma_22 + sigma_9;
    scalar sigma_25 = sigma_2*sigma_6;
    scalar sigma_26 = sigma_10*sigma_25;
    scalar sigma_27 = sigma_16*sigma_25 + sigma_18;
    scalar sigma_28 = pow(sigma_10, 2);
    scalar sigma_29 = sigma_10*sigma_16;
    scalar sigma_30 = sigma_11*sigma_29;
    scalar sigma_31 = sigma_29*sigma_6;
    scalar sigma_32 = pow(sigma_16, 2);

    hessE_update(0,0) = sigma_1*sigma_6 + sigma_13;
    hessE_update(0,1) = sigma_15 + sigma_17;
    hessE_update(0,2) = sigma_1*sigma_11 + sigma_19;
    hessE_update(0,3) = sigma_20 + sigma_21;
    hessE_update(1,0) = sigma_15 + sigma_23;
    hessE_update(1,1) = sigma_24 + sigma_3*sigma_6;
    hessE_update(1,2) = sigma_21 + sigma_26;
    hessE_update(1,3) = sigma_11*sigma_3 + sigma_27;
    hessE_update(2,0) = sigma_11*sigma_28 + sigma_19;
    hessE_update(2,1) = sigma_26 + sigma_30;
    hessE_update(2,2) = sigma_13 + sigma_28*sigma_6;
    hessE_update(2,3) = sigma_23 + sigma_31;
    hessE_update(3,0) = sigma_20 + sigma_30;
    hessE_update(3,1) = sigma_11*sigma_32 + sigma_27;
    hessE_update(3,2) = sigma_17 + sigma_31;
    hessE_update(3,3) = sigma_24 + sigma_32*sigma_6;

    Local2Global(i_points, hessE_update, hessE);

}

void ElasticBodySpringForce::addHessVToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, MatrixXs &hessE) {
    assert(x.size() == v.size());
    assert(x.size() == m.size());
    assert(x.size() == hessE.rows());
    assert(x.size() == hessE.cols());
    assert(x.size() % 2 == 0);
    assert(m_idx1 >= 0);
    assert(m_idx1 < x.size() / 2);
    assert(m_idx2 >= 0);
    assert(m_idx2 < x.size() / 2);

}
