#include "ElasticBodyBendingForce.h"
#include "Local2Global.h"
#include <assert.h>

void ElasticBodyBendingForce::addEnergyToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, scalar &E) {
    assert(x.size() == v.size());
    assert(x.size() % 2 == 0);
    assert(m_idx1 >= 0);
    assert(m_idx1 < x.size() / 2);
    assert(m_idx2 >= 0);
    assert(m_idx2 < x.size() / 2);
    assert(m_idx3 >= 0);
    assert(m_idx3 < x.size() / 2);

    const Vector2s xi = x.segment<2>(2 * m_idx1);
    const Vector2s xj = x.segment<2>(2 * m_idx2);
    const Vector2s xk = x.segment<2>(2 * m_idx3);

    const scalar B = m_alpha / (2.0 * (m_eb1n + m_eb2n));

    scalar sigma_0 = -xi(0) + xj(0);
    scalar sigma_1 = -xj(1) + xk(1);
    scalar sigma_2 = -xi(1) + xj(1);
    scalar sigma_3 = -xj(0) + xk(0);

    E += B * pow(-m_theta0 + atan2(sigma_0 * sigma_1 - sigma_2 * sigma_3, sigma_0 * sigma_3 + sigma_1 * sigma_2), 2);

}

void
ElasticBodyBendingForce::addGradEToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, VectorXs &gradE) {
    assert(x.size() == v.size());
    assert(x.size() == gradE.size());
    assert(x.size() % 2 == 0);
    assert(m_idx1 >= 0);
    assert(m_idx1 < x.size() / 2);
    assert(m_idx2 >= 0);
    assert(m_idx2 < x.size() / 2);
    assert(m_idx3 >= 0);
    assert(m_idx3 < x.size() / 2);

    const Vector2s xi = x.segment<2>(2 * m_idx1);
    const Vector2s xj = x.segment<2>(2 * m_idx2);
    const Vector2s xk = x.segment<2>(2 * m_idx3);

    const scalar B = m_alpha / (2.0 * (m_eb1n + m_eb2n));

    VectorXs gradE_cross(6);
    gradE_cross(0) = xj(1) - xk(1);
    gradE_cross(1) = -xj(0) + xk(0);
    gradE_cross(2) = -xi(1) + xk(1);
    gradE_cross(3) = xi(0) - xk(0);
    gradE_cross(4) = xi(1) - xj(1);
    gradE_cross(5) = -xi(0) + xj(0);

    VectorXs gradE_dot(6);
    gradE_dot(0) = xj(0) - xk(0);
    gradE_dot(1) = xj(1) - xk(1);
    gradE_dot(2) = xi(0) - 2 * xj(0) + xk(0);
    gradE_dot(3) = xi(1) - 2 * xj(1) + xk(1);
    gradE_dot(4) = -xi(0) + xj(0);
    gradE_dot(5) = -xi(1) + xj(1);

    VectorXs gradE_update(6);

    scalar sigma_0 = mathutils::crossTwoD(xj - xi, xk - xj);
    scalar sigma_1 = (xj - xi).dot(xk - xj);
    scalar sigma_2 = 2 / (pow(sigma_0, 2) + pow(sigma_1, 2));
    scalar sigma_3 = B * (-m_theta0 + atan2(sigma_0, sigma_1)) * sigma_2;

    gradE_update = (sigma_3 * sigma_1) * gradE_cross - (sigma_3 * sigma_0) * gradE_dot;

    gradE.segment<2>(2 * m_idx1) += gradE_update.segment<2>(0);
    gradE.segment<2>(2 * m_idx2) += gradE_update.segment<2>(2);
    gradE.segment<2>(2 * m_idx3) += gradE_update.segment<2>(4);

}

void
ElasticBodyBendingForce::addHessXToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, MatrixXs &hessE) {
    assert(x.size() == v.size());
    assert(x.size() == m.size());
    assert(x.size() == hessE.rows());
    assert(x.size() == hessE.cols());
    assert(x.size() % 2 == 0);
    assert(m_idx1 >= 0);
    assert(m_idx1 < x.size() / 2);
    assert(m_idx2 >= 0);
    assert(m_idx2 < x.size() / 2);
    assert(m_idx3 >= 0);
    assert(m_idx3 < x.size() / 2);

    const Vector2s xi = x.segment<2>(2 * m_idx1);
    const Vector2s xj = x.segment<2>(2 * m_idx2);
    const Vector2s xk = x.segment<2>(2 * m_idx3);

    std::vector<int> i_points;
    i_points.push_back(m_idx1);
    i_points.push_back(m_idx2);
    i_points.push_back(m_idx3);

    const scalar B = m_alpha / (2.0 * (m_eb1n + m_eb2n));

    // Gradient of (xj - xi ) x (xk - xj)
    VectorXs gradE_cross(6);
    gradE_cross(0) = xj(1) - xk(1);
    gradE_cross(1) = -xj(0) + xk(0);
    gradE_cross(2) = -xi(1) + xk(1);
    gradE_cross(3) = xi(0) - xk(0);
    gradE_cross(4) = xi(1) - xj(1);
    gradE_cross(5) = -xi(0) + xj(0);

    // Gradient of (xj - xi ) . (xk - xj)
    VectorXs gradE_dot(6);
    gradE_dot(0) = xj(0) - xk(0);
    gradE_dot(1) = xj(1) - xk(1);
    gradE_dot(2) = xi(0) - 2 * xj(0) + xk(0);
    gradE_dot(3) = xi(1) - 2 * xj(1) + xk(1);
    gradE_dot(4) = -xi(0) + xj(0);
    gradE_dot(5) = -xi(1) + xj(1);

    // Hessian of (xj - xi ) x (xk - xj)
    MatrixXs hessE_cross(6, 6);
    hessE_cross
            << 0, 0, 0, 1, 0, -1, 0, 0, -1, 0, 1, 0, 0, -1, 0, 0, 0, 1, 1, 0, 0, 0, -1, 0, 0, 1, 0, -1, 0, 0, -1, 0, 1, 0, 0, 0;

    // Hessian of (xj - xi ) . (xk - xj)
    MatrixXs hessE_dot(6, 6);
    hessE_dot
            << 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, 0, -1, 1, 0, -2, 0, 1, 0, 0, 1, 0, -2, 0, 1, -1, 0, 1, 0, 0, 0, 0, -1, 0, 1, 0, 0;

    MatrixXs hessE_update(6, 6);

    scalar sigma_0 = mathutils::crossTwoD(xj - xi, xk - xj);
    scalar sigma_1 = (xj - xi).dot(xk - xj);
    scalar sigma_2 = 1.0 / (pow(sigma_0, 2) + pow(sigma_1, 2));
    scalar sigma_3 = 2 * sigma_2;
    scalar sigma_t = m_theta0 - atan2(sigma_0, sigma_1);

    VectorXs vector_4 = sigma_0 * gradE_dot;
    VectorXs vector_5 = sigma_1 * gradE_cross;
    VectorXs vector_6 = vector_4 - vector_5;

    MatrixXs matrix_1 = sigma_0 * hessE_dot + sigma_1 * hessE_cross;

    for (int i = 0; i < 6; ++i) {
        scalar sigma_4 = vector_6(i);
        for (int j = 0; j < 6; ++j) {
            if (j <= i) {
                scalar sigma_6 = sigma_3 * (gradE_cross(j) * sigma_0 + gradE_dot(j) * sigma_1);
                hessE_update(i, j) = B * sigma_3 *
                                     (sigma_2 * sigma_4 *
                                      (-gradE_cross(j) * sigma_1 + gradE_dot(j) * sigma_0) -
                                      sigma_t * (gradE_cross(i) * gradE_dot(j) - gradE_cross(j) * gradE_dot(i) -
                                                 matrix_1(i, j) + sigma_6 * sigma_4));
                if (j < i) hessE_update(j, i) = hessE_update(i, j);
            }
        }
    };

    Local2Global(i_points, hessE_update, hessE);
}

void
ElasticBodyBendingForce::addHessVToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, MatrixXs &hessE) {
    assert(x.size() == v.size());
    assert(x.size() == m.size());
    assert(x.size() == hessE.rows());
    assert(x.size() == hessE.cols());
    assert(x.size() % 2 == 0);
    assert(m_idx1 >= 0);
    assert(m_idx1 < x.size() / 2);
    assert(m_idx2 >= 0);
    assert(m_idx2 < x.size() / 2);
    assert(m_idx3 >= 0);
    assert(m_idx3 < x.size() / 2);


}
