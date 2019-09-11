#include "ElasticBodyCSTForce.h"
#include "../MathUtilities.h"
#include "Local2Global.h"
#include <assert.h>

/*
Python script to produce energy, its gradient and hessian

from sympy import simplify, symbols, numbered_symbols, sqrt, cse, ccode, simplify, atan2


def main():
    # Produce the code for the energy

    x = symbols('x(i:k)((0:2))', real=True)
    print(x), print()

    B, e, f, g, h = symbols('B, e, f, g, h', real=True)
    a, b, c, d = symbols('a, b, c, d', real=True)
    a = (x[2] - x[0]) * e + (x[4] - x[0]) * f
    b = (x[3] - x[1]) * e + (x[5] - x[1]) * f
    c = (x[2] - x[0]) * g + (x[4] - x[0]) * h
    d = (x[3] - x[1]) * g + (x[5] - x[1]) * h

    m_E, m_nu = symbols('m_E, m_nu', real=True)

    exx, eyy, exy = symbols('exx, eyy, exy', real=True)
    exx = a * a + b * b - 1.0
    eyy = c * c + d * d - 1.0
    exy = a * c + b * d

    E = B * (exx * exx + 2.0 * m_nu * exx * eyy + eyy * eyy + 2.0 * (1.0 - m_nu) * exy * exy)

    variables_sigma = numbered_symbols('sigma_')
    replacements, reduced = cse(E, symbols=variables_sigma)

    for key, val in replacements:
        print('scalar', key, '=', ccode(val) + ';')

    for i, r in enumerate(reduced):
        print('E_update({})'.format(i), '=', ccode(r) + ';')

    print()

    # Produce the gradient

    partials = [simplify(E.diff(var)) for var in x]

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
        print('hessE_update({row},{column})'.format(row=int(i / 6), column=i % 6), '=', ccode(r) + ';')


main() *
*/

void ElasticBodyCSTForce::addEnergyToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, scalar &E) {
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

    scalar A, B, e, f, g, h;
    ABefgh(A, B, e, f, g, h);

    const scalar m_E = m_youngs_modulus;
    const scalar m_nu = m_poisson_ratio;

    scalar sigma_0 = -xi(0);
    scalar sigma_1 = sigma_0 + xj(0);
    scalar sigma_2 = sigma_0 + xk(0);
    scalar sigma_3 = e * sigma_1 + f * sigma_2;
    scalar sigma_4 = -xi(1);
    scalar sigma_5 = sigma_4 + xj(1);
    scalar sigma_6 = sigma_4 + xk(1);
    scalar sigma_7 = e * sigma_5 + f * sigma_6;
    scalar sigma_8 = pow(sigma_3, 2) + pow(sigma_7, 2) - 1.0;
    scalar sigma_9 = g * sigma_1 + h * sigma_2;
    scalar sigma_10 = g * sigma_5 + h * sigma_6;
    scalar sigma_11 = pow(sigma_10, 2) + pow(sigma_9, 2) - 1.0;
    scalar sigma_12 = 2.0 * m_nu;

    E += B * (pow(sigma_11, 2) + sigma_11 * sigma_12 * sigma_8 + pow(sigma_8, 2) +
              (-sigma_12 + 2.0) * pow(sigma_10 * sigma_7 + sigma_3 * sigma_9, 2));

}

void ElasticBodyCSTForce::addGradEToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, VectorXs &gradE) {
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

    VectorXs gradE_update(6);

    scalar A, B, e, f, g, h;
    ABefgh(A, B, e, f, g, h);

    const scalar m_E = m_youngs_modulus;
    const scalar m_nu = m_poisson_ratio;

    scalar sigma_0 = xi(0) - xj(0);
    scalar sigma_1 = xi(0) - xk(0);
    scalar sigma_2 = e * sigma_0 + f * sigma_1;
    scalar sigma_3 = xi(1) - xj(1);
    scalar sigma_4 = xi(1) - xk(1);
    scalar sigma_5 = e * sigma_3 + f * sigma_4;
    scalar sigma_6 = pow(sigma_2, 2) + pow(sigma_5, 2) - 1.0;
    scalar sigma_7 = 4 * sigma_6;
    scalar sigma_8 = e + f;
    scalar sigma_9 = sigma_2 * sigma_8;
    scalar sigma_10 = g * sigma_0 + h * sigma_1;
    scalar sigma_11 = g * sigma_3 + h * sigma_4;
    scalar sigma_12 = pow(sigma_10, 2) + pow(sigma_11, 2) - 1.0;
    scalar sigma_13 = 4 * sigma_12;
    scalar sigma_14 = g + h;
    scalar sigma_15 = sigma_10 * sigma_14;
    scalar sigma_16 = 4.0 * m_nu;
    scalar sigma_17 = sigma_12 * sigma_16;
    scalar sigma_18 = sigma_16 * sigma_6;
    scalar sigma_19 = -xi(0);
    scalar sigma_20 = sigma_19 + xj(0);
    scalar sigma_21 = sigma_19 + xk(0);
    scalar sigma_22 = 4.0 * (m_nu - 1) * (sigma_10 * sigma_2 + sigma_11 * sigma_5);
    scalar sigma_23 = sigma_5 * sigma_7;
    scalar sigma_24 = sigma_11 * sigma_13;
    scalar sigma_25 = sigma_17 * sigma_5;
    scalar sigma_26 = sigma_11 * sigma_18;
    scalar sigma_27 = -xi(1);
    scalar sigma_28 = sigma_27 + xj(1);
    scalar sigma_29 = sigma_27 + xk(1);
    scalar sigma_30 = sigma_2 * sigma_7;
    scalar sigma_31 = sigma_10 * sigma_13;
    scalar sigma_32 = sigma_17 * sigma_2;
    scalar sigma_33 = sigma_10 * sigma_18;

    gradE_update(0) = B * (sigma_13 * sigma_15 + sigma_15 * sigma_18 + sigma_17 * sigma_9 + \
    sigma_22 * (sigma_14 * (e * sigma_20 + f * sigma_21) + sigma_8 * (g * sigma_20 + h * sigma_21)) + \
    sigma_7 * sigma_9);
    gradE_update(1) = B * (sigma_14 * sigma_24 + sigma_14 * sigma_26 + sigma_22 *
                                                                       (sigma_14 * (e * sigma_28 + f * sigma_29) +
                                                                        sigma_8 * (g * sigma_28 + h * sigma_29)) +
                           sigma_23 * sigma_8 + sigma_25 * sigma_8);
    gradE_update(2) =
            -B * (e * sigma_30 + e * sigma_32 + g * sigma_31 + g * sigma_33 - sigma_22 * (e * sigma_10 + g * sigma_2));
    gradE_update(3) =
            -B * (e * sigma_23 + e * sigma_25 + g * sigma_24 + g * sigma_26 - sigma_22 * (e * sigma_11 + g * sigma_5));
    gradE_update(4) =
            -B * (f * sigma_30 + f * sigma_32 + h * sigma_31 + h * sigma_33 - sigma_22 * (f * sigma_10 + h * sigma_2));
    gradE_update(5) =
            -B * (f * sigma_23 + f * sigma_25 + h * sigma_24 + h * sigma_26 - sigma_22 * (f * sigma_11 + h * sigma_5));

    gradE.segment<2>(2 * m_idx1) += gradE_update.segment<2>(0);
    gradE.segment<2>(2 * m_idx2) += gradE_update.segment<2>(2);
    gradE.segment<2>(2 * m_idx3) += gradE_update.segment<2>(4);

}

void ElasticBodyCSTForce::addHessXToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, MatrixXs &hessE) {
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

    MatrixXs hessE_update(6, 6);

    scalar A, B, e, f, g, h;
    ABefgh(A, B, e, f, g, h);

    const scalar m_E = m_youngs_modulus;
    const scalar m_nu = m_poisson_ratio;

    scalar sigma_0 = -xi(0);
    scalar sigma_1 = sigma_0 + xj(0);
    scalar sigma_2 = sigma_0 + xk(0);
    scalar sigma_3 = e * sigma_1 + f * sigma_2;
    scalar sigma_4 = pow(sigma_3, 2);
    scalar sigma_5 = 2 * e;
    scalar sigma_6 = 2 * f;
    scalar sigma_7 = -sigma_5 - sigma_6;
    scalar sigma_8 = 4 * e;
    scalar sigma_9 = 4 * f;
    scalar sigma_10 = -sigma_8 - sigma_9;
    scalar sigma_11 = sigma_10 * sigma_7;
    scalar sigma_12 = g * sigma_1 + h * sigma_2;
    scalar sigma_13 = pow(sigma_12, 2);
    scalar sigma_14 = 2 * g;
    scalar sigma_15 = 2 * h;
    scalar sigma_16 = -sigma_14 - sigma_15;
    scalar sigma_17 = 4 * g;
    scalar sigma_18 = 4 * h;
    scalar sigma_19 = -sigma_17 - sigma_18;
    scalar sigma_20 = sigma_16 * sigma_19;
    scalar sigma_21 = -e - f;
    scalar sigma_22 = sigma_12 * sigma_21;
    scalar sigma_23 = -g - h;
    scalar sigma_24 = sigma_23 * sigma_3;
    scalar sigma_25 = sigma_22 + sigma_24;
    scalar sigma_26 = 2.0 * m_nu;
    scalar sigma_27 = -sigma_26 + 2.0;
    scalar sigma_28 = sigma_27 * (2 * sigma_22 + 2 * sigma_24);
    scalar sigma_29 = sigma_12 * sigma_3;
    scalar sigma_30 = 4.0 * m_nu;
    scalar sigma_31 = sigma_29 * sigma_30;
    scalar sigma_32 = sigma_16 * sigma_7;
    scalar sigma_33 = -xi(1);
    scalar sigma_34 = sigma_33 + xj(1);
    scalar sigma_35 = sigma_33 + xk(1);
    scalar sigma_36 = e * sigma_34 + f * sigma_35;
    scalar sigma_37 = g * sigma_34 + h * sigma_35;
    scalar sigma_38 = sigma_36 * sigma_37;
    scalar sigma_39 = sigma_27 * (sigma_29 + sigma_38);
    scalar sigma_40 = pow(sigma_36, 2);
    scalar sigma_41 = sigma_4 + sigma_40 - 1.0;
    scalar sigma_42 = sigma_10 * sigma_41;
    scalar sigma_43 = pow(sigma_37, 2);
    scalar sigma_44 = sigma_13 + sigma_43 - 1.0;
    scalar sigma_45 = sigma_19 * sigma_44;
    scalar sigma_46 = sigma_26 * sigma_7;
    scalar sigma_47 = sigma_44 * sigma_46;
    scalar sigma_48 = sigma_16 * sigma_26 * sigma_41;
    scalar sigma_49 = sigma_21 * sigma_42 + sigma_21 * sigma_47 + sigma_23 * sigma_45 + sigma_23 * sigma_48 +
                      sigma_39 * (sigma_16 * sigma_21 + sigma_23 * sigma_7);
    scalar sigma_50 = sigma_21 * sigma_37;
    scalar sigma_51 = sigma_23 * sigma_36;
    scalar sigma_52 = sigma_50 + sigma_51;
    scalar sigma_53 = sigma_3 * sigma_36;
    scalar sigma_54 = sigma_12 * sigma_37;
    scalar sigma_55 = sigma_16 * sigma_46;
    scalar sigma_56 = sigma_3 * sigma_37;
    scalar sigma_57 = sigma_12 * sigma_36;
    scalar sigma_58 = 2 * pow(sigma_16, 2) * sigma_54 + 2 * sigma_53 * pow(sigma_7, 2) + sigma_55 * sigma_56 +
                      sigma_55 * sigma_57;
    scalar sigma_59 = e * sigma_12;
    scalar sigma_60 = g * sigma_3;
    scalar sigma_61 = sigma_59 + sigma_60;
    scalar sigma_62 = sigma_10 * sigma_4;
    scalar sigma_63 = sigma_13 * sigma_19;
    scalar sigma_64 = e * sigma_16;
    scalar sigma_65 = g * sigma_7;
    scalar sigma_66 = sigma_31 * sigma_64 + sigma_31 * sigma_65;
    scalar sigma_67 = e * sigma_42 + e * sigma_47 + g * sigma_45 + g * sigma_48 + sigma_39 * (sigma_64 + sigma_65);
    scalar sigma_68 = e * sigma_37;
    scalar sigma_69 = g * sigma_36;
    scalar sigma_70 = sigma_68 + sigma_69;
    scalar sigma_71 = sigma_30 * sigma_57;
    scalar sigma_72 = sigma_30 * sigma_56;
    scalar sigma_73 = sigma_7 * sigma_8;
    scalar sigma_74 = sigma_16 * sigma_17;
    scalar sigma_75 = sigma_53 * sigma_73 + sigma_54 * sigma_74;
    scalar sigma_76 = sigma_64 * sigma_71 + sigma_65 * sigma_72 + sigma_75;
    scalar sigma_77 = f * sigma_12;
    scalar sigma_78 = h * sigma_3;
    scalar sigma_79 = sigma_77 + sigma_78;
    scalar sigma_80 = f * sigma_16;
    scalar sigma_81 = h * sigma_7;
    scalar sigma_82 = sigma_31 * sigma_80 + sigma_31 * sigma_81;
    scalar sigma_83 = f * sigma_42 + f * sigma_47 + h * sigma_45 + h * sigma_48 + sigma_39 * (sigma_80 + sigma_81);
    scalar sigma_84 = f * sigma_37;
    scalar sigma_85 = h * sigma_36;
    scalar sigma_86 = sigma_84 + sigma_85;
    scalar sigma_87 = sigma_7 * sigma_9;
    scalar sigma_88 = sigma_16 * sigma_18;
    scalar sigma_89 = sigma_53 * sigma_87 + sigma_54 * sigma_88;
    scalar sigma_90 = sigma_71 * sigma_80 + sigma_72 * sigma_81 + sigma_89;
    scalar sigma_91 = sigma_27 * (2 * sigma_50 + 2 * sigma_51);
    scalar sigma_92 = sigma_30 * sigma_38;
    scalar sigma_93 = sigma_64 * sigma_72 + sigma_65 * sigma_71 + sigma_75;
    scalar sigma_94 = sigma_10 * sigma_40;
    scalar sigma_95 = sigma_19 * sigma_43;
    scalar sigma_96 = sigma_64 * sigma_92 + sigma_65 * sigma_92;
    scalar sigma_97 = sigma_71 * sigma_81 + sigma_72 * sigma_80 + sigma_89;
    scalar sigma_98 = sigma_80 * sigma_92 + sigma_81 * sigma_92;
    scalar sigma_99 = sigma_27 * (sigma_12 * sigma_5 + sigma_14 * sigma_3);
    scalar sigma_100 = sigma_41 * sigma_8;
    scalar sigma_101 = sigma_17 * sigma_44;
    scalar sigma_102 = sigma_30 * sigma_44;
    scalar sigma_103 = sigma_102 * sigma_21;
    scalar sigma_104 = sigma_30 * sigma_41;
    scalar sigma_105 = sigma_104 * sigma_23;
    scalar sigma_106 = e * sigma_103 + g * sigma_105 + sigma_100 * sigma_21 + sigma_101 * sigma_23 +
                       sigma_39 * (sigma_14 * sigma_21 + sigma_23 * sigma_5);
    scalar sigma_107 = pow(e, 2);
    scalar sigma_108 = 8 * sigma_107;
    scalar sigma_109 = pow(g, 2);
    scalar sigma_110 = 8 * sigma_109;
    scalar sigma_111 = e * m_nu;
    scalar sigma_112 = sigma_111 * sigma_29;
    scalar sigma_113 = 16.0 * g;
    scalar sigma_114 = 4 * sigma_41;
    scalar sigma_115 = 4 * sigma_44;
    scalar sigma_116 = g * sigma_39 * sigma_8 + sigma_102 * sigma_107 + sigma_104 * sigma_109 + sigma_107 * sigma_114 +
                       sigma_109 * sigma_115;
    scalar sigma_117 = 8.0 * m_nu;
    scalar sigma_118 = sigma_117 * sigma_68;
    scalar sigma_119 = sigma_117 * sigma_69;
    scalar sigma_120 = sigma_108 * sigma_53 + sigma_110 * sigma_54 + sigma_118 * sigma_60 + sigma_119 * sigma_59;
    scalar sigma_121 = e * f;
    scalar sigma_122 = 8 * sigma_121;
    scalar sigma_123 = g * h;
    scalar sigma_124 = 8 * sigma_123;
    scalar sigma_125 = 8.0 * h;
    scalar sigma_126 = f * sigma_29;
    scalar sigma_127 = g * sigma_117;
    scalar sigma_128 = f * sigma_100 + h * sigma_101 + sigma_102 * sigma_121 + sigma_104 * sigma_123 +
                       sigma_39 * (g * sigma_6 + h * sigma_5);
    scalar sigma_129 =
            sigma_112 * sigma_125 + sigma_122 * sigma_4 + sigma_124 * sigma_13 + sigma_126 * sigma_127 + sigma_128;
    scalar sigma_130 = sigma_122 * sigma_53 + sigma_124 * sigma_54;
    scalar sigma_131 = sigma_118 * sigma_78 + sigma_119 * sigma_77 + sigma_130;
    scalar sigma_132 = sigma_27 * (sigma_14 * sigma_36 + sigma_37 * sigma_5);
    scalar sigma_133 = sigma_111 * sigma_38;
    scalar sigma_134 = sigma_117 * sigma_85;
    scalar sigma_135 = sigma_117 * sigma_84;
    scalar sigma_136 = sigma_130 + sigma_134 * sigma_59 + sigma_135 * sigma_60;
    scalar sigma_137 = f * sigma_38;
    scalar sigma_138 =
            sigma_122 * sigma_40 + sigma_124 * sigma_43 + sigma_125 * sigma_133 + sigma_127 * sigma_137 + sigma_128;
    scalar sigma_139 = sigma_27 * (sigma_12 * sigma_6 + sigma_15 * sigma_3);
    scalar sigma_140 = f * sigma_103 + h * sigma_105 + sigma_18 * sigma_23 * sigma_44 + sigma_21 * sigma_41 * sigma_9 +
                       sigma_39 * (sigma_15 * sigma_21 + sigma_23 * sigma_6);
    scalar sigma_141 = pow(f, 2);
    scalar sigma_142 = 8 * sigma_141;
    scalar sigma_143 = pow(h, 2);
    scalar sigma_144 = 8 * sigma_143;
    scalar sigma_145 = 16.0 * h * m_nu;
    scalar sigma_146 = h * sigma_39 * sigma_9 + sigma_102 * sigma_141 + sigma_104 * sigma_143 + sigma_114 * sigma_141 +
                       sigma_115 * sigma_143;
    scalar sigma_147 = sigma_134 * sigma_77 + sigma_135 * sigma_78 + sigma_142 * sigma_53 + sigma_144 * sigma_54;
    scalar sigma_148 = sigma_27 * (sigma_15 * sigma_36 + sigma_37 * sigma_6);

    hessE_update(0, 0) =
            B * (sigma_11 * sigma_4 + sigma_13 * sigma_20 + sigma_25 * sigma_28 + sigma_31 * sigma_32 + sigma_49);
    hessE_update(0, 1) = B * (sigma_28 * sigma_52 + sigma_58);
    hessE_update(0, 2) = B * (sigma_14 * sigma_63 + sigma_28 * sigma_61 + sigma_5 * sigma_62 + sigma_66 + sigma_67);
    hessE_update(0, 3) = B * (sigma_28 * sigma_70 + sigma_76);
    hessE_update(0, 4) = B * (sigma_15 * sigma_63 + sigma_28 * sigma_79 + sigma_6 * sigma_62 + sigma_82 + sigma_83);
    hessE_update(0, 5) = B * (sigma_28 * sigma_86 + sigma_90);
    hessE_update(1, 0) = B * (sigma_25 * sigma_91 + sigma_58);
    hessE_update(1, 1) =
            B * (sigma_11 * sigma_40 + sigma_20 * sigma_43 + sigma_32 * sigma_92 + sigma_49 + sigma_52 * sigma_91);
    hessE_update(1, 2) = B * (sigma_61 * sigma_91 + sigma_93);
    hessE_update(1, 3) = B * (sigma_14 * sigma_95 + sigma_5 * sigma_94 + sigma_67 + sigma_70 * sigma_91 + sigma_96);
    hessE_update(1, 4) = B * (sigma_79 * sigma_91 + sigma_97);
    hessE_update(1, 5) = B * (sigma_15 * sigma_95 + sigma_6 * sigma_94 + sigma_83 + sigma_86 * sigma_91 + sigma_98);
    hessE_update(2, 0) = B * (sigma_106 + sigma_13 * sigma_74 + sigma_25 * sigma_99 + sigma_4 * sigma_73 + sigma_66);
    hessE_update(2, 1) = B * (sigma_52 * sigma_99 + sigma_93);
    hessE_update(2, 2) =
            B * (sigma_108 * sigma_4 + sigma_110 * sigma_13 + sigma_112 * sigma_113 + sigma_116 + sigma_61 * sigma_99);
    hessE_update(2, 3) = B * (sigma_120 + sigma_70 * sigma_99);
    hessE_update(2, 4) = B * (sigma_129 + sigma_79 * sigma_99);
    hessE_update(2, 5) = B * (sigma_131 + sigma_86 * sigma_99);
    hessE_update(3, 0) = B * (sigma_132 * sigma_25 + sigma_76);
    hessE_update(3, 1) = B * (sigma_106 + sigma_132 * sigma_52 + sigma_40 * sigma_73 + sigma_43 * sigma_74 + sigma_96);
    hessE_update(3, 2) = B * (sigma_120 + sigma_132 * sigma_61);
    hessE_update(3, 3) = B * (sigma_108 * sigma_40 + sigma_110 * sigma_43 + sigma_113 * sigma_133 + sigma_116 +
                              sigma_132 * sigma_70);
    hessE_update(3, 4) = B * (sigma_132 * sigma_79 + sigma_136);
    hessE_update(3, 5) = B * (sigma_132 * sigma_86 + sigma_138);
    hessE_update(4, 0) = B * (sigma_13 * sigma_88 + sigma_139 * sigma_25 + sigma_140 + sigma_4 * sigma_87 + sigma_82);
    hessE_update(4, 1) = B * (sigma_139 * sigma_52 + sigma_97);
    hessE_update(4, 2) = B * (sigma_129 + sigma_139 * sigma_61);
    hessE_update(4, 3) = B * (sigma_136 + sigma_139 * sigma_70);
    hessE_update(4, 4) =
            B * (sigma_126 * sigma_145 + sigma_13 * sigma_144 + sigma_139 * sigma_79 + sigma_142 * sigma_4 + sigma_146);
    hessE_update(4, 5) = B * (sigma_139 * sigma_86 + sigma_147);
    hessE_update(5, 0) = B * (sigma_148 * sigma_25 + sigma_90);
    hessE_update(5, 1) = B * (sigma_140 + sigma_148 * sigma_52 + sigma_40 * sigma_87 + sigma_43 * sigma_88 + sigma_98);
    hessE_update(5, 2) = B * (sigma_131 + sigma_148 * sigma_61);
    hessE_update(5, 3) = B * (sigma_138 + sigma_148 * sigma_70);
    hessE_update(5, 4) = B * (sigma_147 + sigma_148 * sigma_79);
    hessE_update(5, 5) = B * (sigma_137 * sigma_145 + sigma_142 * sigma_40 + sigma_144 * sigma_43 + sigma_146 +
                              sigma_148 * sigma_86);

    Local2Global(i_points, hessE_update, hessE);


}

void ElasticBodyCSTForce::addHessVToTotal(const VectorXs &x, const VectorXs &v, const VectorXs &m, MatrixXs &hessE) {
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

void ElasticBodyCSTForce::ABefgh(scalar &A, scalar &B, scalar &e, scalar &f, scalar &g, scalar &h) {

    const scalar m_E = m_youngs_modulus;
    const scalar m_nu = m_poisson_ratio;

    Matrix2s mxb, mxb_inverse;
    mxb(0, 0) = m_xb2(0) - m_xb1(0);
    mxb(0, 1) = m_xb3(0) - m_xb1(0);
    mxb(1, 0) = m_xb2(1) - m_xb1(1);
    mxb(1, 1) = m_xb3(1) - m_xb1(1);

    mxb_inverse = mxb.inverse();
    e = mxb_inverse(0, 0);
    f = mxb_inverse(1, 0);
    g = mxb_inverse(0, 1);
    h = mxb_inverse(1, 1);

    A = 0.5 * fabs(mathutils::crossTwoD(m_xb3 - m_xb1, m_xb2 - m_xb1));
    B = (0.5 * m_E * A) / (1.0 - m_nu * m_nu);

}
