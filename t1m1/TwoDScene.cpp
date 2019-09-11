#include "TwoDScene.h"

scalar TwoDScene::computeKineticEnergy() const
{
    return 0.5 * m_v.transpose().cwiseProduct(m_m.transpose()) * m_v;
}
