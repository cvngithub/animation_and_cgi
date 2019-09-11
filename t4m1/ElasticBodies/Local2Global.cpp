#include "Local2Global.h"

void Local2Global(const std::vector<int> &i_points, const MatrixXs &hessE_local, MatrixXs &hessE_global) {

    int ilocal1 = 0;
    for (int iglob1 = 0; iglob1 < (int) i_points.size(); ++iglob1) {
        int ilocal2 = 0;
        for (int iglob2 = 0; iglob2 < (int) i_points.size(); ++iglob2) {

            hessE_global.block(i_points[iglob1] * 2, i_points[iglob2] * 2, 2, 2) +=\
            hessE_local.block(ilocal1 * 2, ilocal2 * 2, 2, 2);
            ilocal2++;
        }
        ilocal1++;
    }
}

/*
void Local2Global( const std::pair<int,int> p_endpoints, const MatrixXs& hessE_local, MatrixXs& hessE_global)
{
    hessE_global.block(p_endpoints.first * 2, p_endpoints.first * 2, 2, 2) += hessE_local.block(0, 0, 2, 2);
    hessE_global.block(p_endpoints.first * 2, p_endpoints.second * 2, 2, 2) += hessE_local.block(0, 2, 2, 2);
    hessE_global.block(p_endpoints.second * 2, p_endpoints.first * 2, 2, 2) += hessE_local.block(2, 0, 2, 2);
    hessE_global.block(p_endpoints.second * 2, p_endpoints.second * 2, 2, 2) += hessE_local.block(2, 2, 2, 2);
}*/
