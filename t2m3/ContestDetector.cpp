#include "ContestDetector.h"
#include <iostream>
#include "TwoDScene.h"
#include <set>

ContestDetector::ContestDetector()
{    
}


// Given particle positions, computes lists of *potentially* overlapping object
// pairs. How exactly to do this is up to you.
// Inputs: 
//   scene:  The scene object. Get edge information, radii, etc. from here. If 
//           for some reason you'd also like to use particle velocities in your
//           algorithm, you can get them from here too.
//   x:      The positions of the particle.
// Outputs:
//   pppairs: A list of (particle index, particle index) pairs of potentially
//            overlapping particles. IMPORTANT: Each pair should only appear
//            in the list at most once. (1, 2) and (2, 1) count as the same 
//            pair.
//   pepairs: A list of (particle index, edge index) pairs of potential
//            particle-edge overlaps.
//   phpairs: A list of (particle index, halfplane index) pairs of potential
//            particle-halfplane overlaps.
void ContestDetector::findCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs)
{
  //naiveFindCollidingPairs(scene,x, pppairs, pepairs, phpairs);
  spatialHashingFindCollidingPairs(scene, x, pppairs, pepairs, phpairs);
}



void ContestDetector::spatialHashingFindCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs)
{
    // reset colliding pairs
    pppairs.clear();
    pepairs.clear();
    phpairs.clear();
    
    // retrieve all the edges
    std::vector<std::pair<int,int> > edges = scene.getEdges();
    std::vector<scalar> edgesRads = scene.getEdgeRadii();

    // retrieve all the half planes
    std::vector<std::pair<VectorXs, VectorXs> > halfPlanes = scene.getHalfplanes();
    

    if(start)
    {
        // initialize SpatialHashing algorithm's cells
        SHCells.resize(hashsize);
        start = false;
    }

    for(auto it = SHCells.begin(); it != SHCells.end(); ++it)
    {
        (*it).objects.resize(0);
        (*it).edges.resize(0);
        (*it).halfplanes.resize(0);

    }

    const std::vector<double>& r = scene.getRadii();

    #ifdef DEBUG
    LOG("Radiuses: ");
    for (size_t i = 0; i < r.size(); i++)
    {
        LOG( r[i] << " " ); 
    }
    LOG(std::endl);
    #endif

    //for(int h = 0; h < halfPlanes.size(); h++)
    //{
    //    
    //}

    for(int e = 0; e < edges.size(); e++)
    {
        // object's minimum dimension
        int x1min = (int) floor(   (x[ndim * edges.at(e).first + 0]  - edgesRads.at(e))/cellGrid   );
        int y1min = (int) floor(   (x[ndim * edges.at(e).first + 1]  - edgesRads.at(e))/cellGrid   );
        int x2min = (int) floor(   (x[ndim * edges.at(e).second + 0] - edgesRads.at(e))/cellGrid   );
        int y2min = (int) floor(   (x[ndim * edges.at(e).second + 1] - edgesRads.at(e))/cellGrid   );
        int x1max = (int) floor(   (x[ndim * edges.at(e).first + 0]  + edgesRads.at(e))/cellGrid   );
        int y1max = (int) floor(   (x[ndim * edges.at(e).first + 1]  + edgesRads.at(e))/cellGrid   );
        int x2max = (int) floor(   (x[ndim * edges.at(e).second + 0] + edgesRads.at(e))/cellGrid   );
        int y2max = (int) floor(   (x[ndim * edges.at(e).second + 1] + edgesRads.at(e))/cellGrid   );

        int x_idx_min = std::min(x1min, x2min);
        int y_idx_min = std::min(y1min, y2min);
        int x_idx_max = std::max(x1max, x2max);
        int y_idx_max = std::max(y1max, y2max);


        LOG( "esges min max idx:" << abs(x_idx_max - x_idx_min) << " and " << abs(y_idx_max - y_idx_min) << std::endl ); 


        for(int i = 0; i <= abs(x_idx_max - x_idx_min); i++)
        {
            for(int j = 0; j <= abs(y_idx_max - y_idx_min); j++)
            {
                int cellIdx = abs(hash(x_idx_min+i, y_idx_min+j, p1, p2, hashsize));
                LOG( "i=" << i << ",j=" << j << " | cell idx = " << cellIdx );
                LOG( " size there: " << SHCells.at( cellIdx ).edges.size() << std::endl );

                if(SHCells.at(cellIdx).edges.size())
                {
                    if(SHCells.at(cellIdx).edges.back()!=e)
                    {
                        SHCells.at( cellIdx ).edges.push_back(e);
                        LOG( "->inserted edge #" << e << std::endl );
                    }
                }else{
                        SHCells.at( cellIdx ).edges.push_back(e);
                        LOG( "->inserted edge #" << e << std::endl );
                }


            }
        }


    }

    for(int o = 0; o < x.size(); o+=ndim)
    {
        LOG( "***********object #" << o/ndim << ", out of " << x.size()/ndim << "objects"<< std::endl );
                
        int x_idx_min = (int) floor((x[o+0] - r[o/ndim]) / cellGrid);
        int y_idx_min = (int) floor((x[o+1] - r[o/ndim]) / cellGrid);
        int x_idx_max = (int) floor((x[o+0] + r[o/ndim]) / cellGrid);
        int y_idx_max = (int) floor((x[o+1] + r[o/ndim]) / cellGrid); 

        LOG( "ixmin = " << x_idx_min << "," << "ximax = " <<x_idx_max << "," << "iymin = " <<y_idx_min << "," << "yimax = " <<y_idx_max << std::endl );
        LOG( "    abs(ixmax-ixmin) = " <<abs(x_idx_max - x_idx_min) << std::endl);

        for(int i = 0; i <= abs(x_idx_max - x_idx_min); i++)
        {
            for(int j = 0; j <= abs(y_idx_max - y_idx_min); j++)
            {
                int cellIdx = abs(hash(x_idx_min+i, y_idx_min+j, p1, p2, hashsize));
                LOG( "i=" << i << ",j=" << j << " | cell idx = " << cellIdx );
                LOG( " size there: " << SHCells.at( cellIdx ).objects.size() << std::endl );

                if(SHCells.at(cellIdx).objects.size())
                {
                    if(SHCells.at(cellIdx).objects.back()!=o/2)
                    {
                        SHCells.at( cellIdx ).objects.push_back(o/2);
                        LOG( "->inserted object #" << o/2 << std::endl );
                    }
                }else{
                        SHCells.at( cellIdx ).objects.push_back(o/2);
                        LOG( "->inserted object #" << o/2 << std::endl );
                }


            }
        }
    }



    for (int c = 0; c < SHCells.size(); c++)
    {

        std::vector<int>& objects   = SHCells.at(c).objects;
        std::vector<int>& cellEdges = SHCells.at(c).edges;


        LOG( "**************cell " << c << ", objects.size() = " << objects.size() << std::endl );
        LOG( "objects = " );

#ifdef DEBUG
        for(auto it = objects.begin(); it != objects.end(); ++it)
        {
            LOG( (*it)<< " ");
        }
        LOG(std::endl);
#endif
        LOG( "edges = " );
#ifdef DEBUG
        for(auto it = cellEdges.begin(); it != cellEdges.end(); ++it)
        {
            LOG( (*it)<< " ");
        }
        LOG(std::endl);
#endif
        for(int i = 0; i < objects.size(); i++)
        {
            for(int j = i + 1; j < objects.size(); j++)
            {
                LOG( objects.at(i) << "," << objects.at(j) << "? :"); 
                std::pair<int, int> lol(objects.at(i), objects.at(j));
                if(lol.first != lol.second)
                {
                    pppairs.insert(lol);
                    LOG(" inserted PP pair." << std::endl);
                }
            }  
        }

        for(int i = 0; i < objects.size(); i++)
        {
            for(int j = 0; j < cellEdges.size(); j++)
            {
                LOG( objects.at(i) << "," << cellEdges.at(j) << "? :"); 
                std::pair<int, int> lol(objects.at(i), cellEdges.at(j));
                pepairs.insert(lol);
                LOG(" inserted PE pair." << std::endl);
            }  
        }


        for(int i = 0; i < objects.size(); i++)
        {   
            for (int j = 0; j < scene.getNumHalfplanes(); j++)
            {
                std::pair<int, int> lol(objects.at(i), j);
                phpairs.insert(lol);
            }
        }

        LOG(std::endl);
    }

#ifdef DEBUG
    for (auto it = pppairs.begin(); it != pppairs.end(); ++it )
    {
        LOG( "(" << (*it).first << "," << (*it).second << ")-");
    }
    LOG(std::endl);
#endif

}



//void ContestDetector::lagrangianFindCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs)
//{
//    // find the BB's min and max for each axis
//    // sort them using insertion sort or bubble sort
//    //boundingBoxesXcoords;
//
//    if(start)
//    {
//        // prepare the first
//        boundingBoxesXcoords.resize(x.size()/ndim);
//        boundingBoxesYcoords.resize(x.size()/ndim);
//    }
//
//    for(int o = 0; o < x.size(); o+=ndim)
//    {
//        LOG( "***********object #" << o/ndim << ", out of " << x.size()/ndim << "objects"<< std::endl );
//                
//        ((x[o+0] - r[o/ndim]) / cellGrid);
//        ((x[o+1] - r[o/ndim]) / cellGrid);
//        ((x[o+0] + r[o/ndim]) / cellGrid);
//        ((x[o+1] + r[o/ndim]) / cellGrid); 
//    }
//    
//
//}
//



void ContestDetector::naiveFindCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs)
{
    int ndim = 2;
    // O(n^2) algorithm
    for (int n = 0; n < x.size()/ndim; n++)
    {
        for (int m = n+1; m < x.size()/ndim; m++)
        {
            // pair n,m
            std::pair<int, int> pair(n,m);
            pppairs.insert(pair);
        }
    }
}


