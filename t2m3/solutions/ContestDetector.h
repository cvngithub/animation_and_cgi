#ifndef CONTEST_DETECTOR_H
#define CONTEST_DETECTOR_H

#define aDEBUG

#ifdef DEBUG
#define LOG(x) std::cout << x 
#else
#define LOG(x)
#endif

#include "CollisionDetector.h"
#include <vector>
#include <set>

typedef std::set<std::pair<int, int> > PPList;
typedef std::set<std::pair<int, int> > PEList;
typedef std::set<std::pair<int, int> > PHList;


class ContestDetector : public CollisionDetector
{
 public:
  ContestDetector();

  virtual void performCollisionDetection(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, DetectionCallback &dc);

 private:
  void findCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs);

  void naiveFindCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs);



public:
  // Fixed grid with spatial hash (See: Teschner et al. 2003)
  struct SHCell
  {
    std::vector<int> objects;
    std::vector<int> edges;  
    std::vector<int> halfplanes;  
  };
  void spatialHashingFindCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs);
  inline int hash(int i, int j, int p1, int p2, int hashsize)
  {
    return ((i * p1) ^ (j * p2)) % hashsize;   
  }


//  struct BoundingBox
//  {
//    scalar xmin;
//    scalar xmax;
//    scalar ymin;
//    scalar ymax;
//    int    id;
//  };
//  void lagrangianFindCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs);


};

//static std::vector<ContestDetector::BoundingBox> boundingBoxesXcoords;
//static std::vector<ContestDetector::BoundingBox> boundingBoxesYcoords;



static std::vector<ContestDetector::SHCell> SHCells;
static bool start = true;
static const int ndim = 2;

static scalar cellGrid=0.075;
static int p1 = 73856093;
static int p2 = 19349663;
static int p3 = 83492791;
static int hashsize = 25000;

//static scalar cellGrid=0.5;
//static int p1 = 73856093;
//static int p2 = 19349663;
//static int p3 = 83492791;
//static int hashsize = 10;
#endif
