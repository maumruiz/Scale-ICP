#ifndef _ICP_H_
#define _ICP_H_

#include "Mesh.h"
#include <vector>
#include <iostream>
using namespace std;

#include "Eigen/Dense"
using namespace Eigen;

class ICP
{
public:
    vector<pair<Vector4f, float **>> transformations;
    ICP(){};
    
    void align(Mesh *mesh1, Mesh *mesh2, int nMaxIters, bool oneToOne, float minDisplacement, bool prescale);
    void transform(Mesh*);

private:
    float calculateMaxEucDist();
};

#endif