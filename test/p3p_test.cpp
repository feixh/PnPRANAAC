//
// Created by feixh on 8/29/15.
//

#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>

#include <pnpsolvers/P3p.h>

#include <Eigen/Core>
#include <Eigen/StdVector>

using namespace std;
using namespace Eigen;

int main()
{
    ifstream ifs("../test/data.txt", ifstream::in );
    assert( ifs.is_open() );
    vector< Vector3d > pc;
    vector< Vector3d > pts;
    int n, n_inliers;
    ifs >> n >> n_inliers;
    for ( int i = 0; i < n; ++i )
    {
        float x, y;
        ifs >> x >> y;
        pc.emplace_back( x, y, 1.0 );
        pc.back().normalize();
    }
    for ( int i = 0; i < n; ++i )
    {
        float x, y, z;
        ifs >> x >> y >> z;
        pts.emplace_back(x, y, z);
        pts.back().normalize();
    }
    float gwc[3][4];
    for ( int i = 0; i < 3; ++i ) {
        for (int j = 0; j < 4; ++j) {
            ifs >> gwc[i][j];
        }
    }
    ifs.close();
    P3p p3p;
    Matrix3d featureVectors;
    Matrix3d worldPoints;
    for ( int i = 0; i < 3; ++i )
    {
        featureVectors.col(i) = pc[i];
        worldPoints.col(i) = pts[i];
    }
    vector< Matrix<double,3,4> > solutions;
    p3p.computePoses( featureVectors, worldPoints, solutions );
    for ( int i = 0; i < solutions.size(); ++i )
    {
        cout << solutions[i] << endl;
    }
}