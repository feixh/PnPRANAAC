//
// Created by feixh on 8/29/15.
//

#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>

#include "pnpsolvers/P3P_Kneip.h"
#include "pnpsolvers/P3p.h"

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <TooN/TooN.h>
#include <pnpsolvers/P3p.h>

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
    }
    float gwc[3][4];
    for ( int i = 0; i < 3; ++i ) {
        for (int j = 0; j < 4; ++j) {
            ifs >> gwc[i][j];
        }
    }
    ifs.close();
    P3p p3p;
    TooN::Matrix<3,3> f;
    TooN::Matrix<3,3> w;
    TooN::Matrix<3,16> sol;
    for ( int i = 0; i < 3; ++i )
    {
        for ( int j = 0; j < 3; ++j ){
            f(j,i) = pc[i](j);
            w(j,i) = pts[i](j);
        }
    }
    p3p.computePoses( f, w, sol );
    cout << sol << endl;



    // eigen version of Kneip's P3P algorithm
    cout << "===" << endl;
    P3P_Kneip kp3p;
    Matrix3d featureVectors;
    Matrix3d worldPoints;
    for ( int i = 0; i < 3; ++i )
    {
        featureVectors.col(i) = pc[i];
        worldPoints.col(i) = pts[i];
    }
    vector< Matrix<double,3,4> > solutions;
    kp3p.computePoses( featureVectors, worldPoints, solutions );
    for ( int i = 0; i < solutions.size(); ++i )
    {
        cout << solutions[i] << endl;
    }
}