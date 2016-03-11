//
// Created by feixh on 8/29/15.
//

// STL
#include <iostream>
#include <fstream>
#include <vector>

// theia
#include <theia/solvers/ransac.h>
#include <theia/solvers/estimator.h>
#include <theia/solvers/random_sampler.h>
#include <theia/util/timer.h>


// eigen
#include <Eigen/Core>
#include <Eigen/StdVector>

// P3P
#include "ransac_estimators.h"

using namespace std;

int main()
{
    // load data
    ifstream ifs("../test/data.txt", ifstream::in );
    assert( ifs.is_open() );
    vector< Eigen::Vector3d > pc;
    vector< Eigen::Vector3d > pts;
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
    vector< ransac_estimators::Match2D3D > data( n );
    for ( size_t i = 0; i < n; ++i )
    {
        data[i].featureVector = pc[i];
        data[i].worldPoint    = pts[i];
    }
    cout << "got " << n << " matches" << endl;

    // setup ransac parameters
    ransac_estimators::RansacParameters ransac_params;
    ransac_params.error_thresh = 1e-2;
    ransac_params.failure_probability = 0.05;
    ransac_params.max_iterations = 3000;
    ransac_params.min_inlier_ratio = 0.1;
    ransac_params.use_mle = true;


    ransac_estimators::P3PEstimator estimator;
     ransac_estimators::Ransac< ransac_estimators::P3PEstimator > ransac(ransac_params, estimator);
    ransac.Initialize();
    ransac_estimators::RansacSummary summary;
    Eigen::Matrix< double, 3, 4 > best_model;
    ransac_estimators::Timer tt;
    ransac.Estimate( data, &best_model, &summary );
    double duration( tt.ElapsedTimeInSeconds() );
    cout << duration << " s" << endl;

    // checkout results
    cout << best_model << endl;
    cout << "iterations:" << summary.num_iterations << endl;
    for ( size_t i = 0; i < summary.inliers.size(); ++i )
    {
        cout << summary.inliers[i] << " ";
    }
    cout << endl;
}
