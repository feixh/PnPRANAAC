//
// Created by feixh on 8/29/15.
//

#pragma once
// STL
#include <iostream>
#include <fstream>
#include <vector>

// theia
#include <theia/solvers/ransac.h>
#include <theia/solvers/arrsac.h>
#include <theia/solvers/estimator.h>
#include <theia/solvers/random_sampler.h>
#include <theia/util/timer.h>


// eigen
#include <Eigen/Core>
#include <Eigen/StdVector>

// P3P
#include "pnpsolvers/P3P_Kneip.h"

namespace ransac_estimators
{
	using namespace theia;
	using namespace Eigen;
	using namespace std;
	
struct Match2D3D
{
    Eigen::Vector3d featureVector;  // unitary bearing vectors
    Eigen::Vector3d worldPoint;  // points in world coordinate system
};


class P3PEstimator : public Estimator< Match2D3D, Matrix<double, 3, 4 > > {
public:
    P3PEstimator():
        Estimator< Match2D3D, Matrix<double, 3, 4> >(),
        solver(){}

// Get the minimum number of samples needed to generate a model.
    virtual double SampleSize() const {
        return 3;
    }

    // Given a set of data points, estimate the model. Users should implement this
    // function appropriately for the task being solved. Returns true for
    // successful model estimation (and outputs model), false for failed
    // estimation. Typically, this is a minimal set, but it is not required to be.
    virtual bool EstimateModel(const std::vector<Datum> &data, std::vector<Model> *model) const {
        assert(data.size() >= 3);
        Matrix3d featureVectors;
        Matrix3d worldPoints;
        for (size_t i = 0; i < 3; ++i) {
            featureVectors.col(i) = data[i].featureVector;
            worldPoints.col(i)    = data[i].worldPoint;
        }
        int success = solver.computePoses(featureVectors, worldPoints, *model);
        for (auto it = model->begin(); it != model->end();) {
            if ( !it->allFinite() ) {
                it = model->erase(it);
            } else {
                ++it;
            }
        }
        if ( model->empty() ){
            success = -1;
        }
        if ( success == -1 ){
            return false;
        }else{
            return true;
        }
    }
  // Given a model and a data point, calculate the error. Users should implement
  // this function appropriately for the task being solved.
  virtual double Error(const Datum& data, const Model& model) const {
      // model is gwc
      const Vector3d &worldPoint( data.worldPoint );
      Vector3d proj( model.block<3,3>(0,0).transpose()*( worldPoint - model.block<3,1>(0,3) ) );
      if ( proj(2) < 0 ){
          return 1000000;
      }
      const Vector3d &featureVector( data.featureVector );
      double dx( featureVector(0)/featureVector(2) - proj(0)/proj(2) );
      double dy( featureVector(1)/featureVector(2) - proj(1)/proj(2) );
      return dx*dx + dy*dy;
  }
private:
    P3P_Kneip solver;
};

//     // setup ransac parameters
//     RansacParameters ransac_params;
//     ransac_params.error_thresh = 1e-2;
//     ransac_params.failure_probability = 0.05;
//     ransac_params.max_iterations = 3000;
//     ransac_params.min_inlier_ratio = 0.1;
//     ransac_params.use_mle = true;
// 
// 
//     P3PEstimator estimator;
//      Ransac< P3PEstimator > ransac(ransac_params, estimator);
// //    Arrsac< P3PEstimator > ransac( ransac_params, estimator, 500, 30);
//     ransac.Initialize();
//     RansacSummary summary;
//     Matrix< double, 3, 4 > best_model;
//     Timer tt;
//     ransac.Estimate( data, &best_model, &summary );
//     double duration( tt.ElapsedTimeInSeconds() );
//     cout << duration << " s" << endl;
// 
//     // checkout results
//     cout << best_model << endl;
//     cout << "iterations:" << summary.num_iterations << endl;
//     for ( size_t i = 0; i < summary.inliers.size(); ++i )
//     {
//         cout << summary.inliers[i] << " ";
//     }
//     cout << endl;

}