//
// Created by feixh on 8/29/15.
//

// STL
#include <iostream>

// theia
#include <theia/solvers/ransac.h>
#include <theia/solvers/estimator.h>

// eigen
#include <Eigen/Core>
#include <Eigen/StdVector>

using namespace theia;
using namespace std;
using namespace Eigen;

class PnPEstimator : public Estimator< float, float >
{
public:
// Get the minimum number of samples needed to generate a model.
  virtual double SampleSize() const{
        return 3;
    }

  // Given a set of data points, estimate the model. Users should implement this
  // function appropriately for the task being solved. Returns true for
  // successful model estimation (and outputs model), false for failed
  // estimation. Typically, this is a minimal set, but it is not required to be.
  virtual bool EstimateModel(const std::vector<Datum>& data,
                             std::vector<Model>* model) const {
      return true;

  }
  // Given a model and a data point, calculate the error. Users should implement
  // this function appropriately for the task being solved.
  virtual double Error(const Datum& data, const Model& model) const {
      return 0.;
  }
};
int main()
{
    RansacParameters ransac_params;
    PnPEstimator estimator;
    Ransac< PnPEstimator > ransac(ransac_params, estimator);
}
