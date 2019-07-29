//
// Created by Tianjiao on 2019-07-28.
//

#include <opengv/trifocal/TrifocalSacProblem.hpp>
#include <opengv/trifocal/methods.hpp>

bool opengv::sac_problems::
TrifocalSacProblem::computeModelCoefficients(
  const std::vector<int> &indices,
  opengv::sac_problems::TrifocalSacProblem::model_t &outModel) const
{
  //Construct the design matrix for selected points.
  opengv::trifocal_design_t A = opengv::trifocal::computeDesignMatrix(
    _adapter1, _adapter2, indices);

  //Find the least singular vector
  Eigen::JacobiSVD< Eigen::MatrixXd > SVD(A, Eigen::ComputeFullU);
  Eigen::Matrix<double,27,1> t = SVD.matrixU().col(27);

  outModel.t = t;

  if (_targetModel == targetModel_t::VectorT)
  {
    fprintf(stdout,"Model type: VectorT. Exiting from computeModelCoefficients.");
    return true;
  }

  Eigen::Matrix<double,27,1> t_corrected =
          opengv::trifocal::algebraicMinimization(t, A);
  outModel.t_corrected = t_corrected;

  fprintf(stdout,"Model more than corrected t is not implemented!");

  return true;
}

void opengv::sac_problems::
TrifocalSacProblem::getSelectedDistancesToModel(
  const opengv::sac_problems::TrifocalSacProblem::model_t &model,
  const std::vector<int> &indices,
  std::vector<double> &scores) const
{
  switch (_errorMetric)
  {
    case AlgebraicDist:
    {
      Eigen::VectorXd errors = opengv::trifocal::computeAlgebraicError(
              _designMatrix, model.t);
      for (size_t i=0; i<scores.size(); i++)
      {
        scores[i] = errors[i];
      }
      break;
    }
    case TransferDist:
    {
      
    }

  }
}

void opengv::sac_problems::
TrifocalSacProblem::optimizeModelCoefficients(const std::vector<int> &inliers,
  const opengv::sac_problems::TrifocalSacProblem::model_t &model,
  opengv::sac_problems::TrifocalSacProblem::model_t &optimized_model)
{

}

int opengv::sac_problems::TrifocalSacProblem::getSampleSize() const
{
  return 0;
}
