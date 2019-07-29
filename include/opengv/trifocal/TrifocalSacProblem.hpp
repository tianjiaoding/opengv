//
// Created by Tianjiao on 2019-07-28.
//

#ifndef OPENGV_TRIFOCALSACPROBLEM_HPP
#define OPENGV_TRIFOCALSACPROBLEM_HPP

#include <opengv/sac/SampleConsensusProblem.hpp>
#include <opengv/types.hpp>
#include <opengv/relative_pose/RelativeAdapterBase.hpp>
#include <opengv/trifocal/methods.hpp>

namespace opengv
{
namespace sac_problems
{
class TrifocalSacProblem :
  public sac::SampleConsensusProblem<trifocalModel_t>
{
  typedef trifocalModel_t model_t;
  typedef opengv::relative_pose::RelativeAdapterBase adapter_t;

  typedef enum TargetModel
  {
    VectorR27 = 0,
    VectorT = 1,
    RelativePoses = 2
  } targetModel_t;

  typedef enum ErrorMetric
  {
    AlgebraicDist = 0,
    EpipolarDist = 1,
    TransferDist = 2,
    ReprojectionDist = 3
  } errorMetric_t;

  TrifocalSacProblem(adapter_t & adapter1,
                     adapter_t & adapter2,
                     targetModel_t targetModel,
                     ErrorMetric errorMetric,
                     bool randomSeed = true) :
    sac::SampleConsensusProblem<model_t> (randomSeed),
    _adapter1(adapter1),
    _adapter2(adapter2),
    _targetModel(targetModel),
    _errorMetric(errorMetric)
  {
    setUniformIndices(adapter1.getNumberCorrespondences());
    _designMatrix = opengv::trifocal::computeDesignMatrix(
            _adapter1, _adapter2);
  };

  TrifocalSacProblem(adapter_t & adapter1,
                     adapter_t & adapter2,
                     targetModel_t targetModel,
                     ErrorMetric errorMetric,
                     const std::vector<int> & indices,
                     bool randomSeed = true) :
    sac::SampleConsensusProblem<model_t> (randomSeed),
    _adapter1(adapter1),
    _adapter2(adapter2),
    _targetModel(targetModel),
    _errorMetric(errorMetric)
  {
    setIndices(indices);
    _designMatrix = opengv::trifocal::computeDesignMatrix(
            _adapter1, _adapter2);
  };

  virtual ~TrifocalSacProblem() {};

  /**
   * \brief See parent-class.
   */
  virtual bool computeModelCoefficients(
    const std::vector<int> & indices,
    model_t & outModel) const;

  /**
   * \brief See parent-class.
   */
  virtual void getSelectedDistancesToModel(
    const model_t & model,
    const std::vector<int> & indices,
    std::vector<double> & scores) const;

  /**
   * \brief See parent-class.
   */
  virtual void optimizeModelCoefficients(
    const std::vector<int> & inliers,
    const model_t & model,
    model_t & optimized_model);

  /**
   * \brief See parent-class.
   */
  virtual int getSampleSize() const;

protected:
  adapter_t & _adapter1;
  adapter_t & _adapter2;
  targetModel_t _targetModel;
  errorMetric_t _errorMetric;
  trifocal_design_t _designMatrix;
};
}
}

#endif //OPENGV_TRIFOCALSACPROBLEM_HPP
