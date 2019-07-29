//
// Created by Richard Ding on 2019-07-13.
//

#ifndef OPENGV_TRIFOCAL_METHODS_HPP
#define OPENGV_TRIFOCAL_METHODS_HPP

#include <opengv/types.hpp>
#include <opengv/relative_pose/RelativeAdapterBase.hpp>
#include <opengv/Indices.hpp>


/**
 * \brief The namespace of this library.
 */
namespace opengv
{
/**
* \brief The namespace for the trifocal methods.
*/
namespace trifocal
{
trifocal_design_t computeDesignMatrix(
  const relative_pose::RelativeAdapterBase &adapter1,
  const relative_pose::RelativeAdapterBase &adapter2,
  const Indices &indices
);

trifocal_design_t computeDesignMatrix(
  const relative_pose::RelativeAdapterBase &adapter1,
  const relative_pose::RelativeAdapterBase &adapter2,
  const std::vector<int> & indices
);

trifocal_design_t computeDesignMatrix(
  const relative_pose::RelativeAdapterBase &adapter1,
  const relative_pose::RelativeAdapterBase &adapter2
);


template <typename derived>
void computeTrilinearEmbeddingsPPP(
        const bearingVector_t & P1,
        const bearingVector_t & P2,
        const bearingVector_t & P3,
        Eigen::MatrixBase<derived> &M
        );

trifocalTensor_t algebraicMinimization(
        trifocalTensor_t t,
        const trifocal_design_t & A);

double computeTransferError(const trifocalTensor_t & t,
                            const bearingVector_t & P1,
                            const bearingVector_t & P2,
                            const bearingVector_t & P3);

Eigen::VectorXd computeAlgebraicError(const trifocal_design_t & A,
                             const trifocalTensor_t & t);

bearingVector_t pointTransfer(const trifocalTensor_t & t,
                              const bearingVector_t & P1,
                              const bearingVector_t & P2);



//trifocalTensor_t transformTrifocalTensor(
//        trifocalTensor_t T,
//        Eigen::Matrix3d M1,
//        Eigen::Matrix3d M2,
//        Eigen::Matrix3d M3,
//        bool normalize
//        );





}
}
#endif //OPENGV_TRIFOCAL_METHODS_HPP
