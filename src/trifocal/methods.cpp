//
// Created by Tianjiao on 2019-07-24.
//

#include <opengv/trifocal/methods.hpp>
#include <Eigen/KroneckerProduct>
#include <iostream>


opengv::trifocal_design_t
opengv::trifocal::computeDesignMatrix(
  const relative_pose::RelativeAdapterBase &adapter1,
  const relative_pose::RelativeAdapterBase &adapter2,
  const Indices &indices)
{
  assert(adapter1.getNumberCorrespondences() == adapter2.getNumberCorrespondences());
  trifocal_design_t A(27, 4*indices.size());
  for (size_t i=0; i< indices.size(); i++){
    auto ASlice = A.block<27,4>(0, i*4);
    opengv::trifocal::computeTrilinearEmbeddingsPPP(
      adapter1.getBearingVector1(indices[i]),
      adapter1.getBearingVector2(indices[i]),
      adapter2.getBearingVector2(indices[i]),
      ASlice
    );
  }
  return A;
}

opengv::trifocal_design_t
opengv::trifocal::computeDesignMatrix(
  const relative_pose::RelativeAdapterBase &adapter1,
  const relative_pose::RelativeAdapterBase &adapter2,
  const std::vector<int> & indices)
{
  Indices idx(indices);
  return opengv::trifocal::computeDesignMatrix(adapter1, adapter2, idx);
}

opengv::trifocal_design_t
opengv::trifocal::computeDesignMatrix(
  const relative_pose::RelativeAdapterBase &adapter1,
  const relative_pose::RelativeAdapterBase &adapter2)
{
  Indices idx(adapter1.getNumberCorrespondences());
  return opengv::trifocal::computeDesignMatrix(adapter1, adapter2, idx);
}



template <typename derived>
void
opengv::trifocal::computeTrilinearEmbeddingsPPP(
        const bearingVector_t & P1,
        const bearingVector_t & P2,
        const bearingVector_t & P3,
        Eigen::MatrixBase<derived> &M)
{
  double x1=P1(0);
  double y1=P1(1);
  double z1=P1(2);
  double x2=P2(0);
  double y2=P2(1);
  double z2=P2(2);
  double x3=P3(0);
  double y3=P3(1);
  double z3=P3(2);

  M.col(0) << x1*z2*z3,0,-x1*x2*z3, 0,0,0, -x1*z2*x3,0,x1*x2*x3,
          y1*z2*z3,0,-x2*y1*z3,  0,0,0, -x3*y1*z2,0,x2*x3*y1,
          z1*z2*z3,0,-z1*x2*z3, 0,0,0, -z1*z2*x3,0,z1*x2*x3;

  M.col(1) << 0,x1*z2*z3,-x1*y2*z3, 0,0,0, 0,-x1*z2*x3,x1*x3*y2,
          0,y1*z2*z3,-y1*y2*z3, 0,0,0, 0,-x3*y1*z2,x3*y1*y2,
          0,z1*z2*z3,-z1*y2*z3, 0,0,0, 0,-z1*z2*x3,x3*y2*z1;

  M.col(2) << 0,0,0, x1*z2*z3,0,-x1*x2*z3, -x1*z2*y3,0,x1*x2*y3,
          0,0,0, y1*z2*z3,0,-x2*y1*z3, -y1*z2*y3,0,x2*y1*y3,
          0,0,0, z1*z2*z3,0,-z1*x2*z3, -z1*z2*y3,0,z1*x2*y3;

  M.col(3) << 0,0,0, 0,x1*z2*z3,-x1*y2*z3, 0,-x1*z2*y3,x1*y2*y3,
          0,0,0, 0,y1*z2*z3,-y1*y2*z3, 0,-y1*z2*y3,y1*y2*y3,
          0,0,0, 0,z1*z2*z3,-z1*y2*z3, 0,-z1*z2*y3,z1*y2*y3;
}

template <typename derivedMat>
Eigen::MatrixXd computeLeastSingularVectorSmall(const Eigen::MatrixBase<derivedMat> &M)
{
  Eigen::JacobiSVD< Eigen::MatrixXd > SVD(M, Eigen::ComputeFullV);
  return SVD.matrixV().rightCols<1>();
}

template <typename derivedMat>
Eigen::MatrixXd computeLeastSingularVectorLarge(const Eigen::MatrixBase<derivedMat> &M)
{
  Eigen::BDCSVD< Eigen::MatrixXd > SVD(M, Eigen::ComputeFullV);
  return SVD.matrixV().rightCols<1>();
}

opengv::trifocalTensor_t
opengv::trifocal::algebraicMinimization(opengv::trifocalTensor_t t,
        const opengv::trifocal_design_t &A)
{
  trifocalSlice_t T1, T2, T3;
  T1 = Eigen::Map<Matrix3d>(t.data());
  T2 = Eigen::Map<Matrix3d>(t.data()+9);
  T3 = Eigen::Map<Matrix3d>(t.data()+18);

  //Find the least singular vector
  Eigen::Matrix3d V;
  computeLeastSingularVectorSmall(T1);
  V.row(0) = computeLeastSingularVectorSmall(T1).transpose();
  V.row(1) = computeLeastSingularVectorSmall(T2).transpose();
  V.row(2) = computeLeastSingularVectorSmall(T3).transpose();
  bearingVector_t epi31 = computeLeastSingularVectorSmall(V);

  V.row(0) = computeLeastSingularVectorSmall(T1.transpose()).transpose();
  V.row(1) = computeLeastSingularVectorSmall(T2.transpose()).transpose();
  V.row(2) = computeLeastSingularVectorSmall(T3.transpose()).transpose();
  bearingVector_t epi21 = computeLeastSingularVectorSmall(V);

  Eigen::Matrix<double, 27, 9>  ELeft;
  ELeft = Eigen::kroneckerProduct(
          Eigen::Matrix3d::Identity(),
          Eigen::kroneckerProduct(
                  epi31,
                  Eigen::Matrix3d::Identity()
                  ));

  Eigen::Matrix<double, 27, 9>  ERight;
  ERight = (-1) * Eigen::kroneckerProduct(
          Eigen::Matrix<double, 9, 9>::Identity(),
          epi21);

  Eigen::Matrix<double, 27, 18> E;
  E << ELeft, ERight;

  Eigen::JacobiSVD< Eigen::MatrixXd > SVD(E,
          Eigen::ComputeThinU | Eigen::ComputeThinV);

  int rankE = SVD.rank();
  Eigen::MatrixXd Up = SVD.matrixU().leftCols(rankE);
  Eigen::MatrixXd Vp = SVD.matrixV().leftCols(rankE);
  Eigen::VectorXd Sp = SVD.singularValues().head(rankE);

  Eigen::MatrixXd AUp = A.transpose()*Up;
  Eigen::VectorXd tp = computeLeastSingularVectorLarge(AUp);
  trifocalTensor_t t_corrected = Up * tp;

  return t_corrected;
}

opengv::bearingVector_t
opengv::trifocal::pointTransfer(const opengv::trifocalTensor_t &t,
                                const opengv::bearingVector_t &P1,
                                const opengv::bearingVector_t &P2) {

  bearingVector_t P3 = bearingVector_t::Zero();
  size_t i=0;
  size_t j=1;

  for (size_t l=0; l<3; l++)
  {
    double v1 = 0;
    double v2 = 0;
    for (size_t k=0; k<3; k++)
    {
      v1 += P1[k] * t[j+3*l+9*k];
      v2 += P1[k] * t[i+3*l+9*k];
    }
    P3[l] = P2[i] * v1 - P2[j] * v2;
  }
  return P3;
}

Eigen::VectorXd
opengv::trifocal::computeAlgebraicError(const opengv::trifocal_design_t &A,
        const opengv::trifocalTensor_t &t) {
  size_t N = A.cols();
  Eigen::MatrixXd temp = (t.transpose() * A).array().square();
  temp.resize(4, N / 4);
  return temp.colwise().sum().cwiseSqrt();
}

double opengv::trifocal::computeTransferError(const opengv::trifocalTensor_t &t,
                                              const opengv::bearingVector_t &P1,
                                              const opengv::bearingVector_t &P2,
                                              const opengv::bearingVector_t &P3) {
  bearingVector_t P3_hat = pointTransfer(t, P1, P2);
  bearingVector_t P3_ = P3 / P3(2);
  P3_hat = P3_hat / P3_hat(2);
  return (P3_hat - P3_).norm();
}
