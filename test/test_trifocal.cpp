/******************************************************************************
 * Author:   Laurent Kneip                                                    *
 * Contact:  kneip.laurent@gmail.com                                          *
 * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 * * Redistributions of source code must retain the above copyright           *
 *   notice, this list of conditions and the following disclaimer.            *
 * * Redistributions in binary form must reproduce the above copyright        *
 *   notice, this list of conditions and the following disclaimer in the      *
 *   documentation and/or other materials provided with the distribution.     *
 * * Neither the name of ANU nor the names of its contributors may be         *
 *   used to endorse or promote products derived from this software without   *
 *   specific prior written permission.                                       *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
 * SUCH DAMAGE.                                                               *
 ******************************************************************************/

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <climits>
#include <Eigen/Eigen>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac/Prosac.hpp>
#include <opengv/sac/Lmeds.hpp>
#include <opengv/sac_problems/relative_pose/CentralRelativePoseSacProblem.hpp>
#include <opengv/relative_pose/NoncentralRelativeMultiAdapter.hpp>
#include <opengv/sac/MultiRansac.hpp>
#include <opengv/sac_problems/relative_pose/MultiNoncentralRelativePoseSacProblem.hpp>
#include <opengv/trifocal/methods.hpp>
#include <opengv/trifocal/TrifocalSacProblem.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp>
#include <opengv/triangulation/methods.hpp>
#include <sstream>
#include <fstream>

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"


using namespace std;
using namespace Eigen;
using namespace opengv;

int main( int argc, char** argv )
{
  // initialize random seed
  initializeRandomSeed();

  //set experiment parameters
  double noise = 0.5;
  double outlierFraction = 0.3;
  size_t numberPoints = 100;

  //generate a random pose for viewpoint 1
  translation_t position1 = Eigen::Vector3d::Zero();
  rotation_t rotation1 = Eigen::Matrix3d::Identity();

  //generate a random pose for viewpoint 2
  translation_t position2 = generateRandomTranslation(2.0);
  rotation_t rotation2 = generateRandomRotation(0.5);

  //generate a random pose for viewpoint 3
  translation_t position3 = generateRandomTranslation(2.0);
  rotation_t rotation3 = generateRandomRotation(0.5);

  //create a fake central camera
  translations_t camOffsets;
  rotations_t camRotations;
  generateCentralCameraSystem( camOffsets, camRotations );

  //derive correspondences based on random point-cloud
  bearingVectors_t bearingVectors1;
  bearingVectors_t bearingVectors2;
  bearingVectors_t bearingVectors3;
  Eigen::MatrixXd gt(3,numberPoints);
  generateRandom2D2D2DCorrespondences(
      position1, rotation1, position2, rotation2, position3, rotation3,
      camOffsets, camRotations, numberPoints, noise, outlierFraction,
      bearingVectors1, bearingVectors2, bearingVectors3, gt );

  //Extract the relative pose
  translation_t position12; rotation_t rotation12;
  extractRelativePose(
      position1, position2, rotation1, rotation2, position12, rotation12 );

  translation_t position13; rotation_t rotation13;
  extractRelativePose(
          position1, position3, rotation1, rotation3, position13, rotation13 );

  //print experiment characteristics
  printExperimentCharacteristics( position12, rotation12, noise, outlierFraction );

  printExperimentCharacteristics( position13, rotation13, noise, outlierFraction );


  //compute and print the essential-matrix
  printEssentialMatrix( position12, rotation12 );

  //create a central relative adapter
  relative_pose::CentralRelativeAdapter adapter1(
      bearingVectors1,
      bearingVectors2,
      rotation12);

  relative_pose::CentralRelativeAdapter adapter2(
      bearingVectors1,
      bearingVectors3,
      rotation13);

  //Create a RelativePoseSac problem and Ransac
  //Set algorithm to NISTER, STEWENIUS, SEVENPT, or EIGHTPT
  sac::Ransac<
          sac_problems::TrifocalSacProblem> ransac;
  std::shared_ptr<
          sac_problems::TrifocalSacProblem> relposeproblem_ptr(
          new sac_problems::TrifocalSacProblem(
                  adapter1, adapter2, opengv::sac_problems::TrifocalSacProblem::VectorR27,
                  opengv::sac_problems::TrifocalSacProblem::AlgebraicDist));
  ransac.sac_model_ = relposeproblem_ptr;
  ransac.threshold_ = 0.1;
  ransac.max_iterations_ = 100;


  struct timeval tic;
  struct timeval toc;
  gettimeofday( &tic, 0 );
  ransac.computeModel();
  gettimeofday( &toc, 0 );
  double ransac_time = TIMETODOUBLE(timeval_minus(toc,tic));


  std::cout << "the ransac threshold is: " << ransac.threshold_ << std::endl;
//  std::cout << "the ransac results is: " << std::endl;
//  std::cout << ransac.model_coefficients_ << std::endl << std::endl;
//  std::cout << "the normalized translation is: " << std::endl;
//  std::cout << ransac.model_coefficients_.col(3)/
//               ransac.model_coefficients_.col(3).norm() << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac.iterations_ << " iterations and ";
  std::cout << ransac_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac.inliers_.size(); i++)
    std::cout << ransac.inliers_[i] << " ";
  std::cout << std::endl << std::endl;



  //Create a RelativePoseSac problem and Ransac
  //Set algorithm to NISTER, STEWENIUS, SEVENPT, or EIGHTPT
  sac::Ransac<
          sac_problems::relative_pose::CentralRelativePoseSacProblem> ransac_epipolar;
  std::shared_ptr<
          sac_problems::relative_pose::CentralRelativePoseSacProblem> epipolarproblem_ptr(
          new sac_problems::relative_pose::CentralRelativePoseSacProblem(
                  adapter1,
                  sac_problems::relative_pose::CentralRelativePoseSacProblem::STEWENIUS));
  ransac_epipolar.sac_model_ = epipolarproblem_ptr;
  ransac_epipolar.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
  ransac_epipolar.max_iterations_ = 100;

  //Run the experiment
  gettimeofday( &tic, 0 );
  ransac_epipolar.computeModel();
  gettimeofday( &toc, 0 );
  double ransac_epipolar_time = TIMETODOUBLE(timeval_minus(toc,tic));

  //print the results
  std::cout << "the ransac results is: " << std::endl;
  std::cout << ransac_epipolar.model_coefficients_ << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac_epipolar.iterations_ << " iterations and ";
  std::cout << ransac_epipolar_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac_epipolar.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac_epipolar.inliers_.size(); i++)
    std::cout << ransac_epipolar.inliers_[i] << " ";
  std::cout << std::endl << std::endl;

  points_t points;
  bearingVectors_t bearingVectors3_matched;
  for(size_t i = 0; i < ransac_epipolar.inliers_.size(); i++)
  {
//    std::cout << ransac_epipolar.inliers_[i] << " ";
    points.push_back(opengv::triangulation::triangulate2(
    adapter1,
    ransac_epipolar.inliers_[i]));
    bearingVectors3_matched.push_back(bearingVectors3[ransac_epipolar.inliers_[i]]);
  }

  //create a central absolute adapter
  absolute_pose::CentralAbsoluteAdapter adapter_abs(
          bearingVectors3_matched,
          points,
          rotation1);

  //Create an AbsolutePoseSac problem and Ransac
  //The method can be set to KNEIP, GAO or EPNP
  sac::Ransac<sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac_abs;
  std::shared_ptr<
          sac_problems::absolute_pose::AbsolutePoseSacProblem> absposeproblem_ptr(
          new sac_problems::absolute_pose::AbsolutePoseSacProblem(
                  adapter_abs,
                  sac_problems::absolute_pose::AbsolutePoseSacProblem::KNEIP));
  ransac_abs.sac_model_ = absposeproblem_ptr;
  ransac_abs.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
  ransac_abs.max_iterations_ = 100;

  gettimeofday( &tic, 0 );
  ransac_abs.computeModel();
  gettimeofday( &toc, 0 );
  double ransac_abs_time = TIMETODOUBLE(timeval_minus(toc,tic));

  //print the results
  std::cout << "the ransac results is: " << std::endl;
  std::cout << ransac_abs.model_coefficients_ << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac_abs.iterations_ << " iterations and ";
  std::cout << ransac_abs_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac_abs.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac_abs.inliers_.size(); i++)
    std::cout << ransac_abs.inliers_[i] << " ";
  std::cout << std::endl << std::endl;

}
