/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2009, Willow Garage, Inc.
 *  Copyright (c) 2012-, Open Perception, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id$
 *
 */

#ifndef OPENGV_SAC_PROSAC_HPP_
#define OPENGV_SAC_PROSAC_HPP_

#include <opengv/sac/SampleConsensus.hpp>
#include <boost/math/distributions/binomial.hpp>

namespace opengv {
    namespace sac {
        /** \brief @b RandomSampleConsensus represents an implementation of the RANSAC (RAndom SAmple Consensus) algorithm, as
          * described in: "Matching with PROSAC – Progressive Sample Consensus", Chum, O. and Matas, J.G., CVPR, I: 220-226
          * 2005.
          * \author Vincent Rabaud
          * \ingroup sample_consensus
          */
        template<typename PROBLEM_T>
        class ProgressiveSampleConsensus : public SampleConsensus<PROBLEM_T> {
        public:
            /** A child of SampleConsensusProblem */
            typedef PROBLEM_T problem_t;
            /** The model we trying to fit */
            typedef typename problem_t::model_t model_t;

//      using Ptr = boost::shared_ptr<ProgressiveSampleConsensus>;
//      using ConstPtr = boost::shared_ptr<const ProgressiveSampleConsensus>;

            using SampleConsensus<problem_t>::max_iterations_;
            using SampleConsensus<problem_t>::threshold_;
            using SampleConsensus<problem_t>::iterations_;
            using SampleConsensus<problem_t>::sac_model_;
            using SampleConsensus<problem_t>::model_;
            using SampleConsensus<problem_t>::model_coefficients_;
            using SampleConsensus<problem_t>::inliers_;
            using SampleConsensus<problem_t>::probability_;
            using SampleConsensus<problem_t>::max_time_;

            /**
             * \brief Constructor.
             */
            ProgressiveSampleConsensus(int maxIterations = 1000,
                                       double threshold = 1.0,
                                       double probability = 0.99,
                                       double maxTime = 1.0);

            /**
             * \brief Destructor.
             */
            virtual ~ProgressiveSampleConsensus();

            /**
             * \brief Fit the model.
             */
            bool
            computeModel(int debug_verbosity_level = 0) override;
        };
    }
}

#include "implementation/prosac.hpp"

#endif /* OPENGV_SAC_PROSAC_HPP_ */