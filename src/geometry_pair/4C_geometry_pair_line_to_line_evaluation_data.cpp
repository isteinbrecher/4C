// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometry_pair_line_to_line_evaluation_data.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
GeometryPair::LineToLineEvaluationData::LineToLineEvaluationData() : GeometryEvaluationDataBase()
{
  // Initialize evaluation data structures.
  clear();
}

/**
 *
 */
void GeometryPair::LineToLineEvaluationData::setup(
    const std::shared_ptr<const Core::FE::Discretization>& discret,
    const Core::Conditions::Condition* condition_line,
    const Core::Conditions::Condition* condition_other)
{
  exit(0);
  // // Setup the Gauss rule for integration.
  // gauss_rule_.setup(discret);

  // // Setup the number of integration points.
  // integration_points_circumference_ = discret->get_integration_points_circumference();

  // // Setup the segmentation search.
  // n_search_points_ = discret->get_n_search_points();

  // // Setup the not all Gauss points project valid action.
  // not_all_gauss_points_project_valid_action_ =
  //     discret->get_not_all_gauss_points_project_valid_action();
}

/**
 *
 */
void GeometryPair::LineToLineEvaluationData::clear()
{
  // Call reset on the base method.
  GeometryEvaluationDataBase::clear();

  // Initialize evaluation data structures.
  {
    // Tracker for gauss point projection method.
    gauss_point_projection_tracker_.clear();

    // Segment tracker for segmentation.
    segment_tracker_.clear();
  }
}

/**
 *
 */
void GeometryPair::LineToLineEvaluationData::reset_tracker()
{
  for (auto& data : gauss_point_projection_tracker_)
    std::fill(data.second.begin(), data.second.end(), false);

  for (auto& data : segment_tracker_) data.second.clear();
}

FOUR_C_NAMESPACE_CLOSE
