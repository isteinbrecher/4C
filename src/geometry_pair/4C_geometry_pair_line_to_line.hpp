// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_LINE_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_LINE_HPP


#include "4C_config.hpp"

#include "4C_geometry_pair.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"

#include <memory>
#include <vector>


FOUR_C_NAMESPACE_OPEN

// Forward declarations.
namespace GeometryPair
{
  enum class ProjectionResult;

  template <typename ScalarType>
  class ProjectionPoint1DTo3D;

  class LineToSurfaceEvaluationData;
}  // namespace GeometryPair


namespace GeometryPair
{
  /**
   * \brief Project a point in space to the surface element.
   * @param point (in) Point in space.
   * @param element_data_surface (in) Degrees of freedom for the surface.
   * @param xi (in/out) Parameter coordinates in the volume. The given values are the start values
   * for the Newton iteration.
   * @param projection_result (out) Flag for the result of the projection.
   * @param normal_influence_direction (in) Threshold of normal influence direction for the surface.
   * @param min_one_iteration (in) Flag if at least one NR iteration should be performed, even if
   * the initial residual satisfies the convergence check.
   */
  template <typename ScalarType, typename LineA, typename LineB>
  void line_to_line_closest_point_projection(
      const ElementData<LineA, ScalarType>& element_data_line_a,
      const ElementData<LineB, ScalarType>& element_data_line_b, ScalarType& eta_a,
      ScalarType& eta_b, ProjectionResult& projection_result, const bool min_one_iteration = false);
}  // namespace GeometryPair

FOUR_C_NAMESPACE_CLOSE

#endif
