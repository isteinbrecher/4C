// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_LINE_EVALUATION_DATA_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_LINE_EVALUATION_DATA_HPP


#include "4C_config.hpp"

#include "4C_fem_general_utils_integration.hpp"
#include "4C_geometry_pair_evaluation_data_base.hpp"
#include "4C_geometry_pair_input.hpp"
#include "4C_geometry_pair_utility_classes.hpp"

#include <set>

FOUR_C_NAMESPACE_OPEN

namespace GeometryPair
{
  /**
   * \brief Class to manage inout parameters and evaluation data for line to line interactions.
   */
  class LineToLineEvaluationData : public GeometryEvaluationDataBase
  {
   public:
    /**
     * \brief Constructor (derived).
     */
    LineToLineEvaluationData();

    /**
     * \brief Setup averaged tangent data if needed.
     */
    void setup(const std::shared_ptr<const Core::FE::Discretization>& discret,
        const Core::Conditions::Condition* condition_line,
        const Core::Conditions::Condition* condition_other);

    /**
     * \brief Clear the evaluation data.
     */
    void clear() override;

    /**
     * \brief Reset the evaluation data.
     */
    void reset_tracker();

   private:
    //! Strategy to be used for contact search.
    GeometryPair::LineTo3DStrategy strategy_;

    //! Gauss rule for Gauss point projection method.
    Core::FE::GaussRule1D gauss_rule_;

    //! Number of integration points in the circumferential direction of the line cross section.
    unsigned int integration_points_circumference_;

    //! Gauss point projection tracking vector.
    std::map<int, std::vector<bool>> gauss_point_projection_tracker_;

    //! Number of points for segmentation search.
    unsigned int n_search_points_;

    //! What to do if not all Gauss points of a segment project valid
    GeometryPair::NotAllGaussPointsProjectValidAction not_all_gauss_points_project_valid_action_;

    //! Segment tracking vector for segmentation. We use double in this case, because otherwise the
    //! class would have to be templated on the type of this tracker.
    std::map<int, std::set<LineSegment<double>>> segment_tracker_;
  };
}  // namespace GeometryPair

FOUR_C_NAMESPACE_CLOSE

#endif
