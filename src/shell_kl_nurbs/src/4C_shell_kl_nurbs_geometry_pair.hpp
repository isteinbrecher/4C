// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SHELL_KL_NURBS_GEOMETRY_PAIR_HPP
#define FOUR_C_SHELL_KL_NURBS_GEOMETRY_PAIR_HPP

#include "4C_config.hpp"

#include "4C_geometry_pair_element.hpp"
#include "4C_mat_shell_kl.hpp"
#include "4C_shell_kl_nurbs.hpp"


FOUR_C_NAMESPACE_OPEN

namespace GeometryPair
{
  /**
   * \brief For shells we need additional information for the geometry pairs, this is stored in this
   * struct.
   */
  template <typename ScalarType>
  struct OptionalElementData<t_nurbs9, ScalarType>
  {
    bool is_shell_;
    double shell_thickness_;
  };

  /**
   * \brief Initialize the optional element data for the shell.
   */
  template <typename ScalarType>
  struct InitializeOptionalElementData<t_nurbs9, ScalarType>
  {
    static void initialize(GeometryPair::ElementData<t_nurbs9, ScalarType> element_data,
        const Core::Elements::Element* element)
    {
      exit(0);
      const auto* kl_shell =
          dynamic_cast<const Discret::Elements::KirchhoffLoveShellNurbs*>(element);
      if (kl_shell != nullptr)
      {
        element_data.is_shell_ = true;
        element_data.shell_thickness_ =
            std::dynamic_pointer_cast<const Mat::KirchhoffLoveShell>(kl_shell->material())
                ->thickness();
      }
      else
      {
        element_data.is_shell_ = false;
        element_data.shell_thickness_ = 0.0;
      };
    }
  };

  /**
   * \brief Convert the additional data to double - this is a plain copy, as the original data is
   * already double.
   */
  template <typename ScalarType>
  struct OptionalElementDataToDouble<t_nurbs9, ScalarType>
  {
    static void to_double(GeometryPair::ElementData<t_nurbs9, double>& element_data_double,
        const GeometryPair::ElementData<t_nurbs9, ScalarType>& element_data)
    {
      element_data_double.is_shell_ = element_data.is_shell_;
      element_data_double.shell_thickness_ = element_data.shell_thickness_;
    }
  };

}  // namespace GeometryPair
FOUR_C_NAMESPACE_CLOSE

#endif
