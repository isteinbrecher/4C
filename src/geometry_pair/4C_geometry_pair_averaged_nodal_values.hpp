// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_AVERAGED_NODAL_VALUES_HPP
#define FOUR_C_GEOMETRY_PAIR_AVERAGED_NODAL_VALUES_HPP


#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"

#include <set>
#include <unordered_map>


FOUR_C_NAMESPACE_OPEN

namespace GeometryPair
{
  struct AveragedNodalValuesData
  {
    struct A
    {
      std::set<int> connected_element_gid;
    };

    std::unordered_map<int, std::vector<int>> element_gid_to_element_dof;
    std::unordered_map<int, std::set<int>> element_gid_to_connected_patch_element_gid;
    std::unordered_map<int, std::vector<int>> element_gid_to_patch_dof;
    std::unordered_map<int, int> gid_to_count_map;
    std::unordered_map<int, double> gid_to_averaged_value_map;

    std::unordered_map<int, A> averaged_node_data;
  };

  AveragedNodalValuesData compute_averaged_nodal_values_id_data(
      const std::set<int>& element_ids_in_pairs,
      const std::unordered_map<int, const Core::Elements::Element*>& gid_to_evaluation_element_map);
}  // namespace GeometryPair

FOUR_C_NAMESPACE_CLOSE

#endif
