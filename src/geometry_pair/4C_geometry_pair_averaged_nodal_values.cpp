// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometry_pair_averaged_nodal_values.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_utils_exceptions.hpp"

#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
void GeometryPair::compute_averaged_nodal_values_id_data(const std::set<int>& element_ids_in_pairs,
    const std::unordered_map<int, const Core::Elements::Element*>& gid_to_evaluation_element_map)
{
  // Get all elements that are needed for evaluation.
  std::set<int> all_elements = element_ids_in_pairs;
  for (const auto& element_in_pair_gid : element_ids_in_pairs)
  {
    const auto find_in_condition = gid_to_evaluation_element_map.find(element_in_pair_gid);
    if (find_in_condition != gid_to_evaluation_element_map.end())
    {
      // Loop over the nodes of the element.
      const Core::Nodes::Node* const* nodes = find_in_condition->second->nodes();
      for (int i_node = 0; i_node < find_in_condition->second->num_node(); i_node++)
      {
        // Loop over the elements connected to that node and check if they are in this condition.
        for (auto ele : nodes[i_node]->adjacent_elements())
        {
          const int element_id = ele.global_id();
          const auto find_adjacent_element_in_condition =
              gid_to_evaluation_element_map.find(element_id);
          if (find_adjacent_element_in_condition != gid_to_evaluation_element_map.end())
          {
            // The element exists in this condition, so lets add it to the set of all elements.
            all_elements.insert(element_id);
          }
          else
          {
            // The element is not part of this condition, i.e. it will not be used for the
            // calculation of averaged normals. This allows for 'sharp' corners.
          }
        }
      }
    }
    else
    {
      FOUR_C_THROW(
          "Could not find element with global id {} in the given gid_to_evaluation_element_map "
          "map.",
          element_in_pair_gid);
    }
  }

  // Get the GIDs of all elements needed for evaluation.
  std::unordered_map<int, std::vector<int>> element_gid_to_dof_map;
  for (const auto& element_id : all_elements)
  {
    const auto* element = gid_to_evaluation_element_map.at(element_id);
    std::vector<int> dummy_1, dummy_2;
    element->location_vector(
        *element->discretization(), element_gid_to_dof_map[element_id], dummy_1, dummy_2);
  }


  // Print the size of element_gid_to_dof_map
  std::cout << "Size of element_gid_to_dof_map: " << element_gid_to_dof_map.size() << std::endl;
}


FOUR_C_NAMESPACE_CLOSE
