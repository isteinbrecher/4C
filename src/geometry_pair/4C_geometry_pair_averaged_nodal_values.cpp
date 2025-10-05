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
GeometryPair::AveragedNodalValuesData GeometryPair::compute_averaged_nodal_values_id_data(
    const std::set<int>& element_ids_in_pairs,
    const std::unordered_map<int, const Core::Elements::Element*>& gid_to_evaluation_element_map)
{
  AveragedNodalValuesData averaged_nodal_data{};

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
        const Core::Nodes::Node* node = nodes[i_node];
        averaged_nodal_data.averaged_node_data[node->id()].connected_element_gid.insert(
            element_in_pair_gid);

        // Loop over the elements connected to that node and check if they are in the condition.
        for (auto ele : node->adjacent_elements())
        {
          const int element_id = ele.global_id();
          const auto find_adjacent_element_in_condition =
              gid_to_evaluation_element_map.find(element_id);
          if (find_adjacent_element_in_condition != gid_to_evaluation_element_map.end())
          {
            // The element exists in this condition, so lets add it to the set of all elements.
            all_elements.insert(element_id);

            // Add the element to the node since it is connected to the node.
            averaged_nodal_data.averaged_node_data[node->id()].connected_element_gid.insert(
                element_id);

            // Also add the element ID to the connected patch element IDs for the current element.
            if (element_in_pair_gid != element_id)
              averaged_nodal_data.element_gid_to_connected_patch_element_gid[element_in_pair_gid]
                  .insert(element_id);
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
  for (const auto& element_id : all_elements)
  {
    const auto* element = gid_to_evaluation_element_map.at(element_id);
    std::vector<int> dummy_1, dummy_2;
    element->location_vector(*element->discretization(),
        averaged_nodal_data.element_gid_to_element_dof[element_id], dummy_1, dummy_2);
  }

  // Get the GIDs of each patch.
  for (const auto& element_id : element_ids_in_pairs)
  {
    auto& this_element_patch_dof = averaged_nodal_data.element_gid_to_patch_dof[element_id];
    this_element_patch_dof = averaged_nodal_data.element_gid_to_element_dof[element_id];
    for (const auto& connected_element_id :
        averaged_nodal_data.element_gid_to_connected_patch_element_gid[element_id])
    {
      for (const auto& item : averaged_nodal_data.element_gid_to_element_dof[connected_element_id])
      {
        if (std::find(this_element_patch_dof.begin(), this_element_patch_dof.end(), item) ==
            this_element_patch_dof.end())
          this_element_patch_dof.push_back(item);
      }
    }
  }



  // Print full contents of averaged_nodal_data for debugging
  std::cout << "AveragedNodalValuesData contents:\n";

  // Print averaged_node_data
  std::cout << "averaged_node_data:\n";
  for (const auto& [node_id, a_struct] : averaged_nodal_data.averaged_node_data)
  {
    std::cout << "  Node ID: " << node_id << ", connected_element_gid   : {";
    for (const auto& gid : a_struct.connected_element_gid)
    {
      std::cout << gid << ", ";
    }
    std::cout << "}\n";
  }

  // Print element_gid_to_connected_patch_element_gid
  std::cout << "element_gid_to_connected_patch_element_gid:\n";
  for (const auto& [element_gid, connected_gids] :
      averaged_nodal_data.element_gid_to_connected_patch_element_gid)
  {
    std::cout << "  Element GID: " << element_gid << ", connected_patch_element_gid: {";
    size_t count = 0;
    for (const auto& gid : connected_gids)
    {
      std::cout << gid;
      if (++count < connected_gids.size()) std::cout << ", ";
    }
    std::cout << "}\n";
  }

  // Print element_gid_to_element_dof
  std::cout << "element_gid_to_element_dof:\n";
  for (const auto& [element_gid, dofs] : averaged_nodal_data.element_gid_to_element_dof)
  {
    std::cout << "  Element GID: " << element_gid << ", dofs: [";
    for (size_t i = 0; i < dofs.size(); ++i)
    {
      std::cout << dofs[i];
      if (i + 1 < dofs.size()) std::cout << ", ";
    }
    std::cout << "]\n";
  }

  // Print element_gid_to_patch_dof
  std::cout << "element_gid_to_patch_dof:\n";
  for (const auto& [element_gid, patch_dofs] : averaged_nodal_data.element_gid_to_patch_dof)
  {
    std::cout << "  Element GID: " << element_gid << ", patch_dofs: [";
    for (size_t i = 0; i < patch_dofs.size(); ++i)
    {
      std::cout << patch_dofs[i];
      if (i + 1 < patch_dofs.size()) std::cout << ", ";
    }
    std::cout << "]\n";
  }

  return averaged_nodal_data;
}


FOUR_C_NAMESPACE_CLOSE
