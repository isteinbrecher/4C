// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_AVERAGED_NODAL_VALUES_HPP
#define FOUR_C_GEOMETRY_PAIR_AVERAGED_NODAL_VALUES_HPP


#include "4C_config.hpp"

#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"

#include <map>


FOUR_C_NAMESPACE_OPEN


namespace GeometryPair
{

  /**
   *
   */
  template <typename Element>
  void setup(const std::shared_ptr<const Core::FE::Discretization>& discretization,
      const Core::Conditions::Condition* condition)
  {
    std::map<int, std::vector<int>> element_id_to_element_gids;
    for (const auto& [element_id, element] : condition->geometry())
    {
      element_id_to_element_gids[element->id()] = {1, 2, 3, 4};
    }


    // // Loop over all pairs and add the needed face elements.
    // std::unordered_map<int, std::shared_ptr<GeometryPair::FaceElement>> pair_face_elements;
    // pair_face_elements.clear();
    // for (const auto& pair : condition_contact_pairs_)
    // {
    //   const int solid_id = pair->element2()->id();
    //   auto find_in_condition = surface_ids_.find(solid_id);
    //   if (find_in_condition != surface_ids_.end())
    //   {
    //     // Check if the face is already in the pair_face_elements map.
    //     auto find_in_pair = pair_face_elements.find(solid_id);
    //     if (find_in_pair == pair_face_elements.end())
    //     {
    //       // The face element has to be created and added to the contact pair.
    //       std::shared_ptr<GeometryPair::FaceElement> new_face_element =
    //           GeometryPair::face_element_factory(find_in_condition->second, fad_order,
    //               line_to_surface_evaluation_data->get_surface_normal_strategy());
    //       new_face_element->set_part_of_pair(true);
    //       pair_face_elements[solid_id] = new_face_element;
    //       pair->set_face_element(new_face_element);
    //     }
    //     else
    //     {
    //       // Add the existing face element to the contact pair.
    //       pair->set_face_element(find_in_pair->second);
    //     }
    //   }
    //   else
    //   {
    //     FOUR_C_THROW("The face of the solid element {} is not in the current condition!",
    //         pair->element2()->id());
    //   }
    // }

    // // Now all faces of contact pairs are in pair_face_elements, we still need to add faces that
    // are
    // // needed for averaged normal calculation, but are not contained in any pair.
    // std::unordered_map<int, std::shared_ptr<GeometryPair::FaceElement>> face_elements_needed;
    // face_elements_needed = pair_face_elements;
    // for (const auto& face_element_iterator : pair_face_elements)
    // {
    //   // Loop over the nodes of the face element.
    //   const Core::Nodes::Node* const* nodes =
    //   face_element_iterator.second->get_element()->nodes(); for (int i_node = 0; i_node <
    //   face_element_iterator.second->get_element()->num_node();
    //       i_node++)
    //   {
    //     // Loop over the elements connected to that node and check if they are in this condition.
    //     for (auto ele : nodes[i_node]->adjacent_elements())
    //     {
    //       const int element_id = ele.global_id();
    //       auto find_in_condition = surface_ids_.find(element_id);
    //       if (find_in_condition != surface_ids_.end())
    //       {
    //         // The element exists in this condition, check if it is already in the needed faces
    //         map. auto find_in_needed = face_elements_needed.find(element_id); if (find_in_needed
    //         == face_elements_needed.end())
    //         {
    //           // It is not already in the needed faces -> add it.
    //           face_elements_needed[element_id] =
    //               GeometryPair::face_element_factory(find_in_condition->second, fad_order,
    //                   line_to_surface_evaluation_data->get_surface_normal_strategy());
    //         }
    //       }
    //       else
    //       {
    //         // The element is not part of this condition, i.e. it will not be used for the
    //         // calculation of averaged normals. This allows for 'sharp' corners.
    //       }
    //     }
    //   }
    // }

    // // Setup the geometry data for the surface patch.
    // line_to_surface_evaluation_data->setup(discret, face_elements_needed);
  }

}  // namespace GeometryPair


FOUR_C_NAMESPACE_CLOSE

#endif
