// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_geometric_search_visualization.hpp"

#include "4C_fem_geometric_search_utils.hpp"
#include "4C_io_visualization_utils.hpp"

#include <utility>

#include "ArborX_Box.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::GeometricSearch
{
  GeometricSearchVisualization::GeometricSearchVisualization(
      Core::IO::VisualizationParameters parameters, const Epetra_Comm& comm,
      std::string base_output_name)
      : Core::IO::VisualizationManager(std::move(parameters), comm, std::move(base_output_name))
  {
    my_rank_ = comm.MyPID();
    this->get_visualization_data().register_cell_data<int>("element_id", 1);
    this->get_visualization_data().register_cell_data<int>("element_created_on_rank", 1);
    this->get_visualization_data().register_cell_data<int>("primitive_predicate_flag", 1);
  }

  void GeometricSearchVisualization::write_primitives_and_predicates_to_disk(
      const double visualziation_time, const int visualization_step,
      const std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>>& primitives,
      const std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>>& predicates)
  {
#ifndef FOUR_C_WITH_ARBORX
    FOUR_C_THROW(
        "Core::GeometricSearch::write_primitives_and_predicates_to_disk can only be used with "
        "ArborX."
        "To use it, enable ArborX during the configure process.");
#else
    clear_data();

    auto& visualization_data = get_visualization_data();

    // Function that adds bounding volumes to the visualization data
    const auto add_bounding_volumes =
        [&](const auto& bounding_volumes, const int primitive_predicate_flag)
    {
      for (const auto& bounding_volume : bounding_volumes)
      {
        // const auto& [points, polygons] =
        //     get_k_dop_polyhedron_representation(bounding_volume.second);
        // IO::append_polyhedron_to_visualization_data(visualization_data, points, polygons);

        const auto box = ArborX::Box(bounding_volume.second.bounding_volume_);
        const auto min_corner = box.minCorner();
        const auto max_corner = box.maxCorner();

        for (int i_dim = 0; i_dim < 3; i_dim++)
        {
          visualization_data.get_point_coordinates().push_back(min_corner[i_dim]);
        }

        visualization_data.get_point_coordinates().push_back(max_corner[0]);
        visualization_data.get_point_coordinates().push_back(min_corner[1]);
        visualization_data.get_point_coordinates().push_back(min_corner[2]);

        visualization_data.get_point_coordinates().push_back(max_corner[0]);
        visualization_data.get_point_coordinates().push_back(max_corner[1]);
        visualization_data.get_point_coordinates().push_back(min_corner[2]);

        visualization_data.get_point_coordinates().push_back(min_corner[0]);
        visualization_data.get_point_coordinates().push_back(max_corner[1]);
        visualization_data.get_point_coordinates().push_back(min_corner[2]);


        visualization_data.get_point_coordinates().push_back(min_corner[0]);
        visualization_data.get_point_coordinates().push_back(min_corner[1]);
        visualization_data.get_point_coordinates().push_back(max_corner[2]);

        visualization_data.get_point_coordinates().push_back(max_corner[0]);
        visualization_data.get_point_coordinates().push_back(min_corner[1]);
        visualization_data.get_point_coordinates().push_back(max_corner[2]);

        for (int i_dim = 0; i_dim < 3; i_dim++)
        {
          visualization_data.get_point_coordinates().push_back(max_corner[i_dim]);
        }

        visualization_data.get_point_coordinates().push_back(min_corner[0]);
        visualization_data.get_point_coordinates().push_back(max_corner[1]);
        visualization_data.get_point_coordinates().push_back(max_corner[2]);



        visualization_data.get_cell_types().push_back(12);
        auto& cell_offsets = visualization_data.get_cell_offsets();
        cell_offsets.push_back(visualization_data.get_point_coordinates_number_of_points());

        visualization_data.get_cell_data<int>("element_id").push_back(bounding_volume.first);
        visualization_data.get_cell_data<int>("element_created_on_rank").push_back(my_rank_);
        visualization_data.get_cell_data<int>("primitive_predicate_flag")
            .push_back(primitive_predicate_flag);
      }
    };

    // All bounding volumes to the output data
    add_bounding_volumes(primitives, 0);
    add_bounding_volumes(predicates, 1);

    // Write the data to disk
    write_to_disk(visualziation_time, visualization_step);
#endif
  }

}  // namespace Core::GeometricSearch

FOUR_C_NAMESPACE_CLOSE
