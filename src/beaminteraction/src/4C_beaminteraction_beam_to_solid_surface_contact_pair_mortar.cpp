// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_beam_to_solid_surface_contact_pair_mortar.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_beam_to_solid_mortar_manager_contact.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_contact_params.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "4C_beaminteraction_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_element_faces.hpp"
#include "4C_geometry_pair_element_shape_functions.hpp"
#include "4C_geometry_pair_factory.hpp"
#include "4C_geometry_pair_line_to_line.hpp"
#include "4C_geometry_pair_line_to_surface.hpp"
#include "4C_geometry_pair_scalar_types.hpp"
#include "4C_geometry_pair_utility_classes.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_fad.hpp"

#include <Epetra_FEVector.h>

#include <filesystem>


FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface, typename Mortar>
BeamInteraction::BeamToSolidSurfaceContactPairMortar<ScalarType, Beam, Surface,
    Mortar>::BeamToSolidSurfaceContactPairMortar()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface, typename Mortar>
void BeamInteraction::BeamToSolidSurfaceContactPairMortar<ScalarType, Beam, Surface,
    Mortar>::evaluate_and_assemble_mortar_contributions(const Core::FE::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager,
    Core::LinAlg::SparseMatrix& global_constraint_lin_beam,
    Core::LinAlg::SparseMatrix& global_constraint_lin_solid,
    Core::LinAlg::SparseMatrix& global_force_beam_lin_lambda,
    Core::LinAlg::SparseMatrix& global_force_solid_lin_lambda, Epetra_FEVector& global_constraint,
    Epetra_FEVector& global_kappa, Core::LinAlg::SparseMatrix& global_kappa_lin_beam,
    Core::LinAlg::SparseMatrix& global_kappa_lin_solid, Epetra_FEVector& global_lambda_active,
    const std::shared_ptr<const Core::LinAlg::Vector<double>>& displacement_vector)
{
  // Call Evaluate on the geometry Pair
  this->cast_geometry_pair()->evaluate(
      this->ele1pos_, this->face_element_->get_face_element_data(), this->line_to_3D_segments_);

  // If there are no intersection segments, no contact terms will be assembled
  const unsigned int n_segments = this->line_to_3D_segments_.size();
  if (n_segments == 0) return;


  if (n_segments > 1) FOUR_C_THROW("only one segment at the moment");

  // Check if the beam is at the end of a filament.
  std::array<bool, 2> is_end_node;
  for (unsigned int i_node = 0; i_node < 2; i_node++)
  {
    const Core::Nodes::Node* node = this->element1()->nodes()[i_node];
    is_end_node[i_node] = node->num_element() == 1;
  }
  if (is_end_node[0] and is_end_node[1]) FOUR_C_THROW("only one endpoint possible");


  // Pointer to the contact parameters and input parameters
  const auto contact_parameters = this->params()->beam_to_solid_surface_contact_params();
  const auto contact_defined_on =
      contact_parameters->get_beam_to_solid_surface_contact_mortar_defined_in();

  // Get beam cross-section diameter
  const auto beam_ptr = dynamic_cast<const Discret::Elements::Beam3Base*>(this->element1());
  const double beam_cross_section_radius =
      beam_ptr->get_circular_cross_section_radius_for_interactions();

  // Initialize variables for the shape function matrices
  Core::LinAlg::Matrix<mortar_trial::spatial_dim_, mortar_trial::n_dof_, ScalarType> N_lambda_trial;
  Core::LinAlg::Matrix<Mortar::spatial_dim_, Mortar::n_dof_, ScalarType> N_lambda;
  Core::LinAlg::Matrix<Beam::spatial_dim_, Beam::n_dof_, ScalarType> N_beam;
  Core::LinAlg::Matrix<Surface::spatial_dim_, Surface::n_dof_, ScalarType> N_surface;

  // Initialize variables for local vectors
  Core::LinAlg::Matrix<mortar_trial::n_dof_, 1, ScalarType> constraint_vector(
      Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<mortar_trial::n_dof_, 1, ScalarType> kappa(
      Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<3, Mortar::n_dof_, ScalarType> normal_times_lambda_shape(
      Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<Beam::n_dof_, Mortar::n_dof_, ScalarType>
      beam_shape_times_normal_times_lambda_shape_gp(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<Surface::n_dof_, Mortar::n_dof_, ScalarType>
      surface_shape_times_normal_times_lambda_shape_gp(Core::LinAlg::Initialization::zero);
  beam_shape_times_normal_times_lambda_shape_.put_scalar(0.0);
  surface_shape_times_normal_times_lambda_shape_.put_scalar(0.0);

  // Integrate over segments
  for (const auto& segment : this->line_to_3D_segments_)
  {
    // Gauss point loop
    for (const auto& projected_gauss_point : segment.get_projection_points())
    {
      // Get the projection coordinates
      const auto& xi = projected_gauss_point.get_xi();
      const auto& eta = projected_gauss_point.get_eta();

      // Get the current Gauss integration factor. This includes everything, e.g., Gauss weight,
      // segment Jacobian and beam Jacobian.
      const ScalarType gauss_factor = projected_gauss_point.get_gauss_weight() *  //
                                      0.5 * segment.get_segment_length() *
                                      get_jacobian_for_configuration(eta, contact_defined_on);

      // Evaluate the contact kinematics
      const auto [_1, _2, surface_normal, gap] =
          this->evaluate_contact_kinematics_at_projection_point(
              projected_gauss_point, beam_cross_section_radius);

      // Get the shape function matrices
      GeometryPair::evaluate_shape_function_matrix<mortar_trial>(N_lambda_trial, eta);
      GeometryPair::evaluate_shape_function_matrix<Mortar>(N_lambda, eta);
      GeometryPair::evaluate_shape_function_matrix<Beam>(
          N_beam, eta, this->ele1pos_.shape_function_data_);
      GeometryPair::evaluate_shape_function_matrix<Surface>(
          N_surface, xi, this->face_element_->get_face_element_data().shape_function_data_);

      // Modify shape functions if applicable
      if (is_end_node[0])
      {
        ScalarType tmp = N_lambda(0);
        N_lambda(0) = 0.0;
        N_lambda(1) += tmp;

        tmp = N_lambda_trial(0);
        N_lambda_trial(0) = 0.0;
        N_lambda_trial(1) += tmp;
      }
      else if (is_end_node[1])
      {
        ScalarType tmp = N_lambda(1);
        N_lambda(1) = 0.0;
        N_lambda(0) += tmp;

        tmp = N_lambda_trial(1);
        N_lambda_trial(1) = 0.0;
        N_lambda_trial(0) += tmp;
      }

      // Weighted gap
      constraint_vector.update_t(gauss_factor * gap, N_lambda_trial, 1.0);

      // Force contributions
      normal_times_lambda_shape.multiply(surface_normal, N_lambda);

      beam_shape_times_normal_times_lambda_shape_gp.multiply_tn(N_beam, normal_times_lambda_shape);
      beam_shape_times_normal_times_lambda_shape_gp.scale(1.0 * gauss_factor);
      beam_shape_times_normal_times_lambda_shape_ += beam_shape_times_normal_times_lambda_shape_gp;

      surface_shape_times_normal_times_lambda_shape_gp.multiply_tn(
          N_surface, normal_times_lambda_shape);
      surface_shape_times_normal_times_lambda_shape_gp.scale(-1.0 * gauss_factor);
      surface_shape_times_normal_times_lambda_shape_ +=
          surface_shape_times_normal_times_lambda_shape_gp;

      // Scaling vector
      Core::LinAlg::Matrix<mortar_trial::spatial_dim_, 1, double> ones(
          Core::LinAlg::Initialization::zero);
      ones.put_scalar(1.0);
      Core::LinAlg::Matrix<mortar_trial::n_dof_, 1, ScalarType> N_lambda_trial_flat(
          Core::LinAlg::Initialization::zero);
      N_lambda_trial_flat.multiply_tn(N_lambda_trial, ones);
      N_lambda_trial_flat.scale(gauss_factor);
      kappa += N_lambda_trial_flat;
    }
  }

  // GIDs of the pair and the force vector acting on the pair.
  const auto [beam_centerline_gid, patch_gid] =
      get_beam_to_surface_pair_gid<Beam>(discret, *this->element1(), *this->face_element_);

  // Get the Lagrange multiplier GIDs.
  const auto& [lambda_gid_pos, _2] = mortar_manager->location_vector(*this);

  // Assemble into the matrix in the beam row and lambda column
  for (unsigned int i_beam = 0; i_beam < Beam::n_dof_; i_beam++)
  {
    for (unsigned int i_lambda = 0; i_lambda < Mortar::n_dof_; i_lambda++)
    {
      const double value = Core::FADUtils::cast_to_double(
          beam_shape_times_normal_times_lambda_shape_(i_beam, i_lambda));
      global_force_beam_lin_lambda.fe_assemble(
          value, beam_centerline_gid(i_beam), lambda_gid_pos[i_lambda]);
    }
  }

  // Assemble into the matrix in the surface row and lambda column
  for (unsigned int i_surface = 0; i_surface < Surface::n_dof_; i_surface++)
  {
    for (unsigned int i_lambda = 0; i_lambda < Mortar::n_dof_; i_lambda++)
    {
      const double value = Core::FADUtils::cast_to_double(
          surface_shape_times_normal_times_lambda_shape_(i_surface, i_lambda));
      global_force_solid_lin_lambda.fe_assemble(
          value, patch_gid[i_surface], lambda_gid_pos[i_lambda]);
    }
  }

  // Assemble into the in the lambda row
  for (unsigned int i_lambda = 0; i_lambda < Mortar::n_dof_; i_lambda++)
  {
    // Assemble into the beam column
    for (unsigned int i_beam = 0; i_beam < Beam::n_dof_; i_beam++)
    {
      const double value = constraint_vector(i_lambda).dx(i_beam);
      global_constraint_lin_beam.fe_assemble(
          value, lambda_gid_pos[i_lambda], beam_centerline_gid(i_beam));

      const double value_kappa_linearization = kappa(i_lambda).dx(i_beam);
      global_kappa_lin_beam.fe_assemble(
          value_kappa_linearization, lambda_gid_pos[i_lambda], beam_centerline_gid(i_beam));
    }

    // Assemble into the solid column
    for (unsigned int i_patch = 0; i_patch < patch_gid.size(); i_patch++)
    {
      const double value = constraint_vector(i_lambda).dx(Beam::n_dof_ + i_patch);
      global_constraint_lin_solid.fe_assemble(value, lambda_gid_pos[i_lambda], patch_gid[i_patch]);

      const double value_kappa_linearization = kappa(i_lambda).dx(Beam::n_dof_ + i_patch);
      global_kappa_lin_solid.fe_assemble(
          value_kappa_linearization, lambda_gid_pos[i_lambda], patch_gid[i_patch]);
    }
  }

  // Assemble into global coupling vector
  const auto constraint_vector_double = Core::FADUtils::cast_to_double(constraint_vector);
  global_constraint.SumIntoGlobalValues(
      lambda_gid_pos.size(), lambda_gid_pos.data(), constraint_vector_double.data());

  // Assemble into global kappa vector
  auto kappa_double = Core::FADUtils::cast_to_double(kappa);
  global_kappa.SumIntoGlobalValues(
      lambda_gid_pos.size(), lambda_gid_pos.data(), kappa_double.data());
  kappa_double.put_scalar(1.0);


  // Modify shape functions if applicable
  if (is_end_node[0])
  {
    kappa_double(0) = 0.0;
  }
  else if (is_end_node[1])
  {
    kappa_double(1) = 0.0;
  }

  global_lambda_active.SumIntoGlobalValues(
      lambda_gid_pos.size(), lambda_gid_pos.data(), kappa_double.data());
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface, typename Mortar>
void BeamInteraction::BeamToSolidSurfaceContactPairMortar<ScalarType, Beam, Surface,
    Mortar>::evaluate_and_assemble(const Core::FE::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager,
    const std::shared_ptr<Epetra_FEVector>& force_vector,
    const std::shared_ptr<Core::LinAlg::SparseMatrix>& stiffness_matrix,
    const Core::LinAlg::Vector<double>& global_lambda,
    const Core::LinAlg::Vector<double>& displacement_vector)
{
  is_edge_to_edge_ = false;
  is_endpoint_ = false;
  {  // Add line penalty terms


    const auto beam_ptr = dynamic_cast<const Discret::Elements::Beam3Base*>(this->element1());
    const double beam_cross_section_radius =
        beam_ptr->get_circular_cross_section_radius_for_interactions();

    // for (int i_node = 0; i_node < 8; i_node++)
    // {
    //   std::cout << "\n node " << i_node << " posx " << this->element2()->nodes()[i_node]->x()[0]
    //             << " posy " << this->element2()->nodes()[i_node]->x()[1] << " posz "
    //             << this->element2()->nodes()[i_node]->x()[2];
    // }
    // std::cout << "\n\n";

    // Nodes 0 and 1 are the edge nodes
    using mein_type = typename Core::FADUtils::HigherOrderFadType<1,
        Sacado::ELRFad::SLFad<double, GeometryPair::t_line2::n_dof_ + Beam::n_dof_>>::type;
    GeometryPair::ElementData<Beam, mein_type> beam_position;
    beam_position.shape_function_data_.ref_length_ =
        this->ele1pos_.shape_function_data_.ref_length_;
    for (int i_dof = 0; i_dof < Beam::n_dof_; i_dof++)
    {
      beam_position.element_position_(i_dof) =
          Core::FADUtils::HigherOrderFadValue<ScalarType>::apply(
              GeometryPair::t_line2::n_dof_ + Beam::n_dof_, i_dof,
              Core::FADUtils::cast_to_double(this->ele1pos_.element_position_(i_dof)));
    }

    mein_type eta_a = 0.0;
    mein_type eta_b = 0.0;
    GeometryPair::ProjectionResult result;
    GeometryPair::line_to_line_closest_point_projection(
        beam_position, this->edge_position_, eta_a, eta_b, result);

    if (this->element2()->id() == 47 and
        result == GeometryPair::ProjectionResult::projection_found_valid)
    {
      // std::cout << "\n\n\n\n\n\nFOUND!!!!!!!!!!\n\n\n\n\n";

      Core::LinAlg::Matrix<3, 1, mein_type> pos_beam;
      GeometryPair::evaluate_position(eta_a, beam_position, pos_beam);
      Core::LinAlg::Matrix<3, 1, mein_type> pos_edge;
      GeometryPair::evaluate_position(eta_b, this->edge_position_, pos_edge);


      Core::LinAlg::Matrix<3, 1, mein_type> diff;
      diff = pos_beam;
      diff -= pos_edge;
      mein_type gap = diff.norm2() - beam_cross_section_radius;
      mein_type force = penalty_force(gap, *this->params()->beam_to_solid_surface_contact_params(),
          2000 * beam_cross_section_radius);

      Core::LinAlg::Matrix<3, 1, mein_type> force_vec = diff;
      force_vec.scale(1.0 / diff.norm2());
      force_vec.scale(force);

      std::cout << "\n Gap: " << Core::FADUtils::cast_to_double(gap);



      // Get the shape function matrices.
      Core::LinAlg::Matrix<Beam::n_dof_ + GeometryPair::t_line2::n_dof_, 1, mein_type>
          pair_force_vector;
      Core::LinAlg::Matrix<1, Beam::n_nodes_ * Beam::n_val_, mein_type> N_beam;
      Core::LinAlg::Matrix<1, GeometryPair::t_line2::n_nodes_ * GeometryPair::t_line2::n_val_,
          mein_type>
          N_edge;
      GeometryPair::EvaluateShapeFunction<Beam>::evaluate(
          N_beam, eta_a, beam_position.shape_function_data_);
      GeometryPair::EvaluateShapeFunction<GeometryPair::t_line2>::evaluate(
          N_edge, eta_b, this->edge_position_.shape_function_data_);

      // Calculate the variation of the gap function multiplied with the surface normal vector.
      for (unsigned int i_shape = 0; i_shape < N_beam.num_cols(); i_shape++)
      {
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        {
          pair_force_vector(i_shape * 3 + i_dim) = N_beam(i_shape) * diff(i_dim);
        }
      }
      for (unsigned int i_shape = 0; i_shape < N_edge.num_cols(); i_shape++)
      {
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        {
          pair_force_vector(N_beam.num_cols() * 3 + i_shape * 3 + i_dim) =
              -1.0 * N_edge(i_shape) * diff(i_dim);
        }
      }


      pair_force_vector.scale(force);
      pair_force_vector.scale(-1.0);

      // pos_beam.print(std::cout);
      // pos_edge.print(std::cout);


      const Core::Elements::Element* ele = this->face_element_->get_element();
      // std::cout << "\n is face1 " << ele->is_face_element();
      // std::cout << "\n node1 " << ele->nodes()[0]->id();
      // std::cout << "\n node1 " << ele->nodes()[0]->x()[0];
      // std::cout << "\n node1 " << ele->nodes()[0]->x()[1];
      // std::cout << "\n node1 " << ele->nodes()[0]->x()[2];
      // std::cout << "\n node2 " << ele->nodes()[1]->id();
      // std::cout << "\n node2 " << ele->nodes()[1]->x()[0];
      // std::cout << "\n node2 " << ele->nodes()[1]->x()[1];
      // std::cout << "\n node2 " << ele->nodes()[1]->x()[2];

      is_edge_to_edge_ = true;
      GeometryPair::evaluate_position(Core::FADUtils::cast_to_double(eta_a),
          GeometryPair::ElementDataToDouble<Beam>::to_double(this->ele1posref_),
          edge_to_edge_reference_);
      GeometryPair::evaluate_position(Core::FADUtils::cast_to_double(eta_a),
          GeometryPair::ElementDataToDouble<Beam>::to_double(this->ele1pos_),
          edge_to_edge_displacement_);
      edge_to_edge_displacement_ -= edge_to_edge_reference_;
      edge_to_edge_gap_ = Core::FADUtils::cast_to_double(gap);
      edge_to_edge_force_ = Core::FADUtils::cast_to_double(force_vec);
      edge_to_edge_normal_ = Core::FADUtils::cast_to_double(diff);
      edge_to_edge_normal_.scale(edge_to_edge_normal_.norm2());


      // GIDs of the pair and the force vector acting on the pair.
      const auto pair_gid = get_beam_to_surface_pair_gid_combined<Beam>(
          discret, *this->element1(), *this->face_element_);

      // If given, assemble force terms into the global vector.
      if (force_vector != nullptr)
      {
        std::vector<double> force_pair_double(pair_gid.size(), 0.0);
        for (unsigned int j_dof = 0; j_dof < pair_force_vector.num_rows(); j_dof++)
          force_pair_double[j_dof] = Core::FADUtils::cast_to_double(pair_force_vector(j_dof));
        force_vector->SumIntoGlobalValues(
            pair_gid.size(), pair_gid.data(), force_pair_double.data());
      }

      // If given, assemble force terms into the global stiffness matrix.
      if (stiffness_matrix != nullptr)
        for (unsigned int i_dof = 0; i_dof < pair_force_vector.num_rows(); i_dof++)
          for (unsigned int j_dof = 0; j_dof < Beam::n_dof_ + GeometryPair::t_line2::n_dof_;
              j_dof++)
            stiffness_matrix->fe_assemble(
                Core::FADUtils::cast_to_double(pair_force_vector(i_dof).dx(j_dof)), pair_gid[i_dof],
                pair_gid[j_dof]);
    }
    // Core::FADUtils::cast_to_double(this->ele1pos_.element_position_).print(std::cout);
    // Core::FADUtils::cast_to_double(this->edge_position_.element_position_).print(std::cout);
    // std::cout << "\neta_a " << Core::FADUtils::cast_to_double(eta_a);
    // std::cout << "\neta_b " << Core::FADUtils::cast_to_double(eta_b);
  }



  // At this point the pair is already evaluated in the current deformation state, so we don't have
  // to perform the projections or integration again, we can simply take the values previously
  // computed and multiply them with the Lagrange multipliers.

  // If there are no intersection segments, no contact terms will be assembled.
  const unsigned int n_segments = this->line_to_3D_segments_.size();
  if (n_segments == 0) return;

  if (n_segments > 1) FOUR_C_THROW("only one secment allowed");
  // Add end point penalty terms
  {
    // Get beam cross-section diameter
    const auto beam_ptr = dynamic_cast<const Discret::Elements::Beam3Base*>(this->element1());
    const double beam_cross_section_radius =
        beam_ptr->get_circular_cross_section_radius_for_interactions();

    // Check if the beam is at the end of a filament.
    std::array<bool, 2> is_end_node;
    for (unsigned int i_node = 0; i_node < 2; i_node++)
    {
      const Core::Nodes::Node* node = this->element1()->nodes()[i_node];
      is_end_node[i_node] = node->num_element() == 1;
    }

    for (const GeometryPair::LineSegment<ScalarType>& segment : this->line_to_3D_segments_)
    {
      GeometryPair::ProjectionPoint1DTo3D<ScalarType> point;
      if (abs(segment.get_eta_a() + 1.0) < 1e-7 and is_end_node[0])
      {
        std::cout << "\n end node a ids: " << this->element1()->id() << " "
                  << this->element2()->id();
        point = segment.get_start_point();
      }
      else if (abs(segment.get_eta_b() - 1.0) < 1e-7 and is_end_node[1])
      {
        std::cout << "\n end node b ids: " << this->element1()->id() << " "
                  << this->element2()->id();
        point = segment.get_end_point();
      }
      else
      {
        std::cout << "\n no end point fourn for: " << this->element1()->id() << " "
                  << this->element2()->id();
        continue;
      }

      // Evaluate the contact kinematics
      const auto [_1, _2, surface_normal, gap] =
          this->evaluate_contact_kinematics_at_projection_point(point, beam_cross_section_radius);

      // Get the contact force.
      ScalarType force = penalty_force(gap, *this->params()->beam_to_solid_surface_contact_params(),
          2 * beam_cross_section_radius);

      is_endpoint_ = true;
      GeometryPair::evaluate_position(Core::FADUtils::cast_to_double(point.get_eta()),
          GeometryPair::ElementDataToDouble<Beam>::to_double(this->ele1posref_),
          endpoint_reference_);
      GeometryPair::evaluate_position(Core::FADUtils::cast_to_double(point.get_eta()),
          GeometryPair::ElementDataToDouble<Beam>::to_double(this->ele1pos_),
          endpoint_displacement_);
      endpoint_displacement_ -= endpoint_reference_;
      endpoint_gap_ = Core::FADUtils::cast_to_double(gap);
      endpoint_normal_ = Core::FADUtils::cast_to_double(surface_normal);
      endpoint_normal_.scale(endpoint_normal_.norm2());
      endpoint_force_ = endpoint_normal_;
      endpoint_force_.scale(Core::FADUtils::cast_to_double(force));

      std::cout << "\nendoiint contact active, gap: " << Core::FADUtils::cast_to_double(gap);

      // Get the shape function matrices.
      Core::LinAlg::Matrix<Beam::n_dof_ + Surface::n_dof_, 1, ScalarType>
          gap_variation_times_normal;
      Core::LinAlg::Matrix<1, Beam::n_nodes_ * Beam::n_val_, ScalarType> N_beam;
      Core::LinAlg::Matrix<1, Surface::n_nodes_ * Surface::n_val_, ScalarType> N_surface;
      GeometryPair::EvaluateShapeFunction<Beam>::evaluate(
          N_beam, point.get_eta(), this->ele1pos_.shape_function_data_);
      GeometryPair::EvaluateShapeFunction<Surface>::evaluate(N_surface, point.get_xi(),
          this->face_element_->get_face_element_data().shape_function_data_);

      // Calculate the variation of the gap function multiplied with the surface normal vector.
      for (unsigned int i_shape = 0; i_shape < N_beam.num_cols(); i_shape++)
      {
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        {
          gap_variation_times_normal(i_shape * 3 + i_dim) = N_beam(i_shape) * surface_normal(i_dim);
        }
      }
      for (unsigned int i_shape = 0; i_shape < N_surface.num_cols(); i_shape++)
      {
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        {
          gap_variation_times_normal(N_beam.num_cols() * 3 + i_shape * 3 + i_dim) =
              -1.0 * N_surface(i_shape) * surface_normal(i_dim);
        }
      }

      // Add the Gauss point contributions to the pair force vector.
      Core::LinAlg::Matrix<Beam::n_dof_ + Surface::n_dof_, 1, ScalarType> pair_force_vector;
      gap_variation_times_normal.scale(force);
      pair_force_vector -= gap_variation_times_normal;

      // GIDs of the pair and the force vector acting on the pair.
      const auto pair_gid = get_beam_to_surface_pair_gid_combined<Beam>(
          discret, *this->element1(), *this->face_element_);

      // If given, assemble force terms into the global vector.
      if (force_vector != nullptr)
      {
        std::vector<double> force_pair_double(pair_gid.size(), 0.0);
        for (unsigned int j_dof = 0; j_dof < pair_force_vector.num_rows(); j_dof++)
          force_pair_double[j_dof] = Core::FADUtils::cast_to_double(pair_force_vector(j_dof));
        force_vector->SumIntoGlobalValues(
            pair_gid.size(), pair_gid.data(), force_pair_double.data());
      }

      // If given, assemble force terms into the global stiffness matrix.
      if (stiffness_matrix != nullptr)
        for (unsigned int i_dof = 0; i_dof < pair_force_vector.num_rows(); i_dof++)
          for (unsigned int j_dof = 0; j_dof < pair_gid.size(); j_dof++)
            stiffness_matrix->fe_assemble(
                Core::FADUtils::cast_to_double(pair_force_vector(i_dof).dx(j_dof)), pair_gid[i_dof],
                pair_gid[j_dof]);
    }
  }

  // This pair only gives contributions to the stiffness matrix
  if (stiffness_matrix == nullptr) return;

  // Get the Lagrange multipliers DOF vector for this pair
  const auto& [lambda_gid_pos, _] = mortar_manager->location_vector(*this);
  std::vector<double> lambda_pos_vector = Core::FE::extract_values(global_lambda, lambda_gid_pos);
  const auto lambda_pos = Core::LinAlg::Matrix<Mortar::n_dof_, 1, double>(lambda_pos_vector.data());

  // Multiply with the matrices evaluated in evaluate_and_assemble_mortar_contributions
  Core::LinAlg::Matrix<Beam::n_dof_, 1, ScalarType> force_beam(Core::LinAlg::Initialization::zero);
  force_beam.multiply(beam_shape_times_normal_times_lambda_shape_, lambda_pos);
  Core::LinAlg::Matrix<Surface::n_dof_, 1, ScalarType> force_surface(
      Core::LinAlg::Initialization::zero);
  force_surface.multiply(surface_shape_times_normal_times_lambda_shape_, lambda_pos);

  // Assemble the terms to the global stiffness matrix
  Core::LinAlg::Matrix<Beam::n_dof_, 1, int> beam_centerline_gid;
  BeamInteraction::Utils::get_element_centerline_gid_indices(
      discret, this->element1(), beam_centerline_gid);
  const std::vector<int>& patch_gid = this->face_element_->get_patch_gid();

  for (unsigned int i_beam = 0; i_beam < Beam::n_dof_; i_beam++)
  {
    for (unsigned int j_beam = 0; j_beam < Beam::n_dof_; j_beam++)
    {
      const double value = force_beam(i_beam).dx(j_beam);
      stiffness_matrix->fe_assemble(
          value, beam_centerline_gid(i_beam), beam_centerline_gid(j_beam));
    }
    for (unsigned int j_patch = 0; j_patch < patch_gid.size(); j_patch++)
    {
      const double value = force_beam(i_beam).dx(Beam::n_dof_ + j_patch);
      stiffness_matrix->fe_assemble(value, beam_centerline_gid(i_beam), patch_gid[j_patch]);
    }
  }
  for (unsigned int i_surface = 0; i_surface < Surface::n_dof_; i_surface++)
  {
    for (unsigned int j_beam = 0; j_beam < Beam::n_dof_; j_beam++)
    {
      const double value = force_surface(i_surface).dx(j_beam);
      stiffness_matrix->fe_assemble(value, patch_gid[i_surface], beam_centerline_gid(j_beam));
    }
    for (unsigned int j_patch = 0; j_patch < patch_gid.size(); j_patch++)
    {
      const double value = force_surface(i_surface).dx(Beam::n_dof_ + j_patch);
      stiffness_matrix->fe_assemble(value, patch_gid[i_surface], patch_gid[j_patch]);
    }
  }
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface, typename Mortar>
void BeamInteraction::BeamToSolidSurfaceContactPairMortar<ScalarType, Beam, Surface,
    Mortar>::get_pair_visualization(std::shared_ptr<BeamToSolidVisualizationOutputWriterBase>
                                        visualization_writer,
    Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base method.
  base_class::get_pair_visualization(visualization_writer, visualization_params);
  if (is_endpoint_)
  {
    std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization> visualization =
        visualization_writer->get_visualization_writer("btss-contact-end-point-forces");
    auto& visualization_data = visualization->get_visualization_data();
    std::vector<double>& point_coordinates = visualization_data.get_point_coordinates(3);
    std::vector<double>& displacement =
        visualization_data.get_point_data<double>("displacement", 3);
    std::vector<double>& force_beam = visualization_data.get_point_data<double>("force_beam", 3);
    std::vector<double>& normal = visualization_data.get_point_data<double>("normal", 3);
    std::vector<double>& gap = visualization_data.get_point_data<double>("gap", 1);

    gap.push_back(endpoint_gap_);
    for (unsigned int dim = 0; dim < 3; dim++)
    {
      point_coordinates.push_back(endpoint_reference_(dim));
      displacement.push_back(endpoint_displacement_(dim));
      force_beam.push_back(endpoint_force_(dim));
      normal.push_back(endpoint_normal_(dim));
    }
    // "btss-contact-edge-to-edge-forces"
  }
  if (is_edge_to_edge_)
  {
    std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization> visualization =
        visualization_writer->get_visualization_writer("btss-contact-edge-to-edge-forces");
    auto& visualization_data = visualization->get_visualization_data();
    std::vector<double>& point_coordinates = visualization_data.get_point_coordinates(3);
    std::vector<double>& displacement =
        visualization_data.get_point_data<double>("displacement", 3);
    std::vector<double>& force_beam = visualization_data.get_point_data<double>("force_beam", 3);
    std::vector<double>& normal = visualization_data.get_point_data<double>("normal", 3);
    std::vector<double>& gap = visualization_data.get_point_data<double>("gap", 1);

    gap.push_back(edge_to_edge_gap_);
    for (unsigned int dim = 0; dim < 3; dim++)
    {
      point_coordinates.push_back(edge_to_edge_reference_(dim));
      displacement.push_back(edge_to_edge_displacement_(dim));
      force_beam.push_back(edge_to_edge_force_(dim));
      normal.push_back(edge_to_edge_normal_(dim));
    }
    // "btss-contact-edge-to-edge-forces"
  }
  std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization> visualization_continuous =
      visualization_writer->get_visualization_writer("btss-contact-mortar-continuous");
  if (visualization_continuous == nullptr) return;

  const std::shared_ptr<const BeamToSolidSurfaceVisualizationOutputParams>& output_params_ptr =
      visualization_params.get<std::shared_ptr<const BeamToSolidSurfaceVisualizationOutputParams>>(
          "btss-output_params_ptr");
  const bool write_unique_ids = output_params_ptr->get_write_unique_ids_flag();

  if (visualization_continuous != nullptr)
  {
    // Get beam cross-section diameter
    const auto beam_ptr = dynamic_cast<const Discret::Elements::Beam3Base*>(this->element1());
    const double beam_cross_section_radius =
        beam_ptr->get_circular_cross_section_radius_for_interactions();

    // Setup variables.
    Core::LinAlg::Matrix<3, 1, ScalarType> r0_beam;
    Core::LinAlg::Matrix<3, 1, ScalarType> u;
    Core::LinAlg::Matrix<1, 1, ScalarType> lambda_scalar;

    // Get the mortar manager and the global lambda vector, those objects will be used to get the
    // discrete Lagrange multiplier values for this pair.
    std::shared_ptr<const BeamInteraction::BeamToSolidMortarManager> mortar_manager =
        visualization_params.get<std::shared_ptr<const BeamInteraction::BeamToSolidMortarManager>>(
            "mortar_manager");
    std::shared_ptr<Core::LinAlg::Vector<double>> lambda =
        visualization_params.get<std::shared_ptr<Core::LinAlg::Vector<double>>>("lambda");

    // Get the lambda GIDs of this pair.
    auto q_lambda = GeometryPair::InitializeElementData<Mortar, double>::initialize(nullptr);
    const auto& [lambda_row_pos, _] = mortar_manager->location_vector(*this);
    std::vector<double> lambda_pair = Core::FE::extract_values(*lambda, lambda_row_pos);
    for (unsigned int i_dof = 0; i_dof < Mortar::n_dof_; i_dof++)
      q_lambda.element_position_(i_dof) = lambda_pair[i_dof];

    // Add the continuous values for the Lagrange multipliers.
    if (visualization_continuous != nullptr and this->line_to_3D_segments_.size() > 0)
    {
      unsigned int n_points = 0;
      for (const auto& segment : this->line_to_3D_segments_)
      {
        n_points += 2 + segment.get_number_of_projection_points();
      }
      auto& visualization_data = visualization_continuous->get_visualization_data();
      std::vector<double>& point_coordinates =
          visualization_data.get_point_coordinates(n_points * 3);
      std::vector<double>& displacement =
          visualization_data.get_point_data<double>("displacement", n_points * 3);
      std::vector<double>& lambda_vis =
          visualization_data.get_point_data<double>("lambda", n_points);
      std::vector<double>& lambda_times_normal_vis =
          visualization_data.get_point_data<double>("lambda_times_normal", n_points * 3);
      std::vector<double>& surface_normal_vis =
          visualization_data.get_point_data<double>("surface_normal", n_points * 3);
      std::vector<double>& gap_vis = visualization_data.get_point_data<double>("gap", n_points);
      std::vector<uint8_t>& cell_types = visualization_data.get_cell_types();
      std::vector<int32_t>& cell_offsets = visualization_data.get_cell_offsets();

      std::vector<int>* pair_point_beam_id = nullptr;
      std::vector<int>* pair_point_solid_id = nullptr;
      std::vector<int>* pair_cell_beam_id = nullptr;
      std::vector<int>* pair_cell_solid_id = nullptr;
      if (write_unique_ids)
      {
        pair_point_beam_id = &(visualization_data.get_point_data<int>("uid_0_pair_beam_id"));
        pair_point_solid_id = &(visualization_data.get_point_data<int>("uid_1_pair_solid_id"));
        pair_cell_beam_id = &(visualization_data.get_cell_data<int>("uid_0_pair_beam_id"));
        pair_cell_solid_id = &(visualization_data.get_cell_data<int>("uid_1_pair_solid_id"));
      }

      // Check if the beam is at the end of a filament.
      std::array<bool, 2> is_end_node;
      for (unsigned int i_node = 0; i_node < 2; i_node++)
      {
        const Core::Nodes::Node* node = this->element1()->nodes()[i_node];
        is_end_node[i_node] = node->num_element() == 1;
      }
      if (is_end_node[0] and is_end_node[1]) FOUR_C_THROW("only one endpoint possible");


      for (const auto& segment : this->line_to_3D_segments_)
      {
        for (const auto& point : segment.get_all_segment_points())
        {
          const auto [r_beam, r_surface, surface_normal, gap] =
              this->evaluate_contact_kinematics_at_projection_point(
                  point, beam_cross_section_radius);

          // Get the beam displacement
          GeometryPair::evaluate_position<Beam>(point.get_eta(), this->ele1posref_, r0_beam);
          u = r_beam;
          u -= r0_beam;

          Core::LinAlg::Matrix<Mortar::spatial_dim_, Mortar::n_dof_, ScalarType> N_lambda;
          GeometryPair::evaluate_shape_function_matrix<Mortar>(N_lambda, point.get_eta());

          // Get the Lagrange multiplier value
          // GeometryPair::evaluate_position<Mortar>(point.get_eta(), q_lambda, lambda_scalar);

          // r.clear();
          // for (unsigned int node = 0; node < ElementType::n_nodes_; node++)
          //   for (unsigned int dim = 0; dim < ElementType::spatial_dim_; dim++)
          //     for (unsigned int val = 0; val < ElementType::n_val_; val++)
          //       r(dim) += element_data.element_position_(
          //                     ElementType::spatial_dim_ * ElementType::n_val_ * node +
          //                     ElementType::spatial_dim_ * val + dim) *
          //                 N(ElementType::n_val_ * node + val);

          if (is_end_node[0])
          {
            ScalarType tmp = N_lambda(0);
            N_lambda(0) = 0.0;
            N_lambda(1) += tmp;

            // tmp = N_lambda_trial(0);
            // N_lambda_trial(0) = 0.0;
            // N_lambda_trial(1) += tmp;
          }
          else if (is_end_node[1])
          {
            ScalarType tmp = N_lambda(1);
            N_lambda(1) = 0.0;
            N_lambda(0) += tmp;

            // tmp = N_lambda_trial(1);
            // N_lambda_trial(1) = 0.0;
            // N_lambda_trial(0) += tmp;
          }

          lambda_scalar.multiply_nn(N_lambda, q_lambda.element_position_);

          // Add to output data.
          lambda_vis.push_back(Core::FADUtils::cast_to_double(lambda_scalar(0)));
          gap_vis.push_back(Core::FADUtils::cast_to_double(gap));
          for (unsigned int dim = 0; dim < 3; dim++)
          {
            point_coordinates.push_back(Core::FADUtils::cast_to_double(r0_beam(dim)));
            displacement.push_back(Core::FADUtils::cast_to_double(u(dim)));
            lambda_times_normal_vis.push_back(
                Core::FADUtils::cast_to_double(surface_normal(dim) * lambda_scalar(0)));
            surface_normal_vis.push_back(Core::FADUtils::cast_to_double(surface_normal(dim)));
          }
          if (write_unique_ids)
          {
            pair_point_beam_id->push_back(this->element1()->id());
            pair_point_solid_id->push_back(this->element2()->id());
          }
        }

        // Add the cell for this segment (poly line).
        cell_types.push_back(4);
        cell_offsets.push_back(point_coordinates.size() / 3);

        if (write_unique_ids)
        {
          pair_cell_beam_id->push_back(this->element1()->id());
          pair_cell_solid_id->push_back(this->element2()->id());
        }
      }
    }
  }
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface, typename Mortar>
ScalarType BeamInteraction::BeamToSolidSurfaceContactPairMortar<ScalarType, Beam, Surface,
    Mortar>::get_jacobian_for_configuration(const ScalarType& eta,
    const Inpar::BeamToSolid::BeamToSolidSurfaceContactMortarDefinedIn mortar_configuration) const
{
  Core::LinAlg::Matrix<3, 1, ScalarType> dr_beam;
  switch (mortar_configuration)
  {
    case Inpar::BeamToSolid::BeamToSolidSurfaceContactMortarDefinedIn::reference_configuration:
    {
      GeometryPair::evaluate_position_derivative1<Beam>(eta, this->ele1posref_, dr_beam);
      return Core::FADUtils::vector_norm(dr_beam);
    }
    case Inpar::BeamToSolid::BeamToSolidSurfaceContactMortarDefinedIn::current_configuration:
    {
      GeometryPair::evaluate_position_derivative1<Beam>(eta, this->ele1pos_, dr_beam);
      return Core::FADUtils::vector_norm(dr_beam);
    }
    default:
      FOUR_C_THROW("Got unexpected mortar configuration");
  }
}

/**
 * @brief Factory function templated on the type of beam element and the surface shape
 */
template <typename Beam, typename Surface, typename ScalarType>
std::shared_ptr<BeamInteraction::BeamContactPair>
beam_to_solid_surface_contact_pair_mortar_factory_template_beam_surface(
    const BeamInteraction::BeamToSolidSurfaceContactParams& beam_to_surface_contact_params)
{
  using namespace GeometryPair;

  switch (beam_to_surface_contact_params.get_mortar_shape_function_type())
  {
    case Inpar::BeamToSolid::BeamToSolidMortarShapefunctions::line2:
      return std::make_shared<BeamInteraction::BeamToSolidSurfaceContactPairMortar<ScalarType, Beam,
          Surface, t_line2_scalar>>();
    default:
      FOUR_C_THROW("Got unexpected mortar shape function");
  }
}

/**
 * @brief Factory function templated on the type of beam element
 */
template <typename Beam>
std::shared_ptr<BeamInteraction::BeamContactPair>
beam_to_solid_surface_contact_pair_mortar_factory_template_beam(
    const BeamInteraction::BeamToSolidSurfaceContactParams& beam_to_surface_contact_params,
    const Core::FE::CellType& surface_type)
{
  using namespace GeometryPair;

  switch (surface_type)
  {
    case Core::FE::CellType::quad4:
      return beam_to_solid_surface_contact_pair_mortar_factory_template_beam_surface<Beam, t_quad4,
          line_to_surface_patch_scalar_type_1st_order>(beam_to_surface_contact_params);
    case Core::FE::CellType::nurbs9:
      return beam_to_solid_surface_contact_pair_mortar_factory_template_beam_surface<Beam, t_nurbs9,
          line_to_surface_patch_scalar_type_fixed_size_1st_order<Beam, t_nurbs9>>(
          beam_to_surface_contact_params);
    default:
      FOUR_C_THROW("Got unexpected surface shape");
  }
}

/**
 *
 */
std::shared_ptr<BeamInteraction::BeamContactPair>
BeamInteraction::beam_to_solid_surface_contact_pair_mortar_factory(
    const BeamToSolidSurfaceContactParams& beam_to_surface_contact_params,
    const Core::FE::CellType& surface_type, const bool beam_is_hermite)
{
  using namespace GeometryPair;

  if (beam_is_hermite)
  {
    return beam_to_solid_surface_contact_pair_mortar_factory_template_beam<t_hermite>(
        beam_to_surface_contact_params, surface_type);
  }
  else
  {
    FOUR_C_THROW("Beam-to-solid contact with mortar is not implemented for linear beam elements");
  }
}

FOUR_C_NAMESPACE_CLOSE
