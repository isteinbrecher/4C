// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_beam_to_beam_point_coupling_pair.hpp"

#include "4C_beam3_reissner.hpp"
#include "4C_beam3_triad_interpolation_local_rotation_vectors.hpp"
#include "4C_beaminteraction_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_geometry_pair_access_traits.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_geometry_pair_constants.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_line_to_line.hpp"
#include "4C_linalg_fevector.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_utils_exceptions.hpp"

#include <array>
#include <tuple>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2, n_dof_beam_2>::setup()
{
  // This pair only works for Simo Reissner beam elements.
  const auto check_simo_reissner_beam = [](auto element)
  {
    const bool is_sr_beam = dynamic_cast<const Discret::Elements::Beam3r*>(element) != nullptr;
    if (!is_sr_beam)
      FOUR_C_THROW("The BeamToBeamPointCouplingPair only works for Simo Reissner beams");
  };
  check_simo_reissner_beam(this->element1());
  check_simo_reissner_beam(this->element2());

  this->issetup_ = true;
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::evaluate_and_assemble(const std::shared_ptr<const Core::FE::Discretization>&
                                             discret,
    const std::shared_ptr<Core::LinAlg::FEVector<double>>& force_vector,
    const std::shared_ptr<Core::LinAlg::SparseMatrix>& stiffness_matrix,
    const std::shared_ptr<const Core::LinAlg::Vector<double>>& displacement_vector)
{
  check_init_setup();

  if (use_closest_point_projection_ && !closest_point_projection_evaluated_)
  {
    evaluate_closest_point_projection();
    closest_point_projection_evaluated_ = true;
  }
  if (!evaluate_pair_) return;

  constexpr std::array<unsigned int, 2> n_dof_beam = {n_dof_beam_1, n_dof_beam_2};

  const std::array<const Core::Elements::Element*, 2> beam_elements = {
      this->element1(), this->element2()};

  // Transformation matrices for the coupling constraints
  Core::LinAlg::Matrix<n_dof_total, 12> left_transformation_matrix(
      Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<12, n_dof_total> right_transformation_matrix(
      Core::LinAlg::Initialization::zero);

  // Variables holding the GIDs of the beam element DOFs.
  std::array<int, n_dof_total> pair_gid{-1};

  // Kinematic variables.
  std::array<Core::LinAlg::Matrix<3, 1>, 2> r_ref;
  std::array<Core::LinAlg::Matrix<3, 1>, 2> r;
  std::array<std::vector<double>, 2> element_displacement;
  std::array<Core::LinAlg::Matrix<4, 1, scalar_type_rot>, 2> cross_section_quaternion;
  std::array<Core::LinAlg::Matrix<4, 1, double>, 2> cross_section_quaternion_ref;

  // Evaluate all quantities at the current state for both beams.
  {
    for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
    {
      const auto* beam_element =
          dynamic_cast<const Discret::Elements::Beam3Base*>(beam_elements[i_beam]);

      std::vector<int> lm, lmowner, lmstride;
      beam_element->location_vector(*discret, lm, lmowner, lmstride);
      auto eledisp = Core::FE::extract_values(*displacement_vector, lm);
      auto eledisp_ref = std::vector<double>(eledisp.size(), 0.0);
      element_displacement[i_beam] = eledisp;

      for (unsigned int i_dof = 0; i_dof < n_dof_beam[i_beam]; i_dof++)
        pair_gid[i_dof + i_beam * n_dof_beam_1] = lm[i_dof];

      beam_element->get_pos_at_xi(r[i_beam], position_in_parameterspace_[i_beam], eledisp);
      beam_element->get_pos_at_xi(r_ref[i_beam], position_in_parameterspace_[i_beam], eledisp_ref);

      Core::LinAlg::Matrix<3, 3> triad;
      Core::LinAlg::Matrix<3, 1> psi_double;
      Core::LinAlg::Matrix<3, 1, scalar_type_rot> psi;
      beam_element->get_triad_at_xi(triad, position_in_parameterspace_[i_beam], eledisp);
      Core::LinAlg::Matrix<4, 1, double> quaternion_double;
      Core::LargeRotations::triadtoquaternion(triad, quaternion_double);
      Core::LargeRotations::quaterniontoangle(quaternion_double, psi_double);
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        psi(i_dim) = Core::FADUtils::HigherOrderFadValue<scalar_type_rot>::apply(
            6, i_beam * 3 + i_dim, psi_double(i_dim));
      Core::LargeRotations::angletoquaternion(psi, cross_section_quaternion[i_beam]);

      Core::LinAlg::Matrix<3, 3> triad_ref;
      beam_element->get_triad_at_xi(triad_ref, position_in_parameterspace_[i_beam], eledisp_ref);
      Core::LargeRotations::triadtoquaternion(triad_ref, cross_section_quaternion_ref[i_beam]);

      Core::LinAlg::SerialDenseMatrix trafomatrix_left;
      trafomatrix_left.shape(6, n_dof_beam[i_beam]);
      beam_element->get_generalized_interpolation_matrix_variations_at_xi(
          trafomatrix_left, position_in_parameterspace_[i_beam], eledisp);

      Core::LinAlg::SerialDenseMatrix trafomatrix_right;
      trafomatrix_right.shape(6, n_dof_beam[i_beam]);
      beam_element->get_generalized_interpolation_matrix_increments_at_xi(
          trafomatrix_right, position_in_parameterspace_[i_beam], eledisp);

      for (unsigned int i_dof = 0; i_dof < n_dof_beam[i_beam]; i_dof++)
      {
        for (unsigned int i_dir = 0; i_dir < 6; i_dir++)
        {
          const unsigned int i_row = i_dir + i_beam * 6;
          const unsigned int i_col = i_dof + i_beam * n_dof_beam_1;
          left_transformation_matrix(i_col, i_row) = trafomatrix_left(i_dir, i_dof);
          right_transformation_matrix(i_row, i_col) = trafomatrix_right(i_dir, i_dof);
        }
      }
    }
  }

  // Positional coupling terms
  const auto [constraint_position, constraint_position_lin_kinematic, residuum_position_lin_lambda,
      evaluation_data_position] = evaluate_positional_coupling(r_ref, r,
      cross_section_quaternion_ref, cross_section_quaternion);

  // Rotational coupling terms
  const auto [constraint_rotation, constraint_rotation_lin_kinematic, residuum_rotation_lin_lambda,
      evaluation_data_rotation] =
      evaluate_rotational_coupling(cross_section_quaternion_ref, cross_section_quaternion);

  // Coupling residuum and stiffness
  Core::LinAlg::Matrix<12, 1> residuum(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<12, 12> stiffness(Core::LinAlg::Initialization::zero);

  // Penalty regularization positions
  Core::LinAlg::Matrix<3, 1> lambda_position = constraint_position;
  lambda_position.scale(penalty_parameter_pos_);
  Core::LinAlg::Matrix<12, 1> residuum_position;
  residuum_position.multiply_nn(residuum_position_lin_lambda, lambda_position);
  residuum += residuum_position;

  Core::LinAlg::Matrix<12, 12> stiffness_position;
  stiffness_position.multiply_nn(residuum_position_lin_lambda, constraint_position_lin_kinematic);
  stiffness_position.scale(penalty_parameter_pos_);
  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    Core::LinAlg::Matrix<3, 3, double> skew_lambda;
    Core::LargeRotations::computespin(skew_lambda, lambda_position);
    Core::LinAlg::Matrix<3, 3> temp_matrix;
    temp_matrix.multiply_nn(skew_lambda, evaluation_data_position[i_beam]);
    for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
    {
      for (unsigned int j_dir = 0; j_dir < 3; j_dir++)
      {
        stiffness_position(i_dir + 3 + 6 * i_beam, j_dir + 3 + 6 * i_beam) -=
            temp_matrix(i_dir, j_dir);
      }
    }
  }
  stiffness += stiffness_position;

  // Penalty regularization rotations
  Core::LinAlg::Matrix<3, 1> lambda_rotation = constraint_rotation;
  lambda_rotation.scale(penalty_parameter_rot_);
  Core::LinAlg::Matrix<12, 1> residuum_rotation;
  residuum_rotation.multiply_nn(residuum_rotation_lin_lambda, lambda_rotation);
  residuum += residuum_rotation;

  Core::LinAlg::Matrix<12, 12> stiffness_rot;
  stiffness_rot.multiply_nn(residuum_rotation_lin_lambda, constraint_rotation_lin_kinematic);
  stiffness_rot.scale(penalty_parameter_rot_);
  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3; j++)
      {
        for (unsigned int l = 0; l < 3; l++)
        {
          stiffness_rot(i + 3, l + 3 + 6 * i_beam) -=
              evaluation_data_rotation[i_beam][i][j][l] * lambda_rotation(j);
          stiffness_rot(i + 9, l + 3 + 6 * i_beam) +=
              evaluation_data_rotation[i_beam][i][j][l] * lambda_rotation(j);
        }
      }
    }
  }
  stiffness += stiffness_rot;

  // Map residuum and stiffness to element DOFs
  Core::LinAlg::Matrix<n_dof_total, 1> residuum_pair{Core::LinAlg::Initialization::zero};
  Core::LinAlg::Matrix<n_dof_total, n_dof_total> stiffness_pair{Core::LinAlg::Initialization::zero};
  residuum_pair.multiply(left_transformation_matrix, residuum);
  Core::LinAlg::Matrix<n_dof_total, 12> temp_matrix;
  temp_matrix.multiply(left_transformation_matrix, stiffness);
  stiffness_pair.multiply(temp_matrix, right_transformation_matrix);

  // Add stiffness contrbutions due to the beam forumlation
  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    const auto* beam_element =
        dynamic_cast<const Discret::Elements::Beam3Base*>(beam_elements[i_beam]);

    Core::LinAlg::SerialDenseVector force(6, true);
    for (unsigned int i_dim = 0; i_dim < 6; i_dim++) force(i_dim) = residuum(i_dim + 6 * i_beam);

    Core::LinAlg::SerialDenseMatrix stiffness_beam;
    stiffness_beam.shape(n_dof_beam[i_beam], n_dof_beam[i_beam]);
    beam_element->get_stiffmat_resulting_from_generalized_interpolation_matrix_at_xi(
        stiffness_beam, position_in_parameterspace_[i_beam], element_displacement[i_beam], force);

    for (unsigned int i_dof = 0; i_dof < n_dof_beam[i_beam]; i_dof++)
    {
      for (unsigned int j_dof = 0; j_dof < n_dof_beam[i_beam]; j_dof++)
      {
        stiffness_pair(i_dof + i_beam * n_dof_beam_1, j_dof + i_beam * n_dof_beam_1) +=
            stiffness_beam(i_dof, j_dof);
      }
    }
  }

  // Add the coupling terms into the global vector and matrix.
  if (force_vector != nullptr)
    force_vector->sum_into_global_values(pair_gid.size(), pair_gid.data(), residuum_pair.data());
  if (stiffness_matrix != nullptr)
  {
    for (unsigned int i_dof = 0; i_dof < pair_gid.size(); i_dof++)
    {
      for (unsigned int j_dof = 0; j_dof < pair_gid.size(); j_dof++)
      {
        stiffness_matrix->fe_assemble(
            stiffness_pair(i_dof, j_dof), pair_gid[i_dof], pair_gid[j_dof]);
      }
    }
  }
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::evaluate_closest_point_projection()
{
  const std::array<const Core::Elements::Element*, 2> beam_elements = {
      this->element1(), this->element2()};

  auto get_beam_ref_data = [&]<typename Element, unsigned int i_beam>()
  {
    auto beam_pos_ref =
        GeometryPair::InitializeElementData<Element, double>::initialize(beam_elements[i_beam]);

    const auto* beam_element =
        dynamic_cast<const Discret::Elements::Beam3Base*>(beam_elements[i_beam]);
    std::vector<double> zero_reference_displacement(
        std::get<i_beam>(std::make_tuple(n_dof_beam_1, n_dof_beam_2)), 0.0);
    std::vector<double> reference_dof_centerline(Element::n_dof_, 0.0);
    beam_element->extract_centerline_dof_values_from_element_state_vector(
        zero_reference_displacement, reference_dof_centerline, true);

    // Get the current and reference position.
    for (unsigned int i_dof = 0; i_dof < Element::n_dof_; i_dof++)
    {
      beam_pos_ref.element_position_(i_dof) = reference_dof_centerline[i_dof];
    }

    return beam_pos_ref;
  };

  // Closest point projection between the two curves
  const auto beam_pos_ref_1 = get_beam_ref_data.template operator()<Beam1, 0>();
  const auto beam_pos_ref_2 = get_beam_ref_data.template operator()<Beam2, 1>();
  const auto projection_result = GeometryPair::line_to_line_closest_point_projection(beam_pos_ref_1,
      beam_pos_ref_2, position_in_parameterspace_[0], position_in_parameterspace_[1]);

  if (projection_result != GeometryPair::ProjectionResult::projection_found_valid)
  {
    // No projection was found
    evaluate_pair_ = false;
    return;
  }

  // Check the projection distance
  Core::LinAlg::Matrix<3, 1> diff{Core::LinAlg::Initialization::zero};
  Core::LinAlg::Matrix<3, 1> r;
  GeometryPair::evaluate_position<Beam1>(position_in_parameterspace_[0], beam_pos_ref_1, r);
  diff -= r;
  GeometryPair::evaluate_position<Beam2>(position_in_parameterspace_[1], beam_pos_ref_2, r);
  diff += r;

  double beam_radii = 0.0;
  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    const auto* beam_ptr = dynamic_cast<const Discret::Elements::Beam3Base*>(beam_elements[i_beam]);
    beam_radii += beam_ptr->get_circular_cross_section_radius_for_interactions();
  }
  if (projection_valid_factor_ * beam_radii < diff.norm2())
  {
    evaluate_pair_ = false;
    return;
  }

  // Make sure that we have unique pairs for projections directly on nodes
  if (not line_to_line_evaluation_data_->evaluate_projection_coordinates(
          beam_elements, position_in_parameterspace_))
  {
    evaluate_pair_ = false;
    return;
  }

  // All checks passed, this pair shall be evaluated.
  evaluate_pair_ = true;
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
std::tuple<Core::LinAlg::Matrix<3, 1>, Core::LinAlg::Matrix<3, 12>, Core::LinAlg::Matrix<12, 3>,
    std::array<Core::LinAlg::Matrix<3, 3, double>, 2>>
BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::evaluate_positional_coupling(const std::array<Core::LinAlg::Matrix<3, 1>, 2>&
                                                    r_ref,
    const std::array<Core::LinAlg::Matrix<3, 1>, 2>& r,
    const std::array<Core::LinAlg::Matrix<4, 1, double>, 2>& cross_section_quaternion_ref,
    const std::array<Core::LinAlg::Matrix<4, 1, scalar_type_rot>, 2>& cross_section_quaternion)
{
  // Coupling vectors and matrices
  Core::LinAlg::Matrix<3, 1> constraint_position(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<3, 12> constraint_position_lin_kinematic(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<12, 3> residuum_position_lin_lambda(Core::LinAlg::Initialization::zero);
  std::array<Core::LinAlg::Matrix<3, 3, double>, 2> evaluation_data_position;

  // Note: If the two points are not at the same position in the reference configuration, then even
  // for the positional coupling terms there will be contributions to the rotational equilibrium,
  // i.e., moments. To keep things simple, we first evaluate the terms that will affect the
  // positional DOFs and then we compute the terms that will affect the rotational DOFs.

  // Constraint (pure positional terms)
  constraint_position = r[1];
  constraint_position -= r[0];

  // Linearizations (pure positional terms)
  for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
  {
    constraint_position_lin_kinematic(i_dir, i_dir) = -1.0;
    constraint_position_lin_kinematic(i_dir, 6 + i_dir) = 1.0;

    residuum_position_lin_lambda(i_dir, i_dir) = -1.0;
    residuum_position_lin_lambda(6 + i_dir, i_dir) = 1.0;
  }

  // Terms arising due to non matching points in the reference configurtaion.
  {
    // Evaluate reference offset in the material and spatial frames.
    Core::LinAlg::Matrix<3, 1, double> director_21_ref_spatial = r_ref[1];
    director_21_ref_spatial -= r_ref[0];
    std::array<Core::LinAlg::Matrix<3, 1, double>, 2> director_21_ref_material;
    for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
    {
      Core::LinAlg::Matrix<3, 3, double> rotation_matrix;
      Core::LargeRotations::quaterniontotriad(
          cross_section_quaternion_ref[i_beam], rotation_matrix);
      director_21_ref_material[i_beam].multiply_tn(rotation_matrix, director_21_ref_spatial);
    }

    // Add to constraint vector.
    // We add a half of the reference offset vector in each frame of the two beams. By doing it this
    // way, we obtain a coupling formulation invariant to the ordering of the cross-sections in this
    // pair.
    std::array<Core::LinAlg::Matrix<3, 1, double>, 2> director_21_spatial_half;
    for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
    {
      Core::LinAlg::Matrix<3, 3, double> rotation_matrix;
      Core::LargeRotations::quaterniontotriad(
          Core::FADUtils::cast_to_double(cross_section_quaternion[i_beam]), rotation_matrix);
      director_21_spatial_half[i_beam].multiply_nn(
          rotation_matrix, director_21_ref_material[i_beam]);
      director_21_spatial_half[i_beam].scale(0.5);
      constraint_position -= director_21_spatial_half[i_beam];
    }

    // Add the linearizations of the previous constraint terms.
    for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
    {
      Core::LinAlg::Matrix<3, 3, double> skew_director_21_spatial_half;
      Core::LargeRotations::computespin(
          skew_director_21_spatial_half, director_21_spatial_half[i_beam]);
      evaluation_data_position[i_beam] = skew_director_21_spatial_half;

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
      {
        for (unsigned int j_dir = 0; j_dir < 3; j_dir++)
        {
          constraint_position_lin_kinematic(i_dir, j_dir + 3 + 6 * i_beam) +=
              skew_director_21_spatial_half(i_dir, j_dir);

          residuum_position_lin_lambda(j_dir + 3 + 6 * i_beam, i_dir) -=
              skew_director_21_spatial_half(j_dir, i_dir);
        }
      }
    }
  }

  return {constraint_position, constraint_position_lin_kinematic, residuum_position_lin_lambda,
      evaluation_data_position};
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
std::tuple<Core::LinAlg::Matrix<3, 1>, Core::LinAlg::Matrix<3, 12>, Core::LinAlg::Matrix<12, 3>,
    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 2>>
BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::evaluate_rotational_coupling(const std::array<Core::LinAlg::Matrix<4, 1, double>,
                                                    2>& cross_section_quaternion_ref,
    const std::array<Core::LinAlg::Matrix<4, 1, scalar_type_rot>, 2>& cross_section_quaternion)
{
  // Coupling vectors and matrices
  Core::LinAlg::Matrix<3, 1> constraint_rotation(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<3, 12> constraint_rotation_lin_kinematic(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<12, 3> residuum_rotation_lin_lambda(Core::LinAlg::Initialization::zero);
  std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 2> evaluation_data_rotation{};

  // Get the relative rotation vector between the two cross sections.
  Core::LinAlg::Matrix<4, 1, scalar_type_rot> temp_quaternion_1, temp_quaternion_2, quaternion_rel;
  Core::LinAlg::Matrix<3, 1, scalar_type_rot> psi_rel;
  Core::LinAlg::Matrix<4, 1, scalar_type_rot> quaternion_0_inv =
      Core::LargeRotations::inversequaternion(cross_section_quaternion[0]);
  Core::LinAlg::Matrix<4, 1, double> quaternion_1_ref_inv =
      Core::LargeRotations::inversequaternion(cross_section_quaternion_ref[1]);
  Core::LargeRotations::quaternionproduct(
      quaternion_0_inv, cross_section_quaternion_ref[0], temp_quaternion_1);
  Core::LargeRotations::quaternionproduct(
      temp_quaternion_1, quaternion_1_ref_inv, temp_quaternion_2);
  Core::LargeRotations::quaternionproduct(
      temp_quaternion_2, cross_section_quaternion[1], quaternion_rel);
  Core::LargeRotations::quaterniontoangle(quaternion_rel, psi_rel);

  // Transformation matrix
  const Core::LinAlg::Matrix<3, 3, scalar_type_rot> T_psi_rel =
      Core::LargeRotations::tmatrix(psi_rel);

  constraint_rotation = Core::FADUtils::cast_to_double(psi_rel);
  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    const double load_factor = i_beam == 0 ? -1.0 : 1.0;

    // Get the rotation angle
    Core::LinAlg::Matrix<3, 1, double> psi_beam;
    Core::LargeRotations::quaterniontoangle(
        Core::FADUtils::cast_to_double(cross_section_quaternion[i_beam]), psi_beam);

    const Core::LinAlg::Matrix<3, 3, double> T_psi_beam = Core::LargeRotations::tmatrix(psi_beam);
    Core::LinAlg::Matrix<3, 3, double> d_psi_rel_d_psi_beam;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    {
      for (unsigned int j_dim = 0; j_dim < 3; j_dim++)
      {
        d_psi_rel_d_psi_beam(i_dim, j_dim) = psi_rel(i_dim).dx(j_dim + i_beam * 3);
      }
    }

    Core::LinAlg::Matrix<3, 3, double> constraint_lin_psi_beam;
    constraint_lin_psi_beam.multiply_nn(d_psi_rel_d_psi_beam, T_psi_beam);
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    {
      for (unsigned int j_dim = 0; j_dim < 3; j_dim++)
      {
        constraint_rotation_lin_kinematic(i_dim, j_dim + 3 + i_beam * 6) =
            constraint_lin_psi_beam(i_dim, j_dim);
        // Note: Here we insert the transposed, thus the switched ordering in the target but not
        // the source.
        residuum_rotation_lin_lambda(j_dim + 3 + i_beam * 6, i_dim) =
            load_factor * Core::FADUtils::cast_to_double(T_psi_rel(i_dim, j_dim));
      }
    }
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3; j++)
      {
        for (unsigned int l = 0; l < 3; l++)
        {
          for (unsigned int k = 0; k < 3; k++)
          {
            evaluation_data_rotation[i_beam][i][j][l] +=
                T_psi_rel(j, i).dx(k + 3 * i_beam) * T_psi_beam(k, l);
          }
        }
      }
    }
  }

  return {constraint_rotation, constraint_rotation_lin_kinematic, residuum_rotation_lin_lambda,
      evaluation_data_rotation};
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2, n_dof_beam_2>::print(
    std::ostream& out) const
{
  check_init_setup();

  // Print some general information: Element IDs and dofvecs.
  out << "\n------------------------------------------------------------------------";
  out << "\nInstance of BeamToBeamPenaltyPointCouplingPair"
      << "\nBeam1 EleGID:  " << element1()->id() << "\nBeam2 EleGID: " << element2()->id();
  out << "------------------------------------------------------------------------\n";
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::print_summary_one_line_per_active_segment_pair(std::ostream& out) const
{
  check_init_setup();

  out << "Beam-to-beam point coupling pair, beam1 gid: " << element1()->id()
      << " beam2 gid: " << element2()->id() << ", position in parameter space: ["
      << position_in_parameterspace_[0] << ", " << position_in_parameterspace_[1] << "]\n";
}

/**
 * Explicit template initialization of template class.
 */
namespace BeamInteraction
{
  using namespace GeometryPair;

  template class BeamToBeamPointCouplingPair<t_hermite, 21, t_hermite, 21>;
}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE
