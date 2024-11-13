// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_submodel_evaluator_potential.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_utils.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_crosslinker_handler.hpp"
#include "4C_beaminteraction_potential_pair.hpp"
#include "4C_beaminteraction_potential_params.hpp"
#include "4C_beaminteraction_potential_runtime_visualization_output_params.hpp"
#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"
#include "4C_comm_mpi_utils.hpp"
#include "4C_comm_utils_gid_vector.hpp"
#include "4C_fem_geometric_search_bounding_volume.hpp"
#include "4C_fem_geometric_search_bvh.hpp"
#include "4C_fem_geometric_search_distributed_tree.hpp"
#include "4C_fem_geometric_search_params.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_io_visualization_manager.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"
#include "4C_utils_exceptions.hpp"

#include <NOX_Solver_Generic.H>
#include <Teuchos_TimeMonitor.hpp>

#include <unordered_set>
FOUR_C_NAMESPACE_OPEN

namespace
{
  namespace LengthToEdgeImplementation

  {
    struct DataMaps
    {
      // data maps to determine prior element length for potential reduction strategy for single
      // length specific potential determination (maps are utilized for simple conversion into
      // vectors during communication)
      std::unordered_map<int, double> ele_gid_length_map;
      std::unordered_map<int, int> ele_gid_left_node_gid_map;
      std::unordered_map<int, int> ele_gid_right_node_gid_map;
      std::unordered_multimap<int, int> left_node_gid_ele_gid_map;
      std::unordered_multimap<int, int> right_node_gid_ele_gid_map;
    };

    // recursively determine length from beam element to the fiber's
    // end points for usage within potential reduction strategy
    double determine_length_to_edge(
        const DataMaps& maps, int ele_gid, int connecting_node_gid, double prior_length = 0.0)
    {
      // determine number of elements at connecting node
      // due to ghosting multiple entries of nodes->elements are possible
      std::unordered_set<int> elements_at_connecting_node;

      for (auto it = maps.left_node_gid_ele_gid_map.equal_range(connecting_node_gid).first;
           it != maps.left_node_gid_ele_gid_map.equal_range(connecting_node_gid).second; ++it)
      {
        elements_at_connecting_node.insert(it->second);
      }
      for (auto it = maps.right_node_gid_ele_gid_map.equal_range(connecting_node_gid).first;
           it != maps.right_node_gid_ele_gid_map.equal_range(connecting_node_gid).second; ++it)
      {
        elements_at_connecting_node.insert(it->second);
      }

      // if only one element is present => edge of fiber is reached
      // start recursive length evaluation if two elements are present at connecting node
      if (elements_at_connecting_node.size() == 2)
      {
        // determine neighbor element
        int neighbor_ele_gid;
        for (const int& possible_neighbor_ele_gid : elements_at_connecting_node)
        {
          if (possible_neighbor_ele_gid != ele_gid) neighbor_ele_gid = possible_neighbor_ele_gid;
        }

        // determine next connecting node of neighbor element
        int neighbor_connecting_node_gid;
        if (maps.ele_gid_left_node_gid_map.at(neighbor_ele_gid) != connecting_node_gid)
          neighbor_connecting_node_gid = maps.ele_gid_left_node_gid_map.at(neighbor_ele_gid);
        else if (maps.ele_gid_right_node_gid_map.at(neighbor_ele_gid) != connecting_node_gid)
          neighbor_connecting_node_gid = maps.ele_gid_right_node_gid_map.at(neighbor_ele_gid);
        else
          FOUR_C_THROW("Next connecting node for prior length determination not found!");

        // add neighbor element length to prior length
        prior_length += maps.ele_gid_length_map.at(neighbor_ele_gid);

        // call function recursively for next neighbor
        prior_length = determine_length_to_edge(
            maps, neighbor_ele_gid, neighbor_connecting_node_gid, prior_length);
      }
      else if (elements_at_connecting_node.size() > 2)
      {
        FOUR_C_THROW(
            "More than two beam elements are connected via a single node! Determination of length "
            "to edge for potential reduction strategy is only possible for a maximum number of two "
            "elements per node!");
      }

      return prior_length;
    }
  }  // namespace LengthToEdgeImplementation
}  // namespace

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::BeamPotential()
{
  // clear stl stuff
  nearby_elements_map_.clear();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::setup()
{
  check_init();

  // init and setup beam to beam contact data container
  beam_potential_params_ptr_ = std::make_shared<BEAMINTERACTION::BeamPotentialParams>();
  beam_potential_params().init(g_state().get_time_n());
  beam_potential_params().setup();

  print_console_welcome_message(std::cout);

  // build runtime visualization writer if desired
  if (beam_potential_params().runtime_output()) init_output_runtime_beam_potential();

  // set flag
  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::post_setup()
{
  check_init_setup();

  if (beam_potential_params().potential_reduction_length() != -1.0)
    setup_potential_reduction_strategy();

  nearby_elements_map_.clear();
  find_and_store_neighboring_elements();
  create_beam_potential_element_pairs();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::init_submodel_dependencies(
    std::shared_ptr<Solid::ModelEvaluator::BeamInteraction::Map> const submodelmap)
{
  check_init_setup();
  // no active influence on other submodels
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::reset()
{
  check_init_setup();

  std::vector<std::shared_ptr<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
  {
    std::shared_ptr<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

    std::vector<const Core::Elements::Element*> element_ptr(2);

    element_ptr[0] = elepairptr->element1();
    element_ptr[1] = elepairptr->element2();

    // element Dof values relevant for centerline interpolation
    std::vector<std::vector<double>> element_posdofvec_absolutevalues(2);

    for (unsigned int ielement = 0; ielement < 2; ++ielement)
    {
      // extract the Dof values of this element from displacement vector
      BEAMINTERACTION::Utils::extract_pos_dof_vec_absolute_values(discret(), element_ptr[ielement],
          *beam_interaction_data_state_ptr()->get_dis_col_np(),
          element_posdofvec_absolutevalues[ielement]);
    }

    // update the Dof values in the interaction element pair object
    elepairptr->reset_state(g_state().get_time_np(), element_posdofvec_absolutevalues[0],
        element_posdofvec_absolutevalues[1]);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::evaluate_force()
{
  check_init_setup();

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::evaluate_force");

  // resulting discrete element force vectors of the two interacting elements
  std::vector<Core::LinAlg::SerialDenseVector> eleforce(2);

  // resulting discrete force vectors (centerline DOFs only!) of the two
  // interacting elements
  std::vector<Core::LinAlg::SerialDenseVector> eleforce_centerlineDOFs(2);

  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> dummystiff;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  // are non-zero force values returned which need assembly?
  bool pair_is_active = false;


  std::vector<std::shared_ptr<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
  {
    std::shared_ptr<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

    // conditions applied to the elements of this pair
    std::vector<Core::Conditions::Condition*> conditions_element1;
    std::vector<Core::Conditions::Condition*> conditions_element2;
    get_beam_potential_conditions_applied_to_this_element_pair(
        *elepairptr, conditions_element1, conditions_element2);

    for (auto& k : conditions_element1)
    {
      int npotlaw1 = k->parameters().get<int>("POTLAW");

      for (auto& j : conditions_element2)
      {
        int npotlaw2 = j->parameters().get<int>("POTLAW");

        if (npotlaw1 == npotlaw2 and npotlaw1 > 0)
        {
          std::vector<Core::Conditions::Condition*> currconds;
          currconds.clear();
          currconds.push_back(k);
          currconds.push_back(j);

          // be careful here, as npotlaw =1 corresponds to first entry of ki_/mi_, therefore index 0
          if (npotlaw1 > (int)beam_potential_params().potential_law_prefactors().size())
            FOUR_C_THROW(
                "number of potential law specified in line charge condition exceeds"
                " number of defined potential laws!");

          pair_is_active = elepairptr->evaluate(&(eleforce_centerlineDOFs[0]),
              &(eleforce_centerlineDOFs[1]), nullptr, nullptr, nullptr, nullptr, currconds,
              beam_potential_params().potential_law_prefactors().at(npotlaw1 - 1),
              beam_potential_params().potential_law_exponents().at(npotlaw1 - 1));

          // Todo make this more efficient by summing all contributions from one element pair
          //      before assembly and communication
          if (pair_is_active)
          {
            elegids[0] = elepairptr->element1()->id();
            elegids[1] = elepairptr->element2()->id();

            // assemble force vector affecting the centerline DoFs only
            // into element force vector ('all DoFs' format, as usual)
            BEAMINTERACTION::Utils::assemble_centerline_dof_force_stiff_into_element_force_stiff(
                discret(), elegids, eleforce_centerlineDOFs, dummystiff, &eleforce, nullptr);

            // assemble the contributions into force vector class variable
            // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
            BEAMINTERACTION::Utils::fe_assemble_ele_force_stiff_into_system_vector_matrix(discret(),
                elegids, eleforce, dummystiff, beam_interaction_data_state_ptr()->get_force_np(),
                nullptr);
          }
        }
      }
    }
  }
  return true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::evaluate_stiff()
{
  check_init_setup();

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::evaluate_stiff");

  // linearizations
  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> elestiff(
      2, std::vector<Core::LinAlg::SerialDenseMatrix>(2));

  // linearizations (centerline DOFs only!)
  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> elestiff_centerlineDOFs(
      2, std::vector<Core::LinAlg::SerialDenseMatrix>(2));

  std::vector<Core::LinAlg::SerialDenseVector> dummyforce;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  // are non-zero stiffness values returned which need assembly?
  bool pair_is_active = false;


  std::vector<std::shared_ptr<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
  {
    std::shared_ptr<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

    // conditions applied to the elements of this pair
    std::vector<Core::Conditions::Condition*> conditions_element1;
    std::vector<Core::Conditions::Condition*> conditions_element2;
    get_beam_potential_conditions_applied_to_this_element_pair(
        *elepairptr, conditions_element1, conditions_element2);

    for (unsigned int k = 0; k < conditions_element1.size(); ++k)
    {
      int npotlaw1 = conditions_element1[k]->parameters().get<int>("POTLAW");

      for (unsigned int j = 0; j < conditions_element2.size(); ++j)
      {
        int npotlaw2 = conditions_element2[j]->parameters().get<int>("POTLAW");

        if (npotlaw1 == npotlaw2 and npotlaw1 > 0)
        {
          std::vector<Core::Conditions::Condition*> currconds;
          currconds.clear();
          currconds.push_back(conditions_element1[k]);
          currconds.push_back(conditions_element2[j]);

          // be careful here, as npotlaw =1 corresponds to first entry of ki_/mi_, therefore index 0
          if (npotlaw1 > (int)beam_potential_params().potential_law_prefactors().size())
            FOUR_C_THROW(
                "number of potential law specified in line charge condition exceeds"
                " number of defined potential laws!");


          pair_is_active = elepairptr->evaluate(nullptr, nullptr, &(elestiff_centerlineDOFs[0][0]),
              &(elestiff_centerlineDOFs[0][1]), &(elestiff_centerlineDOFs[1][0]),
              &(elestiff_centerlineDOFs[1][1]), currconds,
              beam_potential_params().potential_law_prefactors().at(npotlaw1 - 1),
              beam_potential_params().potential_law_exponents().at(npotlaw1 - 1));

          // Todo make this more efficient by summing all contributions from one element pair
          //      before assembly and communication
          if (pair_is_active)
          {
            elegids[0] = elepairptr->element1()->id();
            elegids[1] = elepairptr->element2()->id();

            // assemble stiffness matrix affecting the centerline DoFs only
            // into element stiffness matrix ('all DoFs' format, as usual)
            BEAMINTERACTION::Utils::assemble_centerline_dof_force_stiff_into_element_force_stiff(
                discret(), elegids, dummyforce, elestiff_centerlineDOFs, nullptr, &elestiff);

            // assemble the contributions into force vector class variable
            // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
            BEAMINTERACTION::Utils::fe_assemble_ele_force_stiff_into_system_vector_matrix(discret(),
                elegids, dummyforce, elestiff, nullptr,
                beam_interaction_data_state_ptr()->get_stiff());
          }
        }
      }
    }
  }
  return true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::evaluate_force_stiff()
{
  check_init_setup();

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::evaluate_force_stiff");

  // resulting discrete element force vectors of the two interacting elements
  std::vector<Core::LinAlg::SerialDenseVector> eleforce(2);

  // resulting discrete force vectors (centerline DOFs only!) of the two
  // interacting elements
  std::vector<Core::LinAlg::SerialDenseVector> eleforce_centerlineDOFs(2);

  // linearizations
  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> elestiff(
      2, std::vector<Core::LinAlg::SerialDenseMatrix>(2));

  // linearizations (centerline DOFs only!)
  std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> elestiff_centerlineDOFs(
      2, std::vector<Core::LinAlg::SerialDenseMatrix>(2));

  // element gids of interacting elements
  std::vector<int> elegids(2);

  // are non-zero stiffness values returned which need assembly?
  bool pair_is_active = false;


  std::vector<std::shared_ptr<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
  {
    std::shared_ptr<BEAMINTERACTION::BeamPotentialPair> elepairptr = *iter;

    elegids[0] = elepairptr->element1()->id();
    elegids[1] = elepairptr->element2()->id();

    // conditions applied to the elements of this pair
    std::vector<Core::Conditions::Condition*> conditions_element1;
    std::vector<Core::Conditions::Condition*> conditions_element2;
    get_beam_potential_conditions_applied_to_this_element_pair(
        *elepairptr, conditions_element1, conditions_element2);

    for (unsigned int k = 0; k < conditions_element1.size(); ++k)
    {
      int npotlaw1 = conditions_element1[k]->parameters().get<int>("POTLAW");

      for (unsigned int j = 0; j < conditions_element2.size(); ++j)
      {
        int npotlaw2 = conditions_element2[j]->parameters().get<int>("POTLAW");

        if (npotlaw1 == npotlaw2 and npotlaw1 > 0)
        {
          std::vector<Core::Conditions::Condition*> currconds;
          currconds.clear();
          currconds.push_back(conditions_element1[k]);
          currconds.push_back(conditions_element2[j]);

          // be careful here, as npotlaw =1 corresponds to first entry of ki_/mi_, therefore index 0
          if (npotlaw1 > (int)beam_potential_params().potential_law_prefactors().size())
            FOUR_C_THROW(
                "number of potential law specified in line charge condition exceeds"
                " number of defined potential laws!");


          pair_is_active =
              elepairptr->evaluate(&(eleforce_centerlineDOFs[0]), &(eleforce_centerlineDOFs[1]),
                  &(elestiff_centerlineDOFs[0][0]), &(elestiff_centerlineDOFs[0][1]),
                  &(elestiff_centerlineDOFs[1][0]), &(elestiff_centerlineDOFs[1][1]), currconds,
                  beam_potential_params().potential_law_prefactors().at(npotlaw1 - 1),
                  beam_potential_params().potential_law_exponents().at(npotlaw1 - 1));

          // Todo make this more efficient by summing all contributions from one element pair
          //      before assembly and communication
          if (pair_is_active)
          {
            elegids[0] = elepairptr->element1()->id();
            elegids[1] = elepairptr->element2()->id();

            // assemble force vector and stiffness matrix affecting the centerline DoFs only
            // into element force vector and stiffness matrix ('all DoFs' format, as usual)
            BEAMINTERACTION::Utils::assemble_centerline_dof_force_stiff_into_element_force_stiff(
                discret(), elegids, eleforce_centerlineDOFs, elestiff_centerlineDOFs, &eleforce,
                &elestiff);

            // assemble the contributions into force vector class variable
            // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
            BEAMINTERACTION::Utils::fe_assemble_ele_force_stiff_into_system_vector_matrix(discret(),
                elegids, eleforce, elestiff, beam_interaction_data_state_ptr()->get_force_np(),
                beam_interaction_data_state_ptr()->get_stiff());
          }
        }
      }
    }
  }

  //  print_active_beam_potential_set(std::cout);

  return true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::update_step_state(const double& timefac_n)
{
  check_init_setup();

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::pre_update_step_element(bool beam_redist)
{
  check_init_setup();

  /* Fixme
   * writing vtk output needs to be done BEFORE updating (and thus clearing
   * element pairs)
   * move this to runtime_output_step_state as soon as we keep element pairs
   * from previous time step */
  if (visualization_manager_ != nullptr and
      g_state().get_step_n() % beam_potential_params()
                                   .get_beam_potential_visualization_output_params()
                                   ->output_interval_in_steps() ==
          0)
  {
    write_time_step_output_runtime_beam_potential();
  }

  // not repartition of binning discretization necessary
  return false;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::update_step_element(
    bool repartition_was_done)
{
  check_init_setup();

  nearby_elements_map_.clear();
  find_and_store_neighboring_elements();
  create_beam_potential_element_pairs();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::post_update_step_element()
{
  check_init_setup();

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
std::map<Solid::EnergyType, double> BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::get_energy()
    const
{
  check_init_setup();

  std::map<Solid::EnergyType, double> beam_interaction_potential;

  for (auto& elepairptr : beam_potential_element_pairs_)
  {
    beam_interaction_potential[Solid::beam_interaction_potential] += elepairptr->get_energy();
  }

  return beam_interaction_potential;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::output_step_state(
    Core::IO::DiscretizationWriter& iowriter) const
{
  check_init_setup();
  // nothing to do (so far)
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::runtime_output_step_state() const
{
  check_init_setup();
  // nothing to do (so far)
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::reset_step_state() { check_init_setup(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::write_restart(
    Core::IO::DiscretizationWriter& ia_writer, Core::IO::DiscretizationWriter& bin_writer) const
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::pre_read_restart()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::read_restart(
    Core::IO::DiscretizationReader& ia_reader)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::post_read_restart()
{
  check_init_setup();

  if (beam_potential_params().potential_reduction_length() != -1.0)
    setup_potential_reduction_strategy();

  nearby_elements_map_.clear();
  find_and_store_neighboring_elements();
  create_beam_potential_element_pairs();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::run_post_iterate(
    const ::NOX::Solver::Generic& solver)
{
  check_init_setup();

  if (visualization_manager_ != nullptr and beam_potential_params()
                                                .get_beam_potential_visualization_output_params()
                                                ->output_every_iteration())
  {
    write_iteration_output_runtime_beam_potential(solver.getNumIterations());
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::add_bins_to_bin_col_map(
    std::set<int>& colbins)
{
  check_init_setup();
  // nothing to do
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    add_bins_with_relevant_content_for_ia_discret_col_map(std::set<int>& colbins) const
{
  check_init_setup();
  // nothing to do
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::get_half_interaction_distance(
    double& half_interaction_distance)
{
  check_init_setup();

  if (beam_potential_params().cutoff_radius() > 0.0)
  {
    half_interaction_distance = 0.5 * beam_potential_params().cutoff_radius();

    if (g_state().get_my_rank() == 0)
      Core::IO::cout(Core::IO::verbose) << " beam potential half interaction distance "
                                        << half_interaction_distance << Core::IO::endl;
  }
  else
  {
    FOUR_C_THROW(
        "You have to set a cutoff radius for beam-to-? potential-based interactions in order "
        "to use REPARTITIONSTRATEGY = Adaptive!");
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::find_and_store_neighboring_elements()
{
  check_init();

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::find_and_store_neighboring_elements");


  const auto geometric_search_params_ptr_ = Core::GeometricSearch::GeometricSearchParams(
      Global::Problem::instance()->geometric_search_params(),
      Global::Problem::instance()->io_params());

  // Get vector of all beam element bounding boxes.
  int const numroweles = ele_type_map_extractor_ptr()->beam_map()->NumMyElements();
  std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> beam_bounding_boxes;
  for (int rowele_i = 0; rowele_i < numroweles; ++rowele_i)
  {
    int const elegid = ele_type_map_extractor_ptr()->beam_map()->GID(rowele_i);
    Core::Elements::Element* currele = discret().g_element(elegid);

    beam_bounding_boxes.emplace_back(std::make_pair(elegid,
        currele->get_bounding_volume(discret(),
            *beam_interaction_data_state_ptr()->get_dis_col_np(), geometric_search_params_ptr_)));
  }

  // Get vector of the bounding boxes of all possible interacting elements (also including beams
  // if beam-to-beam contact is activated).
  int const numcoleles = discret().num_my_col_elements();
  std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> other_bounding_boxes;
  for (int colele_i = 0; colele_i < numcoleles; ++colele_i)
  {
    // Check if the current element is relevant for beam-to-xxx contact.
    Core::Elements::Element* currele = discret().l_col_element(colele_i);
    const Core::Binstrategy::Utils::BinContentType contact_type =
        BEAMINTERACTION::Utils::convert_element_to_bin_content_type(currele);

    {
      other_bounding_boxes.emplace_back(std::make_pair(currele->id(),
          currele->get_bounding_volume(discret(),
              *beam_interaction_data_state_ptr()->get_dis_col_np(), geometric_search_params_ptr_)));
    }
  }

  // Get colliding pairs.
  const auto& [indices, offsets] = collision_search(other_bounding_boxes, beam_bounding_boxes,
      discret().get_comm(), geometric_search_params_ptr_.verbosity_);

  // Create the beam-to-xxx pair pointers according to the search.
  for (size_t i_beam = 0; i_beam < beam_bounding_boxes.size(); i_beam++)
  {
    const int beam_gid = beam_bounding_boxes[i_beam].first;
    Core::Elements::Element* currele = discret().g_element(beam_gid);
    std::set<Core::Elements::Element*> neighboring_elements;
    for (int j = offsets[i_beam]; j < offsets[i_beam + 1]; j++)
    {
      neighboring_elements.insert(discret().g_element(other_bounding_boxes[indices[j]].first));
    }

    // sort out elements that should not be considered in contact evaluation
    select_eles_to_be_considered_for_potential_evaluation(currele, neighboring_elements);

    nearby_elements_map_[beam_gid] = neighboring_elements;
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    select_eles_to_be_considered_for_potential_evaluation(
        Core::Elements::Element* currele, std::set<Core::Elements::Element*>& neighbors) const
{
  check_init();

  // sort out elements that should not be considered in potential evaluation
  std::set<Core::Elements::Element*>::iterator eiter;
  for (eiter = neighbors.begin(); eiter != neighbors.end();)
  {
    bool toerase = false;

    Core::Elements::Element* currneighborele = *eiter;

    // 1) ensure each interaction is only evaluated once (keep in mind that we are
    //    using FEMatrices and FEvectors -> || (*eiter)->Owner() != myrank not necessary)
    if (not(currele->id() < currneighborele->id()))
    {
      toerase = true;
    }

    // 2) exclude "self-interaction", i.e. a pair of elements on the same physical beam
    //    TODO introduce flag for self-interaction in input file
    else
    {
      // get the conditions applied to both elements of the pair and decide whether they need to be
      // evaluated
      std::vector<Core::Conditions::Condition*> conds1, conds2;

      // since only the nodes know about their conditions, we need this workaround
      // we assume that a linecharge condition is always applied to the entire physical beam, i.e.
      // it is sufficient to check only one node
      Core::Nodes::Node** nodes1;
      Core::Nodes::Node** nodes2;
      nodes1 = currele->nodes();
      nodes2 = currneighborele->nodes();

      FOUR_C_ASSERT(nodes1 != nullptr and nodes2 != nullptr, "pointer to nodes is nullptr!");
      FOUR_C_ASSERT(nodes1[0] != nullptr and nodes2[0] != nullptr, "pointer to nodes is nullptr!");

      nodes1[0]->get_condition("BeamPotentialLineCharge", conds1);

      // get correct condition for beam or rigid sphere element
      if (BEAMINTERACTION::Utils::is_beam_element(*currneighborele))
        nodes2[0]->get_condition("BeamPotentialLineCharge", conds2);
      else if (BEAMINTERACTION::Utils::is_rigid_sphere_element(*currneighborele))
        nodes2[0]->get_condition("RigidspherePotentialPointCharge", conds2);
      else
        FOUR_C_THROW(
            "Only beam-to-beampotential or beam-to-sphere -based interaction is implemented yet. "
            "No other types of elements allowed!");

      // validinteraction == true includes: both eles "loaded" by a charge condition of same
      // potential law
      bool validinteraction = false;

      for (unsigned int i = 0; i < conds1.size(); ++i)
      {
        int npotlaw1 = conds1[i]->parameters().get<int>("POTLAW");

        for (unsigned int j = 0; j < conds2.size(); ++j)
        {
          int npotlaw2 = conds2[j]->parameters().get<int>("POTLAW");

          // here, we also exclude "self-interaction", i.e. a pair of elements on the same physical
          // beam
          // TODO introduce flag for self-interaction in input file
          if (conds1[i] != conds2[j] and npotlaw1 == npotlaw2) validinteraction = true;
        }
      }

      if (not validinteraction) toerase = true;
    }


    if (toerase)
      neighbors.erase(eiter++);
    else
      ++eiter;
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::setup_potential_reduction_strategy()
{
  LengthToEdgeImplementation::DataMaps data_maps;

  // get element data on current proc
  for (int rowele_i = 0; rowele_i < ele_type_map_extractor_ptr()->beam_map()->NumMyElements();
       ++rowele_i)
  {
    const int ele_gid = ele_type_map_extractor_ptr()->beam_map()->GID(rowele_i);
    Core::Elements::Element* ele_ptr = discret_ptr()->g_element(ele_gid);

    data_maps.ele_gid_length_map.insert(std::make_pair(
        ele_gid, dynamic_cast<Discret::Elements::Beam3Base*>(ele_ptr)->ref_length()));

    int left_node_gid = *ele_ptr->node_ids();
    // n_right is the local node-ID of the elements right node (at xi = 1) whereas the elements left
    // node (at xi = -1) allways has the local ID 1
    const int n_right = (ele_ptr->num_node() == 2) ? 1 : (ele_ptr->num_node() - 2);
    int right_node_gid = *(ele_ptr->node_ids() + n_right);

    data_maps.ele_gid_left_node_gid_map.insert(std::make_pair(ele_gid, left_node_gid));
    data_maps.ele_gid_right_node_gid_map.insert(std::make_pair(ele_gid, right_node_gid));
    data_maps.left_node_gid_ele_gid_map.insert(std::make_pair(left_node_gid, ele_gid));
    data_maps.right_node_gid_ele_gid_map.insert(std::make_pair(right_node_gid, ele_gid));
  }

  // broadcast all data maps to all procs
  data_maps.ele_gid_length_map =
      Core::Communication::all_gather(data_maps.ele_gid_length_map, discret().get_comm());
  data_maps.ele_gid_left_node_gid_map =
      Core::Communication::all_gather(data_maps.ele_gid_left_node_gid_map, discret().get_comm());
  data_maps.ele_gid_right_node_gid_map =
      Core::Communication::all_gather(data_maps.ele_gid_right_node_gid_map, discret().get_comm());
  data_maps.left_node_gid_ele_gid_map =
      Core::Communication::all_gather(data_maps.left_node_gid_ele_gid_map, discret().get_comm());
  data_maps.right_node_gid_ele_gid_map =
      Core::Communication::all_gather(data_maps.right_node_gid_ele_gid_map, discret().get_comm());

  // determine length to edge for each element and add to map
  for (const auto& [ele_gid, _] : data_maps.ele_gid_length_map)
  {
    beam_potential_params().ele_gid_prior_length_map_.insert(std::make_pair(
        ele_gid, std::make_pair(LengthToEdgeImplementation::determine_length_to_edge(data_maps,
                                    ele_gid, data_maps.ele_gid_left_node_gid_map.at(ele_gid)),
                     LengthToEdgeImplementation::determine_length_to_edge(
                         data_maps, ele_gid, data_maps.ele_gid_right_node_gid_map.at(ele_gid)))));
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::create_beam_potential_element_pairs()
{
  // Todo maybe keep existing pairs and reuse them ?
  beam_potential_element_pairs_.clear();

  std::map<int, std::set<Core::Elements::Element*>>::const_iterator nearbyeleiter;

  for (nearbyeleiter = nearby_elements_map_.begin(); nearbyeleiter != nearby_elements_map_.end();
       ++nearbyeleiter)
  {
    const int elegid = nearbyeleiter->first;
    std::vector<Core::Elements::Element const*> ele_ptrs(2);
    ele_ptrs[0] = discret_ptr()->g_element(elegid);

    std::set<Core::Elements::Element*>::const_iterator secondeleiter;
    for (secondeleiter = nearbyeleiter->second.begin();
         secondeleiter != nearbyeleiter->second.end(); ++secondeleiter)
    {
      ele_ptrs[1] = *secondeleiter;

      std::shared_ptr<BEAMINTERACTION::BeamPotentialPair> newbeaminteractionpair =
          BEAMINTERACTION::BeamPotentialPair::create(ele_ptrs, beam_potential_params());

      newbeaminteractionpair->init(beam_potential_params_ptr(), ele_ptrs[0], ele_ptrs[1]);

      newbeaminteractionpair->setup();

      beam_potential_element_pairs_.push_back(newbeaminteractionpair);
    }
  }

  if (static_cast<int>(beam_potential_element_pairs_.size()) > 0)
  {
    Core::IO::cout(Core::IO::standard)
        << "PID " << std::setw(2) << std::right << g_state().get_my_rank() << " currently monitors "
        << std::setw(5) << std::right << beam_potential_element_pairs_.size()
        << " beam potential pairs" << Core::IO::endl;
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::print_all_beam_potential_element_pairs(
    std::ostream& out) const
{
  out << "\n\nCurrent BeamPotentialElementPairs: ";
  std::vector<std::shared_ptr<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
    (*iter)->print(out);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::print_active_beam_potential_set(
    std::ostream& out) const
{
  // Todo
  out << "\n    Active BeamToBeam Potential Set (PID " << g_state().get_my_rank()
      << "):-----------------------------------------\n";
  out << "    ID1            ID2              T xi       eta      angle    gap         force\n";

  std::vector<std::shared_ptr<BEAMINTERACTION::BeamPotentialPair>>::const_iterator iter;
  for (iter = beam_potential_element_pairs_.begin(); iter != beam_potential_element_pairs_.end();
       ++iter)
    (*iter)->print_summary_one_line_per_active_segment_pair(out);

  out << std::endl;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    get_beam_potential_conditions_applied_to_this_element_pair(
        BEAMINTERACTION::BeamPotentialPair const& elementpair,
        std::vector<Core::Conditions::Condition*>& conditions_element1,
        std::vector<Core::Conditions::Condition*>& conditions_element2) const
{
  // since only the nodes know about their conditions, we need this workaround
  // we assume that a linecharge condition is always applied to the entire physical beam, i.e. it is
  // sufficient to check only one node
  const Core::Elements::Element* ele1 = elementpair.element1();
  const Core::Elements::Element* ele2 = elementpair.element2();

  const Core::Nodes::Node* const* nodes1;
  const Core::Nodes::Node* const* nodes2;
  nodes1 = ele1->nodes();
  nodes2 = ele2->nodes();

  FOUR_C_ASSERT(nodes1 != nullptr and nodes2 != nullptr, "pointer to nodes is nullptr!");
  FOUR_C_ASSERT(nodes1[0] != nullptr and nodes2[0] != nullptr, "pointer to nodes is nullptr!");

  nodes1[0]->get_condition("BeamPotentialLineCharge", conditions_element1);

  // get correct condition for beam or rigid sphere element
  if (BEAMINTERACTION::Utils::is_beam_element(*ele2))
    nodes2[0]->get_condition("BeamPotentialLineCharge", conditions_element2);
  else if (BEAMINTERACTION::Utils::is_rigid_sphere_element(*ele2))
    nodes2[0]->get_condition("RigidspherePotentialPointCharge", conditions_element2);
  else
    FOUR_C_THROW(
        "Only beam-to-beam or beam-to-sphere potential-based interaction is implemented yet. "
        "No other types of elements allowed!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::print_console_welcome_message(
    std::ostream& out) const
{
  // console welcome message
  if (g_state().get_my_rank() == 0)
  {
    std::cout << "=============== Beam Potential-Based Interaction ===============" << std::endl;

    switch (beam_potential_params().potential_type())
    {
      case Inpar::BEAMPOTENTIAL::beampot_surf:
      {
        std::cout << "Potential Type:      Surface" << std::endl;
        break;
      }
      case Inpar::BEAMPOTENTIAL::beampot_vol:
      {
        std::cout << "Potential Type:      Volume" << std::endl;
        break;
      }
      default:
        FOUR_C_THROW("Potential type not supported!");
    }

    std::cout << "Potential Law:       Phi(r) = ";
    for (unsigned int isummand = 0;
         isummand < beam_potential_params().potential_law_prefactors().size(); ++isummand)
    {
      if (isummand > 0) std::cout << " + ";

      std::cout << "(" << beam_potential_params().potential_law_prefactors().at(isummand)
                << ") * r^(-" << beam_potential_params().potential_law_exponents().at(isummand)
                << ")";
    }
    std::cout << std::endl;

    std::cout << "================================================================\n" << std::endl;
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::init_output_runtime_beam_potential()
{
  check_init();

  visualization_manager_ = std::make_shared<Core::IO::VisualizationManager>(
      beam_potential_params()
          .get_beam_potential_visualization_output_params()
          ->get_visualization_parameters(),
      discret().get_comm(), "beam-potential");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    write_time_step_output_runtime_beam_potential() const
{
  check_init_setup();

  auto [output_time, output_step] = Core::IO::get_time_and_time_step_index_for_output(
      beam_potential_params()
          .get_beam_potential_visualization_output_params()
          ->get_visualization_parameters(),
      g_state().get_time_n(), g_state().get_step_n());
  write_output_runtime_beam_potential(output_step, output_time);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::
    write_iteration_output_runtime_beam_potential(int iteration_number) const
{
  check_init_setup();

  auto [output_time, output_step] = Core::IO::get_time_and_time_step_index_for_output(
      beam_potential_params()
          .get_beam_potential_visualization_output_params()
          ->get_visualization_parameters(),
      g_state().get_time_n(), g_state().get_step_n(), iteration_number);
  write_output_runtime_beam_potential(output_step, output_time);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamPotential::write_output_runtime_beam_potential(
    int timestep_number, double time) const
{
  check_init_setup();

  const unsigned int num_spatial_dimensions = 3;

  // estimate for number of interacting Gauss points = number of row points for writer object
  unsigned int num_row_points = 0;

  if (beam_potential_params()
          .get_beam_potential_visualization_output_params()
          ->is_write_forces_moments_per_element_pair())
  {
    num_row_points = 2 * beam_potential_element_pairs_.size() *
                     beam_potential_params().number_integration_segments() *
                     beam_potential_params().number_gauss_points();
  }
  else
  {
    // Todo: this won't perfectly work in parallel yet since some communication would be required
    //    if ( global_state().get_my_rank() != 0 )
    //      FOUR_C_THROW("visualization of resulting forces not implemented in parallel yet!");

    num_row_points = discret().num_global_elements() *
                     beam_potential_params().number_integration_segments() *
                     beam_potential_params().number_gauss_points();
  }

  /* Note: - each UID set is not unique due to the fact that each GP produces two force
   *         vectors (one on the slave side, one on the master side)
   *       - in case of the single length specific approach (SBIP) the uid for the GP refers
   *         to the slave beam element */
  // get and prepare storage for uid_0_beam_1_gid values
  std::vector<int> uid_0_beam_1_gid(0);
  uid_0_beam_1_gid.reserve(num_row_points);

  // get and prepare storage for uid_1_beam_2_gid values
  std::vector<int> uid_1_beam_2_gid(0);
  uid_1_beam_2_gid.reserve(num_row_points);

  // get and prepare storage for uid_2_gp_id values
  std::vector<int> uid_2_gp_id(0);
  uid_2_gp_id.reserve(num_row_points);

  // get and prepare storage for point coordinate values
  auto& visualization_data = visualization_manager_->get_visualization_data();
  std::vector<double>& point_coordinates = visualization_data.get_point_coordinates();
  point_coordinates.clear();
  point_coordinates.reserve(num_spatial_dimensions * num_row_points);

  // force values: collect data and append to visualization results if desired
  std::vector<double> potential_force_vector(0);
  potential_force_vector.reserve(num_spatial_dimensions * num_row_points);

  // moment values: collect data and append to visualization results if desired
  std::vector<double> potential_moment_vector(0);
  potential_moment_vector.reserve(num_spatial_dimensions * num_row_points);


  // loop over my points and collect the geometry/grid data, i.e. interacting points
  std::vector<Core::LinAlg::Matrix<3, 1, double>> coordinates_ele1_this_pair;
  std::vector<Core::LinAlg::Matrix<3, 1, double>> coordinates_ele2_this_pair;

  std::vector<Core::LinAlg::Matrix<3, 1, double>> potential_forces_ele1_this_pair;
  std::vector<Core::LinAlg::Matrix<3, 1, double>> potential_forces_ele2_this_pair;

  std::vector<Core::LinAlg::Matrix<3, 1, double>> potential_moments_ele1_this_pair;
  std::vector<Core::LinAlg::Matrix<3, 1, double>> potential_moments_ele2_this_pair;

  // loop over contact pairs and retrieve all active contact point coordinates
  std::vector<std::shared_ptr<BEAMINTERACTION::BeamPotentialPair>>::const_iterator pair_iter;
  for (pair_iter = beam_potential_element_pairs_.begin();
       pair_iter != beam_potential_element_pairs_.end(); ++pair_iter)
  {
    // retrieve data for interacting points of element 1 and element 2
    (*pair_iter)->get_all_interacting_point_coords_element1(coordinates_ele1_this_pair);
    (*pair_iter)->get_all_interacting_point_coords_element2(coordinates_ele2_this_pair);
    (*pair_iter)->get_forces_at_all_interacting_points_element1(potential_forces_ele1_this_pair);
    (*pair_iter)->get_forces_at_all_interacting_points_element2(potential_forces_ele2_this_pair);
    (*pair_iter)->get_moments_at_all_interacting_points_element1(potential_moments_ele1_this_pair);
    (*pair_iter)->get_moments_at_all_interacting_points_element2(potential_moments_ele2_this_pair);

    const unsigned int num_interacting_points_per_element =
        (unsigned int)coordinates_ele1_this_pair.size();

    FOUR_C_ASSERT(
        num_interacting_points_per_element == (unsigned int)coordinates_ele2_this_pair.size(),
        "number of interacting points on element 1 does not match number of interacting points "
        "on element 2!");

    FOUR_C_ASSERT(
        num_interacting_points_per_element == (unsigned int)potential_forces_ele1_this_pair.size(),
        "number of interacting points on element 1 does not match number of potential forces!");

    FOUR_C_ASSERT(
        num_interacting_points_per_element == (unsigned int)potential_forces_ele2_this_pair.size(),
        "number of interacting points on element 2 does not match number of potential forces!");


    for (unsigned int ipoint = 0; ipoint < num_interacting_points_per_element; ++ipoint)
    {
      // ignore point pairs with zero forces
      /* (e.g. if no valid point-to-curve projection in master-slave approach or
       * contribution is neglected on element pair level due to cutoff value) */
      if (potential_forces_ele1_this_pair[ipoint].norm2() < 1e-16 and
          potential_forces_ele2_this_pair[ipoint].norm2() < 1e-16 and
          potential_moments_ele1_this_pair[ipoint].norm2() < 1e-16 and
          potential_moments_ele2_this_pair[ipoint].norm2() < 1e-16)
      {
        continue;
      }


      // this is easier, since data is computed and stored in this 'element-pairwise' format
      if (beam_potential_params()
              .get_beam_potential_visualization_output_params()
              ->is_write_forces_moments_per_element_pair())
      {
        uid_0_beam_1_gid.push_back((*pair_iter)->element1()->id());
        uid_1_beam_2_gid.push_back((*pair_iter)->element2()->id());
        uid_2_gp_id.push_back(ipoint);

        for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
        {
          point_coordinates.push_back(coordinates_ele1_this_pair[ipoint](idim));

          potential_force_vector.push_back(potential_forces_ele1_this_pair[ipoint](idim));
          potential_moment_vector.push_back(potential_moments_ele1_this_pair[ipoint](idim));
        }

        uid_0_beam_1_gid.push_back((*pair_iter)->element1()->id());
        uid_1_beam_2_gid.push_back((*pair_iter)->element2()->id());
        uid_2_gp_id.push_back(ipoint);

        for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
        {
          point_coordinates.push_back(coordinates_ele2_this_pair[ipoint](idim));

          potential_force_vector.push_back(potential_forces_ele2_this_pair[ipoint](idim));
          potential_moment_vector.push_back(potential_moments_ele2_this_pair[ipoint](idim));
        }
      }
      /* in this case, we need to identify unique Gauss points based on their coordinate values and
       * compute resulting force/moment at this point by summation of contributions from all element
       * pairs */
      else
      {
        // interacting point on first element
        std::vector<double>::iterator xcoord_iter = point_coordinates.begin();

        // try to find data point with identical coordinates
        while (point_coordinates.size() >= 3 and xcoord_iter != point_coordinates.end() - 2)
        {
          // find identical x-coordinate value
          xcoord_iter = std::find(
              xcoord_iter, point_coordinates.end() - 2, coordinates_ele1_this_pair[ipoint](0));

          // check whether we've reached the end -> no match
          if (xcoord_iter == point_coordinates.end() - 2)
          {
            break;
          }
          // we have a match -> check whether also y- and z-coordinate value are identical
          else if ((xcoord_iter - point_coordinates.begin()) % 3 == 0 and
                   *(xcoord_iter + 1) == coordinates_ele1_this_pair[ipoint](1) and
                   *(xcoord_iter + 2) == coordinates_ele1_this_pair[ipoint](2))
          {
            int offset = xcoord_iter - point_coordinates.begin();

            for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
            {
              *(potential_force_vector.begin() + offset + idim) +=
                  potential_forces_ele1_this_pair[ipoint](idim);
              *(potential_moment_vector.begin() + offset + idim) +=
                  potential_moments_ele1_this_pair[ipoint](idim);
            }

            break;
          }
          // we have a matching value but it's not a point with the identical (x,y,z) coordinates
          else
          {
            xcoord_iter++;
          }
        }

        // add as a new point if not found above
        if (xcoord_iter == point_coordinates.end() - 2 or point_coordinates.empty())
        {
          uid_0_beam_1_gid.push_back((*pair_iter)->element1()->id());
          uid_1_beam_2_gid.push_back((*pair_iter)->element2()->id());
          uid_2_gp_id.push_back(ipoint);

          for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
          {
            point_coordinates.push_back(coordinates_ele1_this_pair[ipoint](idim));

            potential_force_vector.push_back(potential_forces_ele1_this_pair[ipoint](idim));
            potential_moment_vector.push_back(potential_moments_ele1_this_pair[ipoint](idim));
          }
        }


        // interacting point on second element
        xcoord_iter = point_coordinates.begin();

        // try to find data point with identical coordinates
        while (xcoord_iter != point_coordinates.end() - 2)
        {
          // find identical x-coordinate value
          xcoord_iter = std::find(
              xcoord_iter, point_coordinates.end() - 2, coordinates_ele2_this_pair[ipoint](0));

          // check whether we've reached the end -> no match
          if (xcoord_iter == point_coordinates.end() - 2)
          {
            break;
          }
          // we have a match -> check whether also y- and z-coordinate value are identical
          else if ((xcoord_iter - point_coordinates.begin()) % 3 == 0 and
                   *(xcoord_iter + 1) == coordinates_ele2_this_pair[ipoint](1) and
                   *(xcoord_iter + 2) == coordinates_ele2_this_pair[ipoint](2))
          {
            int offset = xcoord_iter - point_coordinates.begin();

            for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
            {
              *(potential_force_vector.begin() + offset + idim) +=
                  potential_forces_ele2_this_pair[ipoint](idim);
              *(potential_moment_vector.begin() + offset + idim) +=
                  potential_moments_ele2_this_pair[ipoint](idim);
            }

            break;
          }
          // we have a matching value but it's not a point with the identical (x,y,z) coordinates
          else
          {
            xcoord_iter++;
          }
        }

        // add as a new point if not found above
        if (xcoord_iter == point_coordinates.end() - 2)
        {
          uid_0_beam_1_gid.push_back((*pair_iter)->element1()->id());
          uid_1_beam_2_gid.push_back((*pair_iter)->element2()->id());
          uid_2_gp_id.push_back(ipoint);

          for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
          {
            point_coordinates.push_back(coordinates_ele2_this_pair[ipoint](idim));

            potential_force_vector.push_back(potential_forces_ele2_this_pair[ipoint](idim));
            potential_moment_vector.push_back(potential_moments_ele2_this_pair[ipoint](idim));
          }
        }
      }
    }
  }


  // append all desired output data to the writer object's storage

  if (beam_potential_params().get_beam_potential_visualization_output_params()->is_write_forces())
  {
    visualization_manager_->get_visualization_data().set_point_data_vector(
        "force", potential_force_vector, num_spatial_dimensions);
  }

  if (beam_potential_params().get_beam_potential_visualization_output_params()->is_write_moments())
  {
    visualization_manager_->get_visualization_data().set_point_data_vector(
        "moment", potential_moment_vector, num_spatial_dimensions);
  }

  if (beam_potential_params().get_beam_potential_visualization_output_params()->is_write_uids())
  {
    visualization_manager_->get_visualization_data().set_point_data_vector(
        "uid_0_beam_1_gid", uid_0_beam_1_gid, 1);

    visualization_manager_->get_visualization_data().set_point_data_vector(
        "uid_1_beam_2_gid", uid_1_beam_2_gid, 1);

    visualization_manager_->get_visualization_data().set_point_data_vector(
        "uid_2_gp_id", uid_2_gp_id, 1);
  }

  // finalize everything and write all required vtk files to filesystem
  visualization_manager_->write_to_disk(time, timestep_number);
}

FOUR_C_NAMESPACE_CLOSE
