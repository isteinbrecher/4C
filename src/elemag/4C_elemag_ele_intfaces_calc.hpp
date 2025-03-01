// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ELEMAG_ELE_INTFACES_CALC_HPP
#define FOUR_C_ELEMAG_ELE_INTFACES_CALC_HPP


#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class SparseMatrix;
}

namespace Core::FE
{
  class Discretization;
  class DiscretizationFaces;
}  // namespace Core::FE

namespace Discret
{
  namespace Elements
  {
    class ElemagIntFace;
    class ElemagEleParameter;
    class ElemagEleParameterTimInt;

    /// Interface base class for ElemagIntFaceImpl
    /*!
      This class exists to provide a common interface for all template
      versions of ElemagIntFaceImpl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of ElemagIntFaceImpl.
     */
    class ElemagIntFaceImplInterface
    {
     public:
      /// Empty constructor
      ElemagIntFaceImplInterface() {}

      /// Empty destructor
      virtual ~ElemagIntFaceImplInterface() = default;
      //! Assemble internal faces integrals using data from both parent elements
      virtual void assemble_internal_faces_using_neighbor_data(
          Discret::Elements::ElemagIntFace* intface,      ///< internal face element
          std::vector<int>& nds_master,                   ///< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,                    ///< nodal dofset w.r.t. slave element
          Teuchos::ParameterList& params,                 ///< parameter list
          Core::FE::DiscretizationFaces& discretization,  ///< faces discretization
          std::shared_ptr<Core::LinAlg::SparseMatrix> systemmatrix,   ///< systemmatrix
          std::shared_ptr<Core::LinAlg::Vector<double>> systemvector  ///< systemvector
          ) = 0;

      //! Evaluate internal faces
      virtual int evaluate_internal_faces(
          Discret::Elements::ElemagIntFace* intface,  ///< internal face element
          Teuchos::ParameterList& params,             ///< parameter list
          Core::FE::Discretization& discretization,   ///< discretization
          std::vector<int>& patchlm,                  ///< patch local map
          std::vector<int>& lm_masterToPatch,         ///< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,          ///< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,           ///< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          std::vector<Core::LinAlg::SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
          std::vector<Core::LinAlg::SerialDenseVector>& elevec_blocks   ///< element vector blocks
          ) = 0;


      /// Internal implementation class for ElemagIntFace elements (the first object is created in
      /// Discret::Elements::ElemagIntFace::Evaluate)
      static ElemagIntFaceImplInterface* impl(const Core::Elements::Element* ele);
    };

    /// Internal ElemagIntFace element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the ElemagIntFace element.

      <h3>Purpose</h3>

      The ElemagIntFace element will allocate exactly one object of this class
      for all ElemagIntFace elements with the same number of nodes in the mesh.
      This allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      (see fluid_ele_intfaces_calc.H)

    */
    template <Core::FE::CellType distype>
    class ElemagIntFaceImpl : public ElemagIntFaceImplInterface
    {
      friend class ElemagEleParameterTimInt;
      friend class ElemagEleParameterStd;

     public:
      /// Singleton access method
      static ElemagIntFaceImpl<distype>* instance(
          Core::Utils::SingletonAction action = Core::Utils::SingletonAction::create);

      /// Constructor
      ElemagIntFaceImpl();


      //! Assemble internal faces integrals using data from both parent elements
      void assemble_internal_faces_using_neighbor_data(
          Discret::Elements::ElemagIntFace* intface,      ///< internal face element
          std::vector<int>& nds_master,                   ///< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,                    ///< nodal dofset w.r.t. slave element
          Teuchos::ParameterList& params,                 ///< parameter list
          Core::FE::DiscretizationFaces& discretization,  ///< faces discretization
          std::shared_ptr<Core::LinAlg::SparseMatrix> systemmatrix,   ///< systemmatrix
          std::shared_ptr<Core::LinAlg::Vector<double>> systemvector  ///< systemvector
          ) override;

      //! Evaluate internal faces
      int evaluate_internal_faces(
          Discret::Elements::ElemagIntFace* intface,  ///< internal face element
          Teuchos::ParameterList& params,             ///< parameter list
          Core::FE::Discretization& discretization,   ///< discretization
          std::vector<int>& patchlm,                  ///< patch local map
          std::vector<int>& lm_masterToPatch,         ///< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,          ///< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,           ///< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          std::vector<Core::LinAlg::SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
          std::vector<Core::LinAlg::SerialDenseVector>& elevec_blocks   ///< element vector blocks
          ) override;

      //! decide which terms have to be assembled and decide the assembly pattern, return if no
      //! assembly required
      bool prepare_assemble(Teuchos::ParameterList& stabparams, Teuchos::ParameterList& faceparams);

     private:
      //! pointer to parameter lists
      Discret::Elements::ElemagEleParameter* elemagpara_;
      //! pointer to parameter list for time integration
      Discret::Elements::ElemagEleParameterTimInt* elemagparatimint_;

    };  // end class ElemagIntFaceImpl

  }  // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
