// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SO3_SURFACE_HPP
#define FOUR_C_SO3_SURFACE_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_elementtype.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret
{
  namespace Elements
  {
    class StructuralSurfaceType : public Core::Elements::ElementType
    {
     public:
      [[nodiscard]] std::string name() const override { return "StructuralSurfaceType"; }

      static StructuralSurfaceType& instance();

      Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

      std::shared_ptr<Core::Elements::Element> create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      Core::LinAlg::SerialDenseMatrix compute_null_space(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        Core::LinAlg::SerialDenseMatrix nullspace;
        FOUR_C_THROW("method ComputeNullSpace not implemented!");
      }

     private:
      static StructuralSurfaceType instance_;
    };

    /*!
    \brief A general structural surface element (tri3/tri6/quad4/quad8/quad9/nurbs9)

    This element implements all common shapes. It also implements all
    common types of loads currently demanded in 4C

    */
    class StructuralSurface : public Core::Elements::FaceElement
    {
     public:
      //! @name Friends
      friend class StructuralSurfaceType;

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner: Processor owning this surface
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param parent: The parent fluid element of this surface
      \param lsurface: the local surface number of this surface w.r.t. the parent element
      */
      explicit StructuralSurface(int id, int owner, int nnode, const int* nodeids,
          Core::Nodes::Node** nodes, Core::Elements::Element* parent, const int lsurface);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      explicit StructuralSurface(const StructuralSurface& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      [[nodiscard]] Core::Elements::Element* clone() const override;

      /*!
      \brief Get shape type of element
      */
      [[nodiscard]] Core::FE::CellType shape() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      [[nodiscard]] inline int unique_par_object_id() const override
      {
        return StructuralSurfaceType::instance().unique_par_object_id();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element

      */
      void unpack(Core::Communication::UnpackBuffer& buffer) override;

      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      [[nodiscard]] inline int num_dof_per_node(const Core::Nodes::Node& node) const override
      {
        return numdofpernode_;
      }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual Core::Elements::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      [[nodiscard]] inline int num_dof_per_element() const override { return 0; }

      /*!
      \brief Print this element
      */
      void print(std::ostream& os) const override;

      [[nodiscard]] Core::Elements::ElementType& element_type() const override
      {
        return StructuralSurfaceType::instance();
      }

      std::vector<std::shared_ptr<Core::Elements::Element>> lines() override;

      int num_line() const override;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the fluid3 element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,

      \return 0 if successful, negative otherwise
      */
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      //! Evaluate method for StructuralSurface-Elements
      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elematrix1,
          Core::LinAlg::SerialDenseMatrix& elematrix2, Core::LinAlg::SerialDenseVector& elevector1,
          Core::LinAlg::SerialDenseVector& elevector2,
          Core::LinAlg::SerialDenseVector& elevector3) override;

      //! Evaluate method for StructuralSurface-Elements
      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Elements::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elematrix1,
          Core::LinAlg::SerialDenseMatrix& elematrix2, Core::LinAlg::SerialDenseVector& elevector1,
          Core::LinAlg::SerialDenseVector& elevector2,
          Core::LinAlg::SerialDenseVector& elevector3) override;

      //! \brief Evaluate trace inequality and return the maximal eigenvalue
      //!
      //! \param[in] parent_disp    displacement values on corresponding parent element
      //! \return maximum eigenvalue
      double estimate_nitsche_trace_max_eigenvalue(const std::vector<double>& parent_disp);

      //! \brief Evaluate trace inequality and return the maximal eigenvalue
      //!
      //! \param[in] parent_disp    displacement values on corresponding parent element
      //! \param[in] parent_scalar  scalar values on corresponding parent element
      //! \return maximum eigenvalue
      double estimate_nitsche_trace_max_eigenvalue(
          const std::vector<double>& parent_disp, const std::vector<double>& parent_scalar);

      //! Evaluate trace inequality and return the maximal eigenvalue
      double estimate_nitsche_trace_max_eigenvalue_tsi(std::vector<double>& parent_disp);

      //! Return a pointer to the parent element of this boundary element
      [[nodiscard]] Core::Elements::Element* parent_element() const override
      {
        return parent_master_element();
      }

      //! Return local surface number
      [[nodiscard]] int l_surf_number() const { return face_master_number(); }

     public:  // changed to public rauch 05/2014
      //! action for surface evaluation
      enum ActionType
      {
        none,
        calc_struct_constrvol,
        calc_struct_volconstrstiff,
        calc_struct_constrarea,
        calc_struct_areaconstrstiff,
        calc_struct_monitarea,
        prestress_update,
        prestress_writerestart,
        prestress_readrestart,
        calc_struct_centerdisp,
        calc_struct_rotation,
        calc_undo_struct_rotation,
        calc_ref_nodal_normals,
        calc_cur_nodal_normals,
        calc_struct_area_poro
      };

     protected:
      // don't want = operator
      StructuralSurface& operator=(const StructuralSurface& old);

      // variables

      //! discretization type
      Core::FE::CellType distype_;

      //! numdofpernode
      int numdofpernode_;

      //! gaussian integration to be used
      Core::FE::GaussRule2D gaussrule_;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates normal and detA at gaussian point

      \param detA   (out) : are at gaussian point
      \param normal (out) : The normal at gaussian point, length is detA!
      \param x      (in)  : nodal coords in either material or spatial frame
      \param deriv  (in)  : derivatives of shape functions
      */
      void surface_integration(double& detA, std::vector<double>& normal,
          const Core::LinAlg::SerialDenseMatrix& x, const Core::LinAlg::SerialDenseMatrix& deriv);

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates the normal at gaussian point

      \param normal (out) : The normal at gaussian point, length is detA!
      \param x      (in)  : nodal coords in either material or spatial frame
      \param deriv  (in)  : derivatives of shape functions
      */
      void surface_integration(std::vector<double>& normal,
          const Core::LinAlg::SerialDenseMatrix& x, const Core::LinAlg::SerialDenseMatrix& deriv);

      /*!
      \brief Linearize a Neumann boundary condition analytically

      This method computes the analytical surface derivative if necessary (for follower loads).
      */
      void analytical_d_surface_integration(Core::LinAlg::SerialDenseMatrix& d_normal,
          const Core::LinAlg::SerialDenseMatrix& x, const Core::LinAlg::SerialDenseMatrix& deriv);

      /*!
      \brief Linearize a Neumann boundary condition with automatic differentiation

      This method computes the automatic (Sacado) surface derivative if necessary (for follower
      loads). This method is for development purposes only, and currently not used.
      */
      void automatic_d_surface_integration(Core::LinAlg::SerialDenseMatrix& d_normal,
          const Core::LinAlg::SerialDenseMatrix& x, const Core::LinAlg::SerialDenseMatrix& deriv);

      /*!
      \brief Create matrix with material configuration

      \param x  (out)  : nodal coords in material frame
      */
      inline void material_configuration(Core::LinAlg::SerialDenseMatrix& x) const
      {
        const int numnode = num_node();
        for (int i = 0; i < numnode; ++i)
        {
          x(i, 0) = nodes()[i]->x()[0];
          x(i, 1) = nodes()[i]->x()[1];
          x(i, 2) = nodes()[i]->x()[2];
        }
      }

      /*!
      \brief Create matrix with spatial configuration

      \param x     (out)  : nodal coords in spatial frame
      \param disp  (int)  : displacements
      */
      inline void spatial_configuration(
          Core::LinAlg::SerialDenseMatrix& x, const std::vector<double>& disp) const
      {
        const int numnode = num_node();
        for (int i = 0; i < numnode; ++i)
        {
          x(i, 0) = nodes()[i]->x()[0] + disp[i * 3 + 0];
          x(i, 1) = nodes()[i]->x()[1] + disp[i * 3 + 1];
          x(i, 2) = nodes()[i]->x()[2] + disp[i * 3 + 2];
        }
      }

      /*!
      \brief Create matrix with spatial configuration

      \param x     (out)  : nodal coords in spatial frame
      \param xrefe (out)  : nodal coords in material frame
      \param disp  (int)  : displacements
      */
      inline void spatial_configuration(Core::LinAlg::SerialDenseMatrix& x,
          Core::LinAlg::SerialDenseMatrix& xrefe, const std::vector<double>& disp) const
      {
        const int numnode = num_node();
        for (int i = 0; i < numnode; ++i)
        {
          x(i, 0) = xrefe(i, 0) + disp[i * 3 + 0];
          x(i, 1) = xrefe(i, 1) + disp[i * 3 + 1];
          x(i, 2) = xrefe(i, 2) + disp[i * 3 + 2];
        }
      }

      //! Submethod to compute the enclosed volume for volume constraint boundary condition
      double compute_constr_vols(
          const Core::LinAlg::SerialDenseMatrix& xc,  ///< current configuration
          const int numnode                           ///< num nodes
      );


      /*!
      \brief Submethod to compute interfacial area and its first and second
             derivatives w.r.t. the displacements (this is needed for surface energy problems)
      */
      void compute_area_deriv(const Core::LinAlg::SerialDenseMatrix& x,  ///< spatial configuration
          const int numnode,                                             ///< number of nodes
          const int ndof,                          ///< number of degrees of freedom
          double& A,                               ///< area
          Core::LinAlg::SerialDenseVector& Adiff,  ///< first derivative
          const std::shared_ptr<Core::LinAlg::SerialDenseMatrix>& Adiff2  ///< second derivative
      );

      //! Submethod to compute constraint volume and its first and second derivatives w.r.t. the
      //! displacements
      void compute_vol_deriv(const Core::LinAlg::SerialDenseMatrix& x,  ///< spatial configuration
          const int numnode,                                            ///< number of nodes
          const int ndof,                          ///< number of degrees of freedom
          double& V,                               ///< volume
          Core::LinAlg::SerialDenseVector& Vdiff,  ///< first derivative
          const std::shared_ptr<Core::LinAlg::SerialDenseMatrix>& Vdiff2,  ///< second derivative
          const int minind = 0,  ///< minimal index to compute enclosed volume with
          const int maxind = 2   ///< maximal index to compute enclosed volume with
      );

      /// Submethod to compute normal vectors at nodes (for reference configuration) which can be
      /// assembled (i.e., added) to obtain averaged nodal normal vectors
      void build_normals_at_nodes(Core::LinAlg::SerialDenseVector& nodenormals,
          const std::vector<double>& mydisp, bool refconfig);

      //! Submethod to compute surface porosity
      void calculate_surface_porosity(Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Elements::LocationArray& la);

     private:
      //! Private Constructor to be called from StructuralSurfaceType
      explicit StructuralSurface(int id, int owner);

      //! Set discretization Type of the Surface Element
      void set_distype();

      //! Set Gaussrule dependent on shape of the structural surface
      void set_gaussrule();


      /// @name private methods for calculation of trace inequality matrices
      //@{

      //! Templated version: parent and surface discretization type
      template <Core::FE::CellType dt_vol, Core::FE::CellType dt_surf>
      double estimate_nitsche_trace_max_eigenvalue(
          const std::vector<double>& parent_disp, const std::vector<double>& parent_scalar);

      //! the volume stiffness matrix
      //! unlike the "full" stiffness matrix we don't use the geometric term here
      template <Core::FE::CellType dt_vol>
      void trace_estimate_vol_matrix(
          const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xrefe,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
          const std::vector<double>& parent_scalar,
          Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol> * 3, Core::FE::num_nodes<dt_vol> * 3>&
              vol);

      //! the surface stiffness matrix
      template <Core::FE::CellType dt_vol, Core::FE::CellType dt_surf>
      void trace_estimate_surf_matrix(
          const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xrefe,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
          const std::vector<double>& parent_scalar,
          Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol> * 3, Core::FE::num_nodes<dt_vol> * 3>&
              surf);

      //! evaluate the kinematics
      template <Core::FE::CellType dt_vol>
      void strains(const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xrefe,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
          const Core::LinAlg::Matrix<3, 1>& xi, double& jac, Core::LinAlg::Matrix<3, 3>& defgrd,
          Core::LinAlg::Matrix<6, 1>& glstrain, Core::LinAlg::Matrix<3, 3>& rcg,
          Core::LinAlg::Matrix<6, Core::FE::num_nodes<dt_vol> * 3>& bop,
          Core::LinAlg::Matrix<3, Core::FE::num_nodes<dt_vol>>& N_XYZ);


      //! setup projector for removing the rigid body modes from the generalized eigenvalue problem
      template <Core::FE::CellType dt_vol>
      void subspace_projector(const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
          Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol> * Core::FE::dim<dt_vol>,
              Core::FE::num_nodes<dt_vol> * Core::FE::dim<dt_vol> -
                  Core::FE::dim<dt_vol>*(Core::FE::dim<dt_vol> + 1) / 2>& proj);


      //! Templated version: parent and surface discretization type
      template <Core::FE::CellType dt_vol, Core::FE::CellType dt_surf>
      double estimate_nitsche_trace_max_eigenvalue_tsi(std::vector<double>& parent_disp);

      //! the volume stiffness matrix
      //! unlike the "full" stiffness matrix we don't use the geometric term here
      template <Core::FE::CellType dt_vol>
      void trace_estimate_vol_matrix_tsi(
          const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xrefe,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
          Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, Core::FE::num_nodes<dt_vol>>& vol);

      //! the surface stiffness matrix
      template <Core::FE::CellType dt_vol, Core::FE::CellType dt_surf>
      void trace_estimate_surf_matrix_tsi(
          const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xrefe,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, 3>& xcurr,
          Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, Core::FE::num_nodes<dt_vol>>& surf);


      //! setup projector for removing the rigid body modes from the generalized eigenvalue problem
      template <Core::FE::CellType dt_vol>
      void subspace_projector_scalar(
          Core::LinAlg::Matrix<Core::FE::num_nodes<dt_vol>, Core::FE::num_nodes<dt_vol> - 1>& proj);
      //@}

    };  // class StructuralSurface : public Core::Elements::Element


  }  // namespace Elements
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
