// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_FLD_FBI_MOVINGBOUNDARY_HPP
#define FOUR_C_ADAPTER_FLD_FBI_MOVINGBOUNDARY_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_moving_boundary.hpp"
#include "4C_fluid_meshtying.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations

namespace Core::LinAlg
{
  class SparseOperator;
}

namespace Adapter
{
  class Fluid;

  /// fluid with moving interfaces
  class FBIFluidMB : public FluidMovingBoundary
  {
   public:
    /// constructor
    explicit FBIFluidMB(const Teuchos::ParameterList& prbdyn, std::string condname);

    /*========================================================================*/
    //! @name Misc
    /*========================================================================*/

    /// fluid field
    const std::shared_ptr<Adapter::Fluid>& fluid_field() override { return fluidadapter_; }

    /// return the boundary discretization that matches the structure discretization
    std::shared_ptr<Core::FE::Discretization> discretization() override;

    /// communication object at the interface
    std::shared_ptr<FLD::Utils::MapExtractor> const& interface() const override;

    //@}

    /*========================================================================*/
    //! @name Time step helpers
    /*========================================================================*/

    /// start new time step
    void prepare_time_step() override;

    /// update at time step end
    void update() override;

    /// output results
    void output() override;

    /// read restart information for given time step
    double read_restart(int step) override;

    /*========================================================================*/
    //! @name Solver calls
    /*========================================================================*/

    /// nonlinear solve
    void nonlinear_solve(std::shared_ptr<Core::LinAlg::Vector<double>> idisp,
        std::shared_ptr<Core::LinAlg::Vector<double>> ivel) override;

    /// relaxation solve
    std::shared_ptr<Core::LinAlg::Vector<double>> relaxation_solve(
        std::shared_ptr<Core::LinAlg::Vector<double>> idisp, double dt) override;
    //@}

    /*========================================================================*/
    //! @name Extract interface forces
    /*========================================================================*/

    /// After the fluid solve we need the forces at the FSI interface.
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_forces() override;
    //@}

    /*========================================================================*/
    //! @name extract helpers
    /*========================================================================*/

    /// extract the interface velocity at time t^(n+1)
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_velnp() override;

    /// extract the interface velocity at time t^n
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_veln() override;
    //@}

    /*========================================================================*/
    //! @name Number of Newton iterations
    /*========================================================================*/

    /// get the maximum number of iterations from the fluid field
    int itemax() const override;

    /// set the maximum number of iterations for the fluid field
    void set_itemax(int itemax) override;

    //@}

    /*========================================================================*/
    //! @name others
    /*========================================================================*/

    /// integrate the interface shape functions
    std::shared_ptr<Core::LinAlg::Vector<double>> integrate_interface_shape() override;

    /// create the testing of fields
    std::shared_ptr<Core::Utils::ResultTest> create_field_test() override;

    /// Get velocity at timestep n+1
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> velnp();

    virtual std::shared_ptr<const FLD::Meshtying> get_meshtying();

    /** \brief Pass in additional contributions from coupling terms for the system matrix
     *
     * To enforce weak dirichlet conditions as they arise from meshtying for example, such
     * contributions can be handed to the fluid, which will store the pointer on the coupling
     * contributions to assemble them into the system matrix in each Newton iteration.
     *
     * \param[in] matrix (size fluid_dof x fluid_dof) matrix containing weak dirichlet entries that
     * need to be assembled into the overall fluid system matrix
     */
    virtual void set_coupling_contributions(
        std::shared_ptr<const Core::LinAlg::SparseOperator> matrix);

    /**
     * \brief Pass additional contributions to the fluid residual to the fluid class
     *
     * \param[in] iforce contribution to the fluid residual
     *
     * \param[in] ivel unused in this implementation     *
     *
     */
    void apply_interface_values(std::shared_ptr<Core::LinAlg::Vector<double>> iforce,
        std::shared_ptr<Core::LinAlg::Vector<double>> ivel = nullptr) override;

    /**
     * \brief Resets the external forces acting on the fluid to zero
     */
    virtual void reset_external_forces();

   private:
    /// fluid base algorithm object
    std::shared_ptr<Adapter::Fluid> fluidadapter_;
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
