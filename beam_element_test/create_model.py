# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Embedded elbow structure."""

import numpy as np
from beamme.core.boundary_condition import BoundaryCondition
from beamme.core.conf import bme
from beamme.core.function import Function
from beamme.core.mesh import Mesh
from beamme.core.rotation import Rotation
from beamme.four_c.element_beam import Beam3rLine2Line2
from beamme.four_c.header_functions import set_runtime_output
from beamme.four_c.input_file import InputFile
from beamme.four_c.material import MaterialReissner
from beamme.mesh_creation_functions.beam_arc import create_beam_mesh_arc_segment_2d


def create_model(n_el: int = 1):
    """Create a single element model."""

    input_file = InputFile()
    set_runtime_output(input_file, btsvmt_output=False, btss_output=False)
    input_file.add(
        {
            "PROBLEM TYPE": {"PROBLEMTYPE": "Structure"},
            "IO": {
                "OUTPUT_BIN": False,
                "STRUCT_DISP": False,
                "VERBOSITY": "Standard",
            },
            "STRUCTURAL DYNAMIC": {
                "LINEAR_SOLVER": 1,
                "INT_STRATEGY": "Standard",
                "DYNAMICTYPE": "Statics",
                "RESULTSEVERY": 0,
                "NLNSOL": "fullnewton",
                "PREDICT": "ConstDis",
                "TIMESTEP": 1.0,
                "NUMSTEP": 1,
                "MAXTIME": 1.0,
                "MAXITER": 10,
                "TOLRES": 1.0e-12,
                "TOLDISP": 1.0e-12,
            },
            "SOLVER 1": {"NAME": "Structure_Solver", "SOLVER": "Superlu"},
        }
    )

    mesh = Mesh()
    ft = Function([{"SYMBOLIC_FUNCTION_OF_TIME": "t"}])
    mesh.add(ft)

    mat = MaterialReissner(youngs_modulus=100.0, radius=0.1)
    beam_set = create_beam_mesh_arc_segment_2d(
        mesh,
        Beam3rLine2Line2,
        mat,
        center=[0.0, 1.0, 0.0],
        radius=1.0,
        phi_start=-np.pi / 2,
        phi_end=0.0,
        n_el=n_el,
    )
    mesh.add(
        BoundaryCondition(
            beam_set["start"],
            {"NUMDOF": 6, "ONOFF": [1] * 6, "VAL": [0] * 6, "FUNCT": [None] * 6},
            bc_type=bme.bc.dirichlet,
        )
    )
    mesh.add(
        BoundaryCondition(
            beam_set["end"],
            {
                "NUMDOF": 6,
                "ONOFF": [0, 0, 1, 0, 0, 0],
                "VAL": [0, 0, 0.001, 0, 0, 0],
                "FUNCT": [None, None, ft, None, None, None],
            },
            bc_type=bme.bc.neumann,
        )
    )

    input_file.add(mesh)
    input_file.dump(
        "single_beam_element.4C.yaml",
        mesh_format="yaml",
        add_footer_application_script=False,
    )


if __name__ == "__main__":
    # Execution part of script.
    create_model()
