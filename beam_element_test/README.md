# Single beam element test

The file `single_beam_element.4C.yaml` in this directory contains a simple test case for a single beam element.
The model is created using the `create_model.py` script.

To run the single element test, execute the following command in the terminal:
```bash
<path_to_4C> single_beam_element.4C.yaml xxx
```
here, `xxx` is the output name.
The debug output showing input and output values for the internal force and stiffness matrix calculation will be printed to the terminal.
The resulting deformed shape of the beam can be visualized using ParaView by opening the output file `xxx-structure-beams.pvd`.
