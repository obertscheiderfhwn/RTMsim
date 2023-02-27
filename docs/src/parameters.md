# Input parameters
RTMsim is executed with a well-defined list of parameters specified in an input text file or in the GUI. The mesh file and all parameters must be specified in SI units.

## GUI 
The following figure shows the GUI with explanation of the parameters.
```@raw html
<img src=../assets/figures/rtmsim_help.png>
```

The buttons in the line on the RHS are used to start the simulation with the parameters from the selected input file. For example, click on `Select input file` and select the input.txt in the inputfiles folder from the package installation and click on `Run with input file`. After the simulation is completed, click on `Plot overview`.

Parameters (fluid properties, patch types and patch properties of cell sets specified in the mesh file) can also be specified in the GUI and a simulation is then started by clicking on `Start simulation`. Also other functionalities are available. The buttons in the first line on the LHS are used for mesh inspection, i.e. `Select mesh file`, `Plot mesh` with bounding box and `Plot sets` for inspecting the defined sets in the mesh file. The buttons `Start simulation` and `Continue simulation` in the second line on the LHS are used for starting and continuing a filling simulation. Every time the Start or Continue simulation button is pressed, a filling simulation is started. The simulated flow time `tmax` is specified in the first field in the second line. Every simulation calculates the flow front propagation during the next `tmax` seconds. If started with the `Start simulation` button, the cavity is empty initially. If started with the `Continue simulation` button, the results from the previous simulation are taken as initial condition. With the buttons `Start interactive` and `Continue interactive` in the third line one can start and continue simulations where manually selected inlet ports are used in addition to sets defined below. The radius of the inlet ports is specified in the first field in the third line and the locations are selected with the mouse after clicking on `Select inlet port`. The buttons `Plot results`, `Plot overview`and `Plot filling` in the forth line are used for creating contour plots, i.e. show filling and pressure distribution of a specified output file (path to the results file in the second cell and can be changed by clicking on `Select results file`; final results are saved in `results.jld2`), plot filling at four equidistant time instances and filling at different time instances which are selected with a slider bar. 

## Input file
The complete set of input parameters can be accessed in the input file. The following paragraph shows an example for such an input file:
```
1    #i_model 
meshfiles/mesh_permeameter1_foursets.bdf    #meshfilename 
200    #tmax 
1.01325e5 1.225 1.4 0.06    #p_ref rho_ref gamma mu_resin_val 
1.35e5 1.0e5    #p_a_val p_init_val 
3e-3 0.7 3e-10 1 1 0 0    #t_val porosity_val K_val alpha_val refdir1_val refdir2_val refdir3_val 
3e-3 0.7 3e-10 1 1 0 0    #t1_val porosity1_val K1_val alpha1_val refdir11_val refdir21_val refdir31_val 
3e-3 0.7 3e-10 1 1 0 0    #t2_val porosity2_val K2_val alpha2_val refdir12_val refdir22_val refdir32_val
3e-3 0.7 3e-10 1 1 0 0    #t3_val porosity3_val K3_val alpha3_val refdir13_val refdir23_val refdir33_val
3e-3 0.7 3e-10 1 1 0 0    #t4_val porosity4_val K4_val alpha4_val refdir14_val refdir24_val refdir34_val 
1 0 0 0    #patchtype1val patchtype2val patchtype3val patchtype4val 
0 results.jld2    #i_restart restartfilename
0 0.01    #i_interactive r_p
16    #n_pics
```

Meaning of the variables:
- `i_model`: Identifier for physical model (Default value is 1)
- `meshfilename`: Mesh filename.
- `tmax`: Maximum simulation time.
- `p_ref rho_ref gamma mu_resin_val`: Parameters for the equation of state and dynamic viscosity of resin used in the Darcy term.
- `p_a_val p_init_val `: Absolut pressure value for injection port and for initial cavity pressure.
- `t_val porosity_val K_val alpha_val refdir1_val refdir2_val refdir3_val`: Properties of the cells in the main preform: The vector `(refdir1_val,refdir2_val,refdir3_val)` is projected onto the cell in order to define the first principal cell direction. The second principal cell direction is perpendicular to the first one in the plane spanned by the cell nodes. The principal cell directions are used as the principal permeabilty directions. The cell properties are defined by the thickness `t_val`, the porosity `porosity_val`, the permeability `K_val` in the first principal cell direction, the permeablity `alpha_val` in the second principal direction.
- `t1_val porosity1_val K1_val alpha1_val refdir11_val refdir21_val refdir31_val` etc.: Properties for up to four additional cell regions if preform. 
- `patchtype1val patchtype2val patchtype3val patchtype4val`: These regions are used to specify the location of the pressure boundary conditions and to specify regions with different permeability, porosity and thickness properties (e.g. for different part thickness and layup or for race tracking which are regions with very high permeability typically at the boundary of the preforms). Vents need not be specified. Parameters `patchtype1val` define the patch type. Numerical values 0, 1, 2 and 3 are allowed with the following interpretation:
    - 0 .. the patch is ignored
    - 1 .. the patch represents an inlet gate, where the specified injection pressure level applies
    - 2 .. the patch specifies a preform region
    - 3 .. the patch represents a vent, where the specified initial pressure level applies
- `i_restart restartfilename`: Start with new simulation if `0` or continue previous simulation if `1`.
- `i_interactive r_p`: Select the inlet ports graphically if i_interactive equal to `1`.
- `n_pics`: Number of intermediate output files. Supposed to be a multiple of `4`.
Entries are separated by one blank.


## Mesh file
The new simulation tool does not include mesh generation. A 3-node triangular shell mesh with the pre-defined regions must be generated with a meshing tool before starting the filling simulation. The authors used Altair HyperWorks but also free software tools such as SALOME-MECA, GMSH or NETGEN can be used.

The shell mesh is created on the part's mid-surface. Mid-surface models are often available in composite manufacturing since computational stress analysis for thin-walled parts is performed on the part’s mid-surface too. 

The prepared shell mesh is imported via a text file where nodes, elements and element sets are described in a format similar to the NASTRAN bulk data format. Every line (except for the set definitions) contains ten fields of eight characters each. The first field contains the name of the item. The following paragraph shows an example for such a mesh file:
```
SET 1 = 1,2,3,4,5,6,
        7,8,9,10,11,12,
        13,14,15,16
GRID           1             0.0     0.0     0.0
GRID           3             0.3     0.3     0.0
GRID           4            -0.3     0.3     0.0
[...]
GRID         331        -2.115-2-.113318     0.0
GRID         332        -.117148.1562872     0.0
GRID         333        .2271322.1925105     0.0
CTRIA3         1       0      15       9      16
CTRIA3         2       0      16      10      19
CTRIA3         3       0      19      11      17
[...]
CTRIA3       586       0     243     302     332
CTRIA3       587       0     262     333     259
CTRIA3       588       0     232     259     333
```
Nodes are described by the keyword `GRID`, followed by a grid number, followed by a blank and three fields with the x, y and z coordinates of the node. The triangular cells are defined by the keyword `CTRIA3`, followed by a cell number, followed by a zero, followed by the three node numbers which constitute the cell. Nodes and cells need not be sorted nor starting with one. Cell sets are defined by the keyword `SET` followed by ` N = ` and the cell numbers separated by commas. Not more than 6 cell numbers per line. Then continue with the cell numbers in the line below after 8 blanks at the beginning of the line. Up to four sets can be defined. 








```@meta
EditURL = "https://github.com/obertscheiderfhwn/RTMsim/blob/main/docs/src/parameters.md"
```