# RTMsim - A finite area method for filling simulations in liquid composite moulding

Resin Transfer Moulding (RTM) is a manufacturing process for producing thin-walled fiber reinforced polymer composites where dry fibers are placed inside a mould and resin is injected under pressure into the fibrous preform. During mould design, filling simulations can study different manufacturing concepts (i.e. placement of inlet ports and vents) to guarantee complete filling of the part and avoid air entrapment where flow fronts converge. 

In the past, numerous models have been implemented in different software packages to perform filling simulations for RTM. The used simulation packages can be divided into three groups: 
- General purpose CFD software packages, such as ANSYS Fluent or OpenFOAM
- Commercially available software packages which are tailored for the simulation of the RTM process, such as PAM-RTM, RTM-Worx or LIMS
- Easy-to-use simulation tools such as myRTM

All packages describe the flow on a macroscopic level. The first group models the flow through the porous cavity using volume-averaged Navier-Stokes equations. The second group makes use of some assumptions and solves in a first step a Laplace equation for the pressure inside the region which is already filled and in a second step calculates the flow velocity field to propagate the flow front. It has been shown that the first and second group render very similar results. myRTM from the third group is easy-to-use and can predict the filling pattern properly but neither predict the filling time correctly nor consider orthotropic preform permeability. Solving conservation laws for fluid flow as in the first group requires a volume mesh of the cavity and consequently the solution is more time-consuming. The second and third group can be solved on a shell mesh where the thickness of the cavity is a property of the cell (similar to porosity and permeability) and slip boundary conditions at the top and bottom walls of the cavity are assumed. 

Analyzing the existing software tools for RTM simulations gives the following functional requirements for a new software tool:
- The simulation model shall give correct results for filling pattern and filling time.
- The simulation tool takes only composite-manufacturing related inputs and the simulation shall be robust independent of numerics-related input.
- The simulation tool takes a shell model of the geometry as input and the location-dependent properties are assigned directly on the shell elements.
- New functionalities can be implemented by either adding equations of the same type or modifying existing equations.

Correct results provided, the requirements are summarized as robustness, easy-to-use and simple extendability. These requirements can be fulfilled by (i) an
appropriate choice of the flow physics model, (ii) robust mesh handling and (iii) a parameter-driven execution of the flow simulation runs as will be described in
the following paragraphs.

From a flow physics point of view the filling process is described by an incompressible resin which is injected under pressure into a cavity which initially
is either evacuated or filled with air. Since the simulation tool shall run robust and independent of numerics-related input like under-relaxation coefficients or time-steps, a one phase compressible model is chosen. Conservation laws for a compressible continuity equation, momentum equations in a local cell coordinate system, an adiabatic law as equation of state and a volume-of-fluid equation to track the resin flow front must be solved. 

GIVE EQUATINS MODEL HERE

All partial differential equations can be solved using the same method. Other conservation laws, e.g. temperature or degree-of-cure equations to include non-constant resin viscosity, can be added with a similar method. Modifications of the existing conservation laws, e.g. introducing a thickness-dependent permeability which is needed for modeling vacuum-assisted resin infusion (VARI), can be implemented by changing the existing discretized equations.

The new simulation tool does not include mesh generation. A mesh with the pre-defined regions must be generated with a meshing tool before starting the
filling simulation. Creating a shell mesh is easier than creating a volume mesh and the number of cells and the computational time can be reduced significantly
compared to a volume mesh which has to fulfill cell quality requirements (e.g. for skewness and orthogonality). Shell meshes are typically created on the part’s
mid-surface. Mid-surface models are often available in composite manufacturing since computational stress analysis for thin-walled parts is performed on the part’s mid-surface too. Since a shell mesh is used in the simulation tool, the conservation laws must be solved on a shell mesh using a generalization of the
so-called finite area method. 

The new simulation tool is executed with a well-defined list of parameters specified in an input text file or in the GUI.  

The Julia implementation of the presented model is validated and verified with the following test cases:
- A permeameter experiment with isotropic in-plane permeability to validate with literature results.
- A permeameter experiment with tilted orthotropic in-plane permeability to verify that the simulation gives the expected results.
- A curved permeameter to verify the implemented model for curved geometry by comparison with the flat permeameter results.
- An annulus filler-like part to verify the implemented model for a curved part with T-junctions by comparison with an ANSYS Fluent simulation and with a myRTM simulation.
- A permeameter experiment with two patches with different in-plane permeability and porosity levels.

In a follow-up paper simulations from a real-world RTM mould for a complex part will be analyzed and compared with experiments. This will be used to
validate the new simulation model for patches with different cavity thickness. 

In order to numerically solve the flow model the computational domain (time and space) is discretized. The temporal domain (i.e. the simulation time) is split into a finite number of time steps where values for the physical quantities on the spatial domain are calculated in a time-marching manner. The spatial domain (i.e. the flow volume) is split into wedge-type cells. Since the walls are assumed to be slip boundaries only one cell is used through the thickness and the cells are bounded by three control surfaces only. Therefore, the spatial domain is defined on the part’s mid-surface and the local thickness of the flow volume is a property of the cell. The mid-surface model can be curved and cells can have edges where more than two cells are connected to each other such as for handling T-junctions. The mid-surface model is divided into a finite number of triangular control areas. The control areas cover the spatial domain completely without overlapping. The cell-centered method where the nodes are placed at the centroid of the control volume is used. 

For discretizing the equations on the shell mesh of the part’s mid-surface it is assumed that the geometry was locally flat there. The neighbouring cells are
rotated about the common edges to lie in the plane of the considered cell. 

ADD SOMETHINK ABOUT THE FINITE AREA METHOD HERE.

The shell mesh is imported via a text file where nodes, elements and element sets are described in a format similar to the NASTRAN bulk data format. Every line contains ten fields of eight characters each. The first field contains the character name of the item. The input file for the permeameter reads:
```
SET 1 = 1,2,3,4,5,6,
7,8,9,10,11,12,
13,14,15,16
GRID 1 0.0 0.0 0.0
GRID 3 0.3 0.3 0.0
GRID 4 -0.3 0.3 0.0
[...]
GRID 331 -2.115-2-.113318 0.0
GRID 332 -.117148.1562872 0.0
GRID 333 .2271322.1925105 0.0
CTRIA3 1 0 15 9 16
CTRIA3 2 0 16 10 19
CTRIA3 3 0 19 11 17
[...]
CTRIA3 586 0 243 302 332
CTRIA3 587 0 262 333 259
CTRIA3 588 0 232 259 333
```
Nodes are described by the keyword GRID, followed by a grid number, followed by a blank and three fileds with the x, y and z coordinates of the node. The triangular cells are defined by the keyword CTRIA3, followed by a cell number, followed by a zero7, followed by the three node numbers which constitute the cell. Nodes and elements need not be sorted nor starting with one. Cell sets are defined by the keyword SET followed by ‘ N = ’ and the cell numbers separated by commas. Not more than 6 cell numbers per line. If another line is required for additional cell numbers, these follow after 8 blanks. Up to four sets can be defined. Mesh files of this type can be created with common meshing tools. The authors used Altair HyperWorks but also free software tools like SALOMEMECA, GMSH or NETGEN can be used.


DESCRIBE THE GUI AND THE INPUT TEXT FILE AND HOW TO START

SHOW RESULT PICS FOR THE VALIDATION CASES WHICH SHOW THAT THE CODE PERFORMS WELL

EXPLAIN HOW TO EXTEND THE CODE FOR i_model=2,3,.. and i_method=2,3,... for numerical differentiation and flux functions
