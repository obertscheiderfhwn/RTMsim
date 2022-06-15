---
title: 'RTMsim - A Julia module for filling simulations in Resin Transfer Moulding with the Finite Area Method'
tags:
  - Computational Fluid Dynamics (CFD)
  - Shell mesh
  - Resin Transfer moulding (RTM)
  - Liquid composite moulding (LCM)
  - Filling simulation
  - Julia
authors:
  - name: Christof Obertscheider^[Corresponding author]
    affiliation: 1
  - name: Ewald Fauster
    affiliation: 2
affiliations:
 - name: Aerospace Engineering Department, University of Applied Sciences Wiener Neustadt, Johannes-Gutenberg-Straße 3, 2700, Wiener Neustadt, Austria
   index: 1
 - name: Processing of Composites Group, Department Polymer Engineering and Science, Montanuniversität Leoben, Otto Glöckl-Straße 2, 8700 Leoben, Austria 
   index: 2
date: 7 April 2022
bibliography: paper.bib
---

# Summary
Resin Transfer Moulding (RTM) is a manufacturing process for producing thin-walled fiber reinforced polymer composites where dry fibers are placed inside a mould and resin is injected under pressure into the fibrous preform. RTMsim is a robust, easy-to-use and simple-to-extend software tool, written in Julia, for RTM filling simulations. A shell mesh, injection pressure, resin viscosity and the parameters describing the preform are required input. The software was validated with results from literature as well as experiments and compared with results from well-established RTM filling simulation tools.

# Statement of need
Resin Transfer Moulding (RTM) is a manufacturing process for producing thin-walled fiber reinforced polymer composites where dry fibers are placed inside a mould and resin is injected under pressure into the fibrous preform. The thickness of the part is much smaller than the overall dimensions as depicted in \autoref{fig:rtm_pic1}. During mould design, filling simulations can study different manufacturing concepts (i.e. placement of inlet ports and vents) to guarantee complete filling of the part and avoid air entrapment where flow fronts converge. 

![Schematic of a thin-shell RTM geometry for manufacturing a flat plate with linear (a) and radial (b) flow.\label{fig:rtm_pic1}](figures/rtm_pic1.png)

In the past, numerous models have been implemented in different software packages to perform filling simulations for RTM. The used simulation packages can be divided into three groups: 

- General purpose CFD software packages, such as ANSYS Fluent or OpenFOAM

- Commercially available software packages which are tailored for the simulation of the RTM process, such as PAM-RTM, RTM-Worx or LIMS

- Easy-to-use simulation tools such as myRTM

All packages describe the flow on a macroscopic level. The first group models the flow through the porous cavity using volume-averaged Navier-Stokes equations [@groessing]. The second group makes use of some assumptions and solves in a first step a Laplace equation for the pressure inside the region which is already filled and in a second step calculates the flow velocity field to propagate the flow front, see e.g. [@rtmworx] or [@lims]. In [@groessing] it has been shown that the first and second group render very similar results. myRTM from the third group is easy-to-use and can predict the filling pattern properly but neither predict the filling time correctly nor consider orthotropic preform permeability, see [@myrtm]. Solving conservation laws for fluid flow as in the first group requires a volume mesh of the cavity and consequently the solution is more time-consuming. The second and third group can be solved on a shell mesh where the thickness of the cavity is a property of the cell (similar to porosity and permeability) and slip boundary conditions at the top and bottom walls of the cavity are assumed. 

Based on the analysis of the existing software tools for RTM filling simulations the following functional requirements for a new software tool were derived: 

- The simulation model shall give correct results for filling pattern and filling time. 

- The simulation tool takes only composite-manufacturing related inputs and the simulation shall be robust independent of numerics-related input. 

- The simulation tool takes a shell model of the geometry as input and the location-dependent properties are assigned directly on the shell elements. 

- New functionalities can be implemented by either adding equations of the same type or modifying existing equations.

RTMsim ia a new software tool for RTM filling simulations which fulfills these requirements: Several test cases were used for successfully validating the implemented model. The simulation run robustly and independent of numerics related input such as under-relaxation coefficients or time step. The porous cavity is fully described by a mesh file with triangular cells on the part’s mid-surface and cell set definitions (for specifying the location of the pressure injection ports and vents and regions with different preforms by assigning different thickness, porosity and permeability values). Additional equations (e.g. for modeling the degree-of-cure) can be added with equations of the same type and modifications of existing equations (e.g. for variable cavity thickness as needed for vacuum assisted resin infusion simulations) is possible with cell properties which vary over simulation time. 

# Modeling

From a flow physics point of view the filling process is described by an incompressible resin which is injected under pressure into a cavity which initially is either evacuated or filled with air. Since the simulation tool shall run robust and independent of numerics-related input, a one phase compressible model is chosen. The governing equations shown in [@groessing] are modified. Conservation laws for a compressible continuity equation, momentum equations in a local cell coordinate system, an adiabatic law as equation of state and a volume-of-fluid equation to track the resin flow front must be solved. The local cell coordinate system is used to describe the flow in the coordinate system of the orthotropic in-plane permeability (permeability $k_1$ in first principal direction and permeability $k_2=\alpha \cdot k_1$ in second principal direction). All partial differential equations can be solved using the same method. 

Since the filling is a result of the pressure difference between injection pressure and initial cavity pressure, gauge pressures are used for the simulation. The initial pressure in the cavity is set equal to zero and the injection pressure is set equal to the pressure difference. This normalization to gauge pressures results in a filling which is independent of the pressure level and only depends on the pressure difference. 

In order to numerically solve the flow model the computational domain (time and space) is discretized. The temporal domain (i.e. the simulation time) is split into a finite number of time steps where values for the physical quantities on the spatial domain are calculated in a time-marching manner. The spatial domain (i.e. the flow volume) is split into wedge-type cells. Since the walls are assumed to be slip boundaries only one cell is used through the thickness and the cells are bounded by three control surfaces only. Therefore, the spatial domain is defined on the part’s mid-surface and the local thickness of the flow volume is a property of the cell. The mid-surface model can be curved and cells can have edges where more than two cells are connected to each other such as for handling T-junctions. The mid-surface model is divided into a finite number of triangular control areas. The control areas cover the spatial domain completely without overlapping. The cell-centered method where the nodes are placed at the centroid of the control volume is used. 

Since a shell mesh is used in the simulation tool, the conservation laws must be solved on a shell mesh using a generalization of the so-called finite area method [@tukovic]. The governing equations are discretized for every cell of the shell mesh. For the evaluation of the numerical flux functions and the pressure gradient, the neighbouring triangular cells are rotated about the common edge to lie in the plane of the considered cell. \autoref{fig:mesh_pic1} shows the shell mesh which is used for the test cases 1 and 2.

![Volume mesh (a) and shell mesh (b) and domain with highlighted inlet region (c) for the simulation of a radial flow experiment.\label{fig:mesh_pic1}](figures/mesh_pic1.png)



# Validation and verification

Four different test cases are available, successfully validating the Julia implementation of the RTM filling model:

1. Validation of the software for radial flow with isotropic in-plane peremablity: The simulated flow front position after 200 s is compared with the calculated flow front postion from literature. 

2. Verification of the software for radial flow with tilted orthotropic in-plane permeablity: The simulated tilted elliptical flow front is analysed and the calculated orthotropic permeablity is compared with the input (permeability $k_1$ and $k_2$, angle $\theta$ of the first principal direction to the horizontal). 

3. Comparison of the simulated flow front position for a complex annulus filler-like part with the simulated flow front position with Ansys Fluent and comparison of the simulated filling pattern with results from a myRTM simulation. 

4. Validation with experimental data from a radial permeameter experiment with two patches with different in-plane permeability and porosity levels. 

Most of the test cases are so-called numerical permeameter experiments, see e.g. [@Fauster.2019]. The third test case is a filling study of a typical composite part which is manufactured using RTM, see [@barandun].  

\autoref{fig:pic1} and \autoref{fig:pic2} show the binary filling fraction for cases 1, 2 and 3. For every cell $i$ the fraction function is either $0$ or $1$, depending on the filling fraction being $\leq$ or $>$ than a threshold value. The threshold for the volume-of-fluid method is determined from the permeameter simulation with isotropic permeabilities (case 1). The circular and elliptical flow fronts for cases 1 and 2 are calculated from the simulated flow front with an image processing algorithm [@Fauster.2019]. For case 1 also a comparison with results from an analytical formula [@advani] is shown. The results of validation case 4 and the new physical model for RTM filling will be discussed in a follow-up paper. 

![Simulated flow front for case 1 after 200 s for a coarse and a fine mesh and flow front positions at different time instances calculated with an analytical formula and calculated from the simulated flow front with an image processing algorithm for a coarse and a fine mesh. Simulation input is porosity $0.70$, permeability $3.0 \cdot 10^{−10}$ m$^2$ in both principal directions which are aligned with the edges of the domain. The values for the fine mesh agree well with the results from the analytical formula. The values for the coarse mesh show an error of $\leq 15$% which decreases significantly in the course of time but the shape is no smooth circle.\label{fig:pic1}](figures/validation_pic1.png)

![Simulated flow front for case 2 after 200 s for a coarse and a fine mesh. Simulation input is porosity $0.70$, permeability $3.0 \cdot 10^{−10}$ m$^2$ in first principal direction, $30^\circ$ to the horizontal and permeability $1.5 \cdot 10^{−10}$ m$^2$ in second principal direction. An algorithm for determining the flow front position from an optical permeameter was adapted to calculate the permeablity: The results are $2.91 \cdot 10^{−10}$ m$^2$, $1.41 \cdot 10^{−10}$ m$^2$ and angle $30^\circ$ for the fine mesh and $2.43 \cdot 10^{−10}$ m$^2$, $1.17 \cdot 10^{−10}$ m$^2$ and $29^\circ$ for the coarse mesh. For the fine mesh (2198 cells for a domain with 600 × 600 mm), this reverse engineering shows very good agreement of the calculated permeability values with the values used as simulation input. For the coarse mesh (588 cells) the agreement is still acceptable.\label{fig:pic2}](figures/validation_pic2.png)

![Results for test case 3. The results for the new simulation tool with the fine mesh are in good agreement with the ANSYS Fluent results for both filling pattern and filling time. myRTM shows similar filling patterns for same filling degrees as the new simulation tool with the coarse mesh. This comparison shows that a coarse mesh can be used for filling simulations where only the flow front propagation is investigated. If the filling time is important a mesh refinement study must be performed. \label{fig:pic3}](figures/validation_pic3.png)


# Outlook
A customized version of the implemented PDE solver is used for a bachelor thesis at ISSE Augsburg. Reinforcement learning (RL) algorithms are used to control the RTM process. The aim is to prevent the formation of dry spots by controlling the pressure values of the individual injection valves. RL is a form of machine learning and learns by interacting with a suitable environment. Simulations are often used as a training environment.  A customized lightweight implementation of RTMsim was used as training environment. 

# In the course of the bachelor thesis, the use of reinforcement learning (RL) algorithms to control the RTM process is to be tested. The aim is to prevent the formation of dry spots by controlling the pressure values of the individual injection valves. The method used for this purpose, RL, is a form of machine learning and learns by interacting with a suitable environment. Simulations are often used as a training environment, as this is usually cheaper and, for example, much more time efficient due to the use of parallelism. This is where Christof Obertscheider's solver comes into play, the source code of which he kindly made available to me. By using my own lightweight implementation, I was able to integrate it into my RL project much better than would have been possible with one of the well-known CFD solutions. In addition, the RL algorithm can now interact with many instances of the simulation running in parallel at the same time, which makes training possible in a reasonable time frame in the first place.

# Im Zuge der Bachelorarbeit soll die Verwendung von Reinforcement-Learning(RL) Algorithmen zur Steuerung des RTM-Prozesses getestet werden. Es soll die Bildung von Trockenstellen verhindert werden, indem die Druckwerte der einzelnen Einspritzventile gesteuert werden. Die dazu verwendete Methode, RL, ist eine Form des Maschinellen Lernens und lernt durch Interaktion mit einer geeigneten Umgebung. Häufig werden Simulationen als Trainingsumgebung genutzt, da dies meist billiger ist und beispielsweise durch den Einsatz von Parallelität deutlich zeiteffizienter ist. An dieser Stelle kommt der Solver von Christof Obertscheider ins Spiel, dessen Quellcode er mir freundlicherweise zur Verfügung gestellt hat. Durch eine eigene, leichtgewichtige Implementierung konnte ich diesen deutlich besser in mein RL-Projekt einbinden, als es mit einer der bekannten CFD-Lösungen möglich gewesen wäre. Außerdem kann der RL-Algorithmus nun mit vielen parallel laufenden Instanzen der Simulation gleichzeitig interagieren, wodurch das Training überhaupt erst in einem sinnvollen Zeitrahmen möglich wird.

The source code is prepared for the following extensions:
- Import of mesh file in different format.
- Additional functionalities, e.g. adding temperature and degree-of-cure equations with variable resin viscosity ar for VARI with variable porosity and permeability.
- Different methods for numerical differentiation and for the calculation of the numerical flux functions.




# References

