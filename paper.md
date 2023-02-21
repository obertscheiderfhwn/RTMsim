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
RTMsim is a robust, easy-to-use and simple-to-extend software tool, implemented in Julia [@doi:10.1137/141000671] for resin transfer moulding (RTM) filling simulations. A shell mesh, injection pressure, resin viscosity and the parameters describing the preform are required input. The software was validated with results from literature [@rudd1997liquid] and compared with results from well-established RTM filling simulation tools [@myrtm],[@fluent].

# Statement of need
Resin Transfer Moulding (RTM) is a manufacturing process for producing thin-walled fiber reinforced polymer composites where dry fibers are placed inside a mould and resin is injected under pressure into the fibrous preform. The thickness of the part is much smaller than the overall dimensions $$H \ll L$$ as depicted in \autoref{fig:rtm_pic1}. During mould design, filling simulations can study different manufacturing concepts (i.e. placement of inlet ports and vents) to guarantee complete filling of the part and avoid air entrapment where flow fronts converge. 

![Schematic of a thin-shell RTM geometry for manufacturing a flat plate with linear (a) and radial (b) flow.\label{fig:rtm_pic1}](figures/rtm_pic1.png)

In the past, numerous models have been implemented in different software packages to perform filling simulations for RTM. Those can be divided into three groups: 

- General purpose CFD software packages, such as ANSYS Fluent or OpenFOAM

- Commercially available software packages which are tailored for the simulation of the RTM process, such as PAM-RTM, RTM-Worx or LIMS

- Easy-to-use simulation tools such as myRTM

All packages describe the flow on a macroscopic level. The first group models the flow through the porous cavity using volume-averaged Navier-Stokes equations [@groessing]. The second group makes use of some assumptions on flow physics and boundary conditions and at first solves a Laplace equation for the pressure inside the filling fluid (resin) and in a second step calculates the flow velocity field to propagate the flow front, see e.g. [@rtmworx] or [@lims]. In [@groessing] it has been shown that the first and second group render very similar results. myRTM from the third group is easy-to-use and can predict the filling pattern properly but neither predict the filling time correctly nor consider orthotropic preform permeability, see [@myrtm]. Solving conservation laws for fluid flow as in the first group requires a volume mesh of the cavity and consequently the solution is often time-consuming. The second and third group can be solved on a shell mesh where the thickness of the cavity is a property of the cell (similar to porosity and permeability) and slip boundary conditions at the top and bottom walls of the cavity are assumed. 

Based on the analysis of the existing software tools for RTM filling simulations the following functional requirements for a new software tool were derived: 

- The simulation model shall give correct results for filling pattern and filling time. 

- The simulation tool takes only composite-manufacturing related inputs and the simulation shall be robust and independent of under-relaxation coefficients or time step. 

- The simulation tool takes a shell model of the geometry as input and the location-dependent properties are assigned directly on the shell elements. 

- New modules can be implemented by either adding equations of the same type or modifying existing equations.

RTMsim ia a new software tool for RTM filling simulations which fulfills these requirements: Several test cases were used for successfully validating the implemented model. The simulations run robustly and independent of numerics related input such as under-relaxation coefficients or time step. The porous cavity is fully described by a mesh file with triangular cells on the part’s mid-surface and cell set definitions (for specifying the location of pressure injection ports and vents and for specifying regions with different properties than the main preform by assigning different thickness, porosity and permeability values). Additional equations (e.g. for modeling the degree-of-cure) can be added with equations of the same type and modifications of existing equations (e.g. for variable cavity thickness as needed for vacuum assisted resin infusion simulations [@strong2008fundamentals]) are possible with cell properties which vary over simulation time. 

# Modeling

From a flow physics point of view the filling process is described by an incompressible resin which is injected under pressure into a cavity which initially is either evacuated or filled with air. Since the simulation tool shall run robust and independent of numerics-related input, a one phase compressible model is chosen. The governing equations shown in [@groessing] are modified. Conservation laws for a compressible continuity equation, momentum equations in a local cell coordinate system, an adiabatic law as equation of state and a volume-of-fluid equation to track the resin flow front must be solved. The local cell coordinate system is used to describe the flow in the coordinate system of the orthotropic in-plane permeability (permeability $k_1$ in first principal direction and permeability $k_2=\alpha \cdot k_1$ in second principal direction). 

The physical quantities in the governing equations are functions of space $\mathbf{x}$ and time $t$ where $\mathbf{x} \in \Omega \subset \mathbb{R}^3$ and $\Omega$ is the fluid body also called the fluid cavity. The governing equations are
$$ \frac{\partial  \rho}{\partial t} + \nabla  \cdot  \left( \rho  \mathbf{u} \right) = 0 $$
$$ \frac{\partial \rho \mathbf{u}}{\partial t} + \nabla \cdot \left( \rho \mathbf{u} \mathbf{u}  \right) = - \nabla p -\mu \mathsf{K}^{-1} \mathbf{u} $$
$$ \frac{\partial \varepsilon c}{\partial t} + \mathbf{u} \cdot \left( \nabla c \right) = 0 $$
$$ p = \kappa \rho^\gamma $$
where

- $\rho$ is the mass density of the fluid, 

- $\mathbf{u}$ is the superficial velocity which is related to the physical velocity via the porosity $\varepsilon$ given by the stationary fibrous media, 

- $\mathbf{u}\mathbf{u}$ is the dyadic product, 

- $p$ is the pressure, 

- $-\mu \mathsf{K}^{-1} \mathbf{u}$ is the pressure loss from flow through a porous medium according to Darcy's law with dynamic viscosity $\mu$ and permeability tensor $\mathsf{K}$, 

- $c$ is the continuous fraction function with $0 \leq c \leq 1$, $c=0$ if empty and $c=1$ if completely filled, 

- $\gamma$ is the adiabatic index, and 

- $\kappa=p_{\rm ref}/\rho_{\rm ref}^\gamma$ is a constant determined from pressure and mass density reference values.

Since the filling is a result of the pressure difference between injection pressure and initial cavity pressure, gauge pressures are used for the simulation. The initial pressure in the cavity is set equal to zero and the injection pressure is set equal to the pressure difference. This normalization results in a filling which is independent of the pressure level and only depends on the pressure difference. 

In order to numerically solve the flow model the computational domain (time and space) is discretized. The temporal domain (i.e. the simulation time) is split into a finite number of time steps where values for the physical quantities on the spatial domain are calculated in a time-marching manner. The spatial domain (i.e. the flow volume) is split into wedge-type cells. Since the walls are assumed to be slip boundaries only one cell is used through the thickness and the cells are bounded by three control surfaces only. Therefore, the spatial domain is defined on the part’s mid-surface and the local thickness of the flow volume is a property of the cell. The mid-surface model can be curved and cells can have edges where more than two cells are connected to each other such as for handling T-junctions. The mid-surface model is divided into a finite number of triangular control areas. The control areas cover the spatial domain completely without overlapping. The cell-centered method where the nodes are placed at the centroid of the control volume is used. 

Since a shell mesh is used in the simulation tool, the conservation laws must be solved on a shell mesh using a generalization of the so-called finite area method [@tukovic]. The governing equations are discretized for every cell of the shell mesh.

the velocity vector has only two non-zero components u and v described in a local coordinate system. The local coordinate system is the result of a projection of a reference vector. The reference vector is chosen such that the local x- and y-directions correspond with the first and second principal permeability directions. For the evaluation of the numerical flux functions and the pressure gradient, the neighbouring triangular cells are rotated about the common edge to lie in the plane of the considered cell and a coordinate transformation for the velocity components must be performed for every neighbouring cell. This is illustrated in \autoref{fig:theory_003}(a) and (b). 

![Sketches illustrating the flattened mesh. \label{fig:theory_005}](figures/theory_005.png)

\autoref{fig:mesh_pic1} shows the shell mesh which is used for two out of three test cases, which are reported in the subsequent section.

![Volume mesh (a) and shell mesh (b) and domain with highlighted inlet region &#40;c) for the simulation of a radial flow experiment.\label{fig:mesh_pic1}](figures/mesh_pic1.png)



# Validation and verification

Different test cases are available, successfully validating the Julia implementation of the RTM filling model. Test cases 1 and 2 are so-called numerical permeameter experiments, see e.g. [@Fauster.2019]. The third test case is a filling study of a typical composite part which is manufactured using RTM, see [@barandun]. The testcases are:

1. Validation of the software for radial flow with isotropic in-plane permeablity: The simulated flow front position after 200 s is compared with the calculated flow front position from literature [@rudd1997liquid]. The preform size is $0.6 \times 0.6 \times 0.003$ m$^3$ with a central injection gate with 0.02 m diameter. The preform has porosity $0.70$ and isotropic permeability $3.0 \cdot 10^{−10}$ m$^2$. Resin with dynamic viscosity 0.06 Pas is injected with 35000 Pa pressure.

2. Verification of the software for radial flow with tilted orthotropic in-plane permeablity: The simulated tilted elliptical flow front is analysed and the calculated orthotropic permeablity is compared with the input (permeability $k_1$ and $k_2$, angle $\theta$ of the first principal direction to the horizontal). The preform size is $0.6 \times 0.6 \times 0.003$ m$^3$ with a central injection gate with 0.02 m diameter. The preform has porosity $0.70$ and permeability $3.0 \cdot 10^{−10}$ m$^2$ in first principal direction $30^\circ$ to the horizontal and permeability $1.5 \cdot 10^{−10}$ m$^2$ in second principal direction. Resin with dynamic viscosity 0.06 Pas is injected with 35000 Pa pressure. 

3. The filling of a complex annulus filler-like part is simulated with different software tools. The bounding box of the part is approximately $0.4 \times 0.15 \times 0.1$ m$^3$. The flow front positions at different time instances simulated with Ansys Fluent and RTMsim are compared in case 3.1. The filling patterns for same filling degrees simulated with myRTM and RTMsim are compared in case 3.2.

\autoref{fig:pic1}, \autoref{fig:pic2} and \autoref{fig:pic3} show the binary filling fraction for cases 1, 2 and 3. For every cell $i$ the fraction function is either $0$ or $1$, depending on the filling fraction being $\leq$ or $>$ than a threshold value. The threshold for the volume-of-fluid method is determined from the permeameter simulation with isotropic permeabilities (case 1). The circular and elliptical flow fronts for cases 1 and 2 are calculated from the simulated filling fraction contour plots with an image processing algorithm [@Fauster.2019]. For case 1 also a comparison with results from an analytical formula [@advani] is shown. 

A follow-up paper will discuss the physical model in detail and show the results of a validation with experimental data from a radial flow experiment with two patches with different in-plane permeability and porosity levels.

![Simulated flow front for case 1 after 200 s for a coarse and a fine mesh and flow front positions at different time instances calculated with an analytical formula and calculated from the simulated flow front with an image processing algorithm for a coarse and a fine mesh. The values for the fine mesh agree well with the results from the analytical formula. The values for the coarse mesh show an error of $\leq 15$% which decreases significantly in the course of time but the shape is no smooth circle.\label{fig:pic1}](figures/validation_pic1.png)

![Simulated flow front for case 2 after 200 s for a coarse and a fine mesh. Simulation input for the permeability was $3.0 \cdot 10^{−10}$ m$^2$ in first principal direction $30^\circ$ to the horizontal and $1.5 \cdot 10^{−10}$ m$^2$ in second principal direction. An algorithm for determining the flow front position from an optical permeameter was adapted to calculate the permeablity. The results are $2.97 \cdot 10^{−10}$ m$^2$, $1.45 \cdot 10^{−10}$ m$^2$ and angle $35^\circ$ for the fine mesh and $2.48 \cdot 10^{−10}$ m$^2$, $1.18 \cdot 10^{−10}$ m$^2$ and $35^\circ$ for the coarse mesh. For the fine mesh (2198 cells for a domain with 600 × 600 mm), this reverse engineering shows very good agreement of the calculated permeability values with the values used as simulation input. For the coarse mesh (588 cells) the agreement is still acceptable.\label{fig:pic2}](figures/validation_pic2.png)

![Results for test case 3. The results for the new simulation tool with the fine mesh are in good agreement with the ANSYS Fluent results for both filling pattern and filling time. myRTM shows similar filling patterns for same filling degrees as the new simulation tool with the coarse mesh. This comparison shows that a coarse mesh can be used for filling simulations where only the flow front propagation is investigated. If the filling time is important a mesh refinement study must be performed. \label{fig:pic3}](figures/validation_pic3.png)


# Outlook
The existing software packages for RTM mould filling simulations perform well for most scientific and engineering applications, but
there are applications where a new adaptable software package are required. For example, a transient pressure drop in the cavity after closing a valve in the hose between pressure tank and test rig can be studied with RTMsim after some extensions. Furthermore, the parts of RTMsim package can be integrated in separate codes. For example, a customized lightweight version of the solver was already used as training environment for Reinforcement learning (RL) algorithms to keep the flow front straight in a rectangular flow with reinforcement patches and three inlet ports [@stieber].

The source code is prepared for the following extensions:

- Import of mesh file in different format.

- Additional functionalities, e.g. adding temperature and degree-of-cure equations with variable resin viscosity or for VARI with variable porosity and permeability.

- Different methods for numerical differentiation and for the calculation of the numerical flux functions.


# Software download
The current release can be downloaded from https://github.com/obertscheiderfhwn/RTMsim/releases/. Installation instruction and use cases are described in the included README.md file.



# References

