# RTMsim - A Julia module for filling simulations in Resin Transfer Moulding

<img src="figures/overview.png"><br>

Resin Transfer Moulding (RTM) is a manufacturing process for producing thin-walled fiber reinforced polymer composites where dry fibers are placed inside a mould and resin is injected under pressure into the fibrous preform. During mould design, filling simulations can study different manufacturing concepts (i.e. placement of inlet ports and vents) to guarantee complete filling of the part and avoid air entrapment where flow fronts converge. 

RTMsim is a new software tool for RTM filling simulations. The porous cavity is fully described by a mesh file with triangular cells on the part’s mid-surface and cell set definitions. The latter can be used for specifying the location of the pressure injection ports and regions with different preforms by assigning different thickness, porosity and permeability values. Additional equations (e.g. for modeling the degree-of-cure) can either be added with equations of the same type or modifications of existing equations (e.g. for variable cavity thickness as needed for vacuum assisted resin infusion simulations). Several [test cases](https://obertscheiderfhwn.github.io/RTMsim/build/tutorials/)  were used for successfully validating the implemented model.


## How to get

### Requirements and installing Julia
The RTMsim module was developed and tested with Julia version 1.8. Julia is a high level open source programming language and it is as easy to use as python or Matlab. 

First of all you need a Julia installation.  Download Julia from https://julialang.org/downloads/ and install.

On a Windows operating systems add an environment variable such that the Julia terminal can be started from the command line.


### Installing using the Julia Package manager
Open a Julia terminal. The only thing you have to do is to add the package to your Julia environment with the following commands:
- `using Pkg`
- `Pkg.add(url="https://github.com/obertscheiderfhwn/RTMsim")` for the current version  or `Pkg.add(url="https://github.com/obertscheiderfhwn/RTMsim",rev="1.0.4")` for the specific release `1.0.4` corresponding to the [JOSS paper](https://joss.theoj.org/papers/ac97b5f0bc886be23981c56fe9673ca2)
- `Pkg.test("rtmsim")` for an automated test.

Alternatively, one can use the package manager with the following commands:
- Change to package manager with `]` 
- `add "https://github.com/obertscheiderfhwn/RTMsim"`
- `test rtmsim`
- Return with the `backspace` key


## How to use

### Preparation
For testing the software create a directory and download mesh- and input-files:
- Create a working directory
- Figure out the location of the package with `using rtmsim` and afterwards `pathof(rtmsim)`
- Copy from the package location the folder with the `meshfiles` and the `inputfiles` into the working directory

Start the GUI in the Julia terminal:
- In the Julia terminal change to the working directory, for example with `cd("C:\\work\\rtmsim")` where the separation is a `\\`
- `using rtmsim`
- Start the GUI with `rtmsim.gui()`

If you are working on a Windows operating system you can start RTMsim by double clicking: 
- Copy `start_rtmsim_gui.bat` and `start_rtmsim_gui.jl` from the package folder to the working directory
- Double click on `start_rtmsim_gui.bat` in the Explorer

### GUI
You can start a simulation in the GUI. 
<br><img src="figures/rtmsim_help.png"><br>

The buttons in the line on the RHS are used to start the simulation with the parameters from the selected input file. For example, click on `Select input file` and select the input.txt in the inputfiles folder from the package installation and click on `Run with input file`. After the simulation is completed, click on `Plot overview`.

Parameters (fluid properties, patch types and patch properties of cell sets specified in the mesh file) can also be specified in the GUI and a simulation is then started by clicking on `Start simulation`. Also other functionalities are available. The buttons in the first line on the LHS are used for mesh inspection, i.e. `Select mesh file`, `Plot mesh` with bounding box and `Plot sets` for inspecting the defined sets in the mesh file. The buttons `Start simulation` and `Continue simulation` in the second line on the LHS are used for starting and continuing a filling simulation. Every time the Start or Continue simulation button is pressed, a filling simulation is started. The simulated flow time `tmax` is specified in the first field in the second line. Every simulation calculates the flow front propagation during the next `tmax` seconds. If started with the `Start simulation` button, the cavity is empty initially. If started with the `Continue simulation` button, the results from the previous simulation are taken as initial condition. With the buttons `Start interactive` and `Continue interactive` in the third line one can start and continue simulations where manually selected inlet ports are used in addition to sets defined below. The radius of the inlet ports is specified in the first field in the third line and the locations are selected with the mouse after clicking on `Select inlet port`. The buttons `Plot results`, `Plot overview`and `Plot filling` in the forth line are used for creating contour plots, i.e. show filling and pressure distribution of a specified output file (path to the results file in the second cell and can be changed by clicking on `Select results file`; final results are saved in `results.jld2`), plot filling at four equidistant time instances and filling at different time instances which are selected with a slider bar. 

Click [here](https://obertscheiderfhwn.github.io/RTMsim/build/parameters) for additional information (for example the meaning of the parameters) and [here](https://obertscheiderfhwn.github.io/RTMsim/build/tutorials/) for tutorials (with typical use cases). 



## How to support and contribute
Suggestions for functionalities to be implemented and user feedback from case studies with this software are appreciated. Please have a look at the contribution item in the [community standards](https://github.com/obertscheiderfhwn/RTMsim/community). Click [here](https://obertscheiderfhwn.github.io/RTMsim/build/functions/) to see the API documentation.


## How to cite
If RTMsim is used for research or development, please use the following entry, shown in BiBTeX format:
```
@misc{RTMsim,
  author =       {Obertscheider, Christof and Fauster, Ewald},
  title =        {RTMsim - A Julia module for filling simulations in Resin Transfer Moulding},
  howpublished = {\url{https://github.com/obertscheiderfhwn/RTMsim}},
  year =         {2022}
}
```
