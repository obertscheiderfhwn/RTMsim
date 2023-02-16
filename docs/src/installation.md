# Installation

## How to get

### Requirements and installing Julia
The RTMsim module was developed and tested with Julia version 1.8. Julia is a high level open source programming language and it is as easy to use as python or Matlab. 

First of all you need a Julia installation.  Download Julia from https://julialang.org/downloads/ and install.

On a Windows operating systems add an environment variable such that the Julia terminal can be started from the command line.


### Installing using the Julia Package manager
Open a Julia terminal. The only thing you have to do is to add the package to your Julia environment with the following commands:
- `using Pkg`
- `Pkg.add(url="https://github.com/obertscheiderfhwn/RTMsim")` or `Pkg.add(url="https://github.com/obertscheiderfhwn/RTMsim",rev="1.0.4")` for specific revision `1.0.4`
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
- In the Julia terminal `cd("my\\workdirectory")` where the separation is a `\\`
- `using rtmsim
- Start the GUI with `rtmsim.gui()`

If you are working on a Windows operating system you can start RTMsim by double clicking: 
- Copy `start_rtmsim_gui.bat` and `start_rtmsim_gui.jl` from the package folder to the working directory
- Double click on `start_rtmsim_gui.bat` in the Explorer

### GUI
You can start a simulation in the GUI. 
```@raw html
<br><img src="../assets/figures/rtmsim_help.png"><br>
```
The buttons in the first line on the LHS are used for mesh inspection, i.e. select a mesh file, plot the mesh with bounding box and plot the defined sets. The buttons in the second line on the LHS are used for starting and continuing a filling simulation. Every time the Start or Continue simulation button is pressed, a filling simulation is started. The simulated flow time `tmax`, the patch types and patch properties must be specified before. Every simulation calculates the flow front propagation during the next `tmax` seconds. If started with the Start simulation button, the cavity is empty initially. If started with the Continue simulation button, the results from the previous simulation are taken as initial condition. With the buttons in the third line one can select inlet ports with specified radius interactively in addition to using the defined sets, and start and continue such a simulation. The buttons in the forth line are used for post-processing, i.e. show filling and pressure distribution of a specified output file (final results are saved in results.jld2), plot filling at four equidistant time instances and filling at different time instances which are selected with a slider bar. The buttons in the line on the RHS are used to start the simulation with the parameters from the selected input file.









```@meta
EditURL = "https://github.com/obertscheiderfhwn/RTMsim/blob/main/docs/src/installation.md"
```