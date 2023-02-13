using Glob
using rtmsim

#-----------------------
# Basic test
#-----------------------
#
# for simulating and plotting results
rtmsim.rtmsim_rev1(1,"..\\meshfiles\\mesh_permeameter1_foursets.bdf",200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2,0,"results.jld2",0,0.01,16)
rtmsim.plot_results("results.jld2")
rm.(glob("*.jld2"))

#=
#-----------------------
# Unit tests
#-----------------------
#
# for plotting the mesh defined in the bdf-file
rtmsim.plot_mesh("..\\meshfiles\\mesh_permeameter1_foursets.bdf",1) 
#
# for plotting the sets specified in the bdf-file
rtmsim.plot_sets("..\\meshfiles\\mesh_permeameter1_foursets.bdf") 
#
# for starting a simulation with one pressure inlet port (sets 2, 3 and 4 are not used and consequently the preform parameters are ignored; since set 1 is a pressure inlet, also the parameters for set 1 are ignored and the only relevant parameter for the specified set is the pressure difference between injection and initial cavity pressure)
rtmsim.rtmsim_rev1(1,"..\\meshfiles\\mesh_permeameter1_foursets.bdf",200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,0,0,0, 0,"results.jld2",0,0.01,16) 
#
# for starting a simulation with different patches and race tracking
rtmsim.rtmsim_rev1(1,"..\\meshfiles\\mesh_permeameter1_foursets.bdf",200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2, 0,"results.jld2",0,0.01,16) 
#
# for continuing the previous simulation
rtmsim.rtmsim_rev1(1,"..\\meshfiles\\mesh_permeameter1_foursets.bdf",200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2, 1,"results.jld2",0,0.01,16) 
#
# for the manual selection of inlet ports with left mouse button click while key p is pressed
#rtmsim.plot_mesh("..\\meshfiles\\mesh_annulusfiller1.bdf",2) 
#
# for starting only with the interactively selected inlet ports
#rtmsim.rtmsim_rev1(1,"..\\meshfiles\\mesh_annulusfiller1.bdf",200, 0.35e5,1.205,1.4,0.06, 0.35e5,0.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 0,0,0,0, 0,"results.jld2",1,0.01,16) 
#
# for plotting the final filling and pressure contours
rtmsim.plot_results("results.jld2") 
#
# for plotting the filling contours at four equidistant time instances
rtmsim.plot_overview(-1,-1) 
#
# for plotting the filling at different time instances selected with a slider bar
rtmsim.plot_filling(-1,-1) 
#
# for starting a simulation with the parameters specified in the text file input.txt
rtmsim.start_rtmsim("..\\inputfiles\\input.txt") 
#
# 
=#