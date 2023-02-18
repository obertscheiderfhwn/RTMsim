using Glob
using rtmsim

#-----------------------
# Basic test
#-----------------------
#
println("Basic testing...")
println("...for running a simulation")
MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]  #path of folder of current package
meshfilename=string(MODULE_ROOT,"\\meshfiles\\mesh_permeameter1_foursets.bdf")   #meshfile in meshfiles folder of package 
rtmsim.rtmsim_rev1(1,meshfilename,200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2,0,"results.jld2",0,0.01,16)
println("...for plotting the results")
rtmsim.plot_results("results.jld2")
rm.(glob("*.jld2"))
