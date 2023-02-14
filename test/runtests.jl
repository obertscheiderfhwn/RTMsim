using Glob
using rtmsim

#-----------------------
# Basic test
#-----------------------
#
println("Basic testing...")
println("...for running a simulation")
rtmsim.rtmsim_rev1(1,"..\\meshfiles\\mesh_permeameter1_foursets.bdf",200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2,0,"results.jld2",0,0.01,16)
println("...for plotting the results")
rtmsim.plot_results("results.jld2")
rm.(glob("*.jld2"))
