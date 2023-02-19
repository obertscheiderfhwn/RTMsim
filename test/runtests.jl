using Glob
using rtmsim

#-----------------------
# Basic test
#-----------------------
#
println("Basic testing...")
println("...for running a simulation")
MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]
meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf")
#print(meshfilename*" \n")
rtmsim.rtmsim_rev1(1,meshfilename,200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2,0,"results.jld2",0,0.01,16)

println("...for plotting the results")
WORK_DIR=pwd();
resultsfilename=joinpath(WORK_DIR,"results.jld2")
#print(resultsfilename*" \n")
rtmsim.plot_results(resultsfilename);

println("...for deleting output files")
rm.(glob("*.jld2",WORK_DIR));
