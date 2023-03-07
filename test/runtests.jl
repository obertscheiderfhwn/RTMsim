using Glob
using rtmsim

#-----------------------
# Basic test
#-----------------------
#
println("Basic testing...")
println("...for the numerical flux functions")
input_struct=rtmsim.input_args_flux(1,[1.0; 1.2; 0.0; 0.0],[1.0; 0.4; 0.0; 0.0],[1.0;0.0;1.0]);return_struct=rtmsim.numerical_flux_function(input_struct);println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add))
input_struct=rtmsim.input_args_flux(1,[1.225; 1.2; 0.4; 0.9],[1.0; 0.4; 1.2; 0.1],[1/sqrt(2);1/sqrt(2);1.0]);return_struct=rtmsim.numerical_flux_function(input_struct);println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma_num1="*string(return_struct.F_gamma_num1_add))
input_struct=rtmsim.input_args_flux_boundary(1,[1.0; 1.2; 0.0; 0.0],[1.0; 0.0; 0.0; 0.0],[1.0;0.0;1.0],1.2);return_struct=rtmsim.numerical_flux_function(input_struct);println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add))
input_struct=rtmsim.input_args_flux_boundary(1,[1.225; 1.2; 0.4; 0.9],[1.0; 0.4; 1.2; 0.1],[1/sqrt(2);1/sqrt(2);1.0],-0.8);return_struct=rtmsim.numerical_flux_function_boundary(input_struct);println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma_num1="*string(return_struct.F_gamma_num1_add))

println("...for running a simulation")
MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]
meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf")
#print(meshfilename*" \n")
param=rtmsim.input_vals(1,meshfilename,200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2,0,"results.jld2",0,0.01,16);
rtmsim.rtmsim_rev1(param)

println("...for plotting the results")
WORK_DIR=pwd();
resultsfilename=joinpath(WORK_DIR,"results.jld2")
#print(resultsfilename*" \n")
rtmsim.plot_results(resultsfilename);

println("...for deleting output files")
rm.(glob("*.jld2",WORK_DIR));
