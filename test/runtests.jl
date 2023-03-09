using Glob
using Test
using rtmsim

#-----------------------
# Basic test
#-----------------------
#
@testset "numerical flux functions" begin
    println(" ")
    println("Basic testing...")
    println()
    println("...for the numerical flux functions")
    #check consistency if boundary face is aligned in principal directions

    rho_val=1.0
    u_val=1.2
    v_val=0.4
    gamma_val=0.9
    n1_val=1.0
    n2_val=0.0
    A_val=1.0
    input_struct=rtmsim.input_args_flux(1,[rho_val; u_val; v_val; gamma_val],[rho_val; u_val; v_val; gamma_val],[n1_val; n2_val; A_val]);return_struct=rtmsim.numerical_flux_function(input_struct);
    #println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma1="*string(return_struct.F_gamma_num1_add))
    @test return_struct.F_rho_num_add==rho_val*u_val 
    @test return_struct.F_u_num_add==rho_val*u_val*u_val 
    @test return_struct.F_gamma_num_add==gamma_val*u_val
    @test return_struct.F_gamma_num1_add==u_val

    rho_val=1.0
    u_val=1.2
    v_val=0.4
    gamma_val=0.9
    n1_val=-1.0
    n2_val=0.0
    A_val=1.0
    input_struct=rtmsim.input_args_flux(1,[rho_val; u_val; v_val; gamma_val],[rho_val; u_val; v_val; gamma_val],[n1_val; n2_val; A_val]);return_struct=rtmsim.numerical_flux_function(input_struct);
    #println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma1="*string(return_struct.F_gamma_num1_add))
    @test return_struct.F_rho_num_add==-rho_val*u_val
    @test return_struct.F_u_num_add==-rho_val*u_val*u_val 
    @test return_struct.F_gamma_num_add==-gamma_val*u_val
    @test return_struct.F_gamma_num1_add==-u_val

    rho_val=1.0
    u_val=1.2
    v_val=0.4
    gamma_val=0.9
    n1_val=0.0
    n2_val=1.0
    A_val=1.0
    input_struct=rtmsim.input_args_flux(1,[rho_val; u_val; v_val; gamma_val],[rho_val; u_val; v_val; gamma_val],[n1_val; n2_val; A_val]);return_struct=rtmsim.numerical_flux_function(input_struct);
    #println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma1="*string(return_struct.F_gamma_num1_add))
    @test return_struct.F_rho_num_add==rho_val*v_val
    @test return_struct.F_v_num_add==rho_val*v_val*v_val
    @test return_struct.F_gamma_num_add==gamma_val*v_val
    @test return_struct.F_gamma_num1_add==v_val

    rho_val=1.0
    u_val=1.2
    v_val=0.4
    gamma_val=0.9
    n1_val=0.0
    n2_val=-1.0
    A_val=1.0
    input_struct=rtmsim.input_args_flux(1,[rho_val; u_val; v_val; gamma_val],[rho_val; u_val; v_val; gamma_val],[n1_val; n2_val; A_val]);return_struct=rtmsim.numerical_flux_function(input_struct);
    #println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma1="*string(return_struct.F_gamma_num1_add))
    @test return_struct.F_rho_num_add==-rho_val*v_val
    @test return_struct.F_v_num_add==-rho_val*v_val*v_val
    @test return_struct.F_gamma_num_add==-gamma_val*v_val
    @test return_struct.F_gamma_num1_add==-v_val

    rho_val=1.0
    u_val=1.2
    v_val=0.4
    gamma_val=0.9
    n1_val=1.0
    n2_val=0.0
    A_val=1.0
    ndotu_val=1.6
    input_struct=rtmsim.input_args_flux_boundary(1,[rho_val; u_val; v_val; gamma_val],[rho_val; u_val; v_val; gamma_val],[n1_val; n2_val; A_val], ndotu_val);return_struct=rtmsim.numerical_flux_function_boundary(input_struct);
    #println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma1="*string(return_struct.F_gamma_num1_add))
    @test return_struct.F_rho_num_add==rho_val*ndotu_val
    @test return_struct.F_u_num_add==rho_val*ndotu_val*u_val
    @test return_struct.F_gamma_num_add==gamma_val*ndotu_val
    @test return_struct.F_gamma_num1_add==ndotu_val

    rho_val=1.0
    u_val=1.2
    v_val=0.4
    gamma_val=0.9
    n1_val=0.0
    n2_val=1.0
    A_val=1.0
    ndotu_val=1.6
    input_struct=rtmsim.input_args_flux_boundary(1,[rho_val; u_val; v_val; gamma_val],[rho_val; u_val; v_val; gamma_val],[n1_val; n2_val; A_val], ndotu_val);return_struct=rtmsim.numerical_flux_function_boundary(input_struct);
    #println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma1="*string(return_struct.F_gamma_num1_add))
    @test return_struct.F_rho_num_add==rho_val*ndotu_val
    @test return_struct.F_v_num_add==rho_val*ndotu_val*v_val 
    @test return_struct.F_gamma_num_add==gamma_val*ndotu_val
    @test return_struct.F_gamma_num1_add==ndotu_val
end

@testset "bugers" begin
    i_method=1
    x,unew=rtmsim.burgers(i_method)  
    N=length(x)
    area=sum(unew[2:N-1])*(x[2]-x[1])   #Target area is 2.24
    println("i_method=$(i_method) area=$(area)")
    @test 2.2 <= area <=2.28
    i_method=2
    x,unew=rtmsim.burgers(i_method)  
    N=length(x)
    area=sum(unew[2:N-1])*(x[2]-x[1])   #Target area is 2.24
    println("i_method=$(i_method) area=$(area)")
    @test 2.2 <= area <=2.28
    i_method=3
    x,unew=rtmsim.burgers(i_method)  
    N=length(x)
    area=sum(unew[2:N-1])*(x[2]-x[1])   #Target area is 2.24
    println("i_method=$(i_method) area=$(area)")
    @test 2.2 <= area <=2.28
    i_method=4
    x,unew=rtmsim.burgers(i_method)  
    N=length(x)
    area=sum(unew[2:N-1])*(x[2]-x[1])   #Target area is 2.24
    println("i_method=$(i_method) area=$(area)")
    @test 2.2 <= area <=2.28
    i_method=5
    x,unew=rtmsim.burgers(i_method)  
    N=length(x)
    area=sum(unew[2:N-1])*(x[2]-x[1])   #Target area is 2.24
    println("i_method=$(i_method) area=$(area)")
    @test 2.2 <= area <=2.28
end

println(" ")
println("...for running a simulation")
MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]
meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf")
#print(meshfilename*" \n")
param=rtmsim.input_vals(1,meshfilename,200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2,0,"results.jld2",0,0.01,16);
rtmsim.rtmsim_rev1(param)

println(" ")
println("...for plotting the results")
WORK_DIR=pwd();
resultsfilename=joinpath(WORK_DIR,"results.jld2")
#print(resultsfilename*" \n")
rtmsim.plot_results(resultsfilename);

println(" ")
println("...for deleting output files")
rm.(glob("*.jld2",WORK_DIR));
println(" ")