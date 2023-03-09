using Glob
using Test
using rtmsim

#-----------------------
# Basic test
#-----------------------
#
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
@test return_struct.F_rho_num_add==rho_val*u_val && return_struct.F_u_num_add==rho_val*u_val*u_val && return_struct.F_gamma_num_add==gamma_val*u_val && return_struct.F_gamma_num1_add==u_val

rho_val=1.0
u_val=1.2
v_val=0.4
gamma_val=0.9
n1_val=-1.0
n2_val=0.0
A_val=1.0
input_struct=rtmsim.input_args_flux(1,[rho_val; u_val; v_val; gamma_val],[rho_val; u_val; v_val; gamma_val],[n1_val; n2_val; A_val]);return_struct=rtmsim.numerical_flux_function(input_struct);
#println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma1="*string(return_struct.F_gamma_num1_add))
@test return_struct.F_rho_num_add==-rho_val*u_val && return_struct.F_u_num_add==-rho_val*u_val*u_val && return_struct.F_gamma_num_add==-gamma_val*u_val && return_struct.F_gamma_num1_add==-u_val

rho_val=1.0
u_val=1.2
v_val=0.4
gamma_val=0.9
n1_val=0.0
n2_val=1.0
A_val=1.0
input_struct=rtmsim.input_args_flux(1,[rho_val; u_val; v_val; gamma_val],[rho_val; u_val; v_val; gamma_val],[n1_val; n2_val; A_val]);return_struct=rtmsim.numerical_flux_function(input_struct);
#println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma1="*string(return_struct.F_gamma_num1_add))
@test return_struct.F_rho_num_add==rho_val*v_val && return_struct.F_v_num_add==rho_val*v_val*v_val && return_struct.F_gamma_num_add==gamma_val*v_val && return_struct.F_gamma_num1_add==v_val

rho_val=1.0
u_val=1.2
v_val=0.4
gamma_val=0.9
n1_val=0.0
n2_val=-1.0
A_val=1.0
input_struct=rtmsim.input_args_flux(1,[rho_val; u_val; v_val; gamma_val],[rho_val; u_val; v_val; gamma_val],[n1_val; n2_val; A_val]);return_struct=rtmsim.numerical_flux_function(input_struct);
#println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma1="*string(return_struct.F_gamma_num1_add))
@test return_struct.F_rho_num_add==-rho_val*v_val && return_struct.F_v_num_add==-rho_val*v_val*v_val && return_struct.F_gamma_num_add==-gamma_val*v_val && return_struct.F_gamma_num1_add==-v_val

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
@test return_struct.F_rho_num_add==rho_val*ndotu_val && return_struct.F_u_num_add==rho_val*ndotu_val*u_val && return_struct.F_gamma_num_add==gamma_val*ndotu_val && return_struct.F_gamma_num1_add==ndotu_val

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
@test return_struct.F_rho_num_add==rho_val*ndotu_val && return_struct.F_v_num_add==rho_val*ndotu_val*v_val && return_struct.F_gamma_num_add==gamma_val*ndotu_val && return_struct.F_gamma_num1_add==ndotu_val

# Solve Burgers equation with the implemented numerical schemes
function burgers(i_method)
    # Parameters
    a=-1
    b=1
    L=b-a
    N=102
    deltax=L/(N-2)
    deltat=1e-2
    u_a=1.2
    u_b=0.4
    tmax=1.0

    #Array initialization
    rhoold=Vector{Float64}(undef, N)
    uold=Vector{Float64}(undef, N)
    vold=Vector{Float64}(undef, N)
    gammaold=Vector{Float64}(undef, N)
    rhonew=Vector{Float64}(undef, N)
    unew=Vector{Float64}(undef, N)
    vnew=Vector{Float64}(undef, N)
    gammanew=Vector{Float64}(undef, N)
    x=Vector{Float64}(undef, N)

    #Mesh
    for ind in 1:N
        x[ind]=a+(ind-3/2)*deltax;
    end

    #Initialization
    for ind in 1:N
        rhoold[ind]=1.0
        rhonew[ind]=1.0
        if x[ind]<=0.0
            uold[ind]=u_a
            unew[ind]=u_a
        else
            uold[ind]=u_b
            unew[ind]=u_b
        end
        vold[ind]=0.0
        vnew[ind]=0.0
        gammaold[ind]=0.0    
        gammanew[ind]=0.0        
    end

    #Time evolution
    iter=1
    t=deltat
    while t<=tmax           
        for ind in 2:N-1   
            if i_method==0;
                F_num_a=0.5*uold[ind-1]*uold[ind-1]
                F_num_b=0.5*uold[ind]*uold[ind]
                unew[ind]=uold[ind]-deltat/deltax*(F_num_b-F_num_a)
            elseif i_method==1
                input_struct=rtmsim.input_args_flux(1,[uold[ind-1]; uold[ind-1]; 0.0; 0.0],[uold[ind-1]; uold[ind-1]; 0.0; 0.0],[-1.0;0.0;1.0])
                return_struct=rtmsim.numerical_flux_function(input_struct)
                F_num_a=0.5*return_struct.F_rho_num_add
                input_struct=rtmsim.input_args_flux(1,[uold[ind]; uold[ind]; 0.0; 0.0],[uold[ind]; uold[ind]; 0.0; 0.0],[1.0;0.0;1.0])
                return_struct=rtmsim.numerical_flux_function(input_struct)
                F_num_b=0.5*return_struct.F_rho_num_add
                unew[ind]=uold[ind]-deltat/deltax*(F_num_b+F_num_a)       
            elseif i_method==2
                input_struct=rtmsim.input_args_flux(1,[1.0; uold[ind]; 0.0; 0.0],[1.0; uold[ind-1]; 0.0; 0.0],[-1.0;0.0;1.0])
                return_struct=rtmsim.numerical_flux_function(input_struct)
                F_num_a=0.5*return_struct.F_u_num_add
                input_struct=rtmsim.input_args_flux(1,[1.0; uold[ind]; 0.0; 0.0],[1.0; uold[ind+1]; 0.0; 0.0],[1.0;0.0;1.0])
                return_struct=rtmsim.numerical_flux_function(input_struct)
                F_num_b=0.5*return_struct.F_u_num_add
                unew[ind]=uold[ind]-deltat/deltax*(F_num_b+F_num_a)  
            elseif i_method==3
                input_struct=rtmsim.input_args_flux(1,[1.0; 0.0; uold[ind]; 0.0],[1.0; 0.0; uold[ind-1]; 0.0],[0.0;-1.0;1.0])
                return_struct=rtmsim.numerical_flux_function(input_struct)
                F_num_a=0.5*return_struct.F_v_num_add
                input_struct=rtmsim.input_args_flux(1,[1.0; 0.0; uold[ind]; 0.0],[1.0; 0.0; uold[ind+1]; 0.0],[0.0;1.0;1.0])
                return_struct=rtmsim.numerical_flux_function(input_struct)
                F_num_b=0.5*return_struct.F_v_num_add
                unew[ind]=uold[ind]-deltat/deltax*(F_num_b+F_num_a)      
            elseif i_method==4
                input_struct=rtmsim.input_args_flux(1,[1.0; uold[ind]; 0.0; uold[ind];],[1.0; uold[ind-1]; 0.0; uold[ind-1]],[-1.0;0.0;1.0])
                return_struct=rtmsim.numerical_flux_function(input_struct)
                F_num_a=0.5*return_struct.F_gamma_num_add
                input_struct=rtmsim.input_args_flux(1,[1.0; uold[ind]; 0.0; uold[ind]],[1.0; uold[ind+1]; 0.0; uold[ind+1]],[1.0;0.0;1.0])
                return_struct=rtmsim.numerical_flux_function(input_struct)
                F_num_b=0.5*return_struct.F_gamma_num_add
                unew[ind]=uold[ind]-deltat/deltax*(F_num_b+F_num_a)      
            elseif i_method==5
                input_struct=rtmsim.input_args_flux(1,[1.0; uold[ind-1]; 0.0; uold[ind-1]],[1.0; uold[ind-1]; 0.0; uold[ind-1]],[-1.0;0.0;1.0])
                return_struct=rtmsim.numerical_flux_function(input_struct)
                F_num_a=0.5*return_struct.F_gamma_num1_add*abs(return_struct.F_gamma_num1_add)
                input_struct=rtmsim.input_args_flux(1,[1.0; uold[ind]; 0.0; uold[ind]],[1.0; uold[ind]; 0.0; uold[ind]],[1.0;0.0;1.0])
                return_struct=rtmsim.numerical_flux_function(input_struct)
                F_num_b=0.5*return_struct.F_gamma_num1_add*abs(return_struct.F_gamma_num1_add)
                unew[ind]=uold[ind]-deltat/deltax*(F_num_b+F_num_a)    
            end
        end
        uold=copy(unew)
        t=t+deltat
        iter=iter+1
    end
    
    return x, unew
end
i_method=1
x,unew=burgers(i_method)  
N=length(x)
area=sum(unew[2:N-1])*(x[2]-x[1])   #Target area is 2.24
println("i_method="*string(i_method)*": area="*string(area))
@test area>=2.2 && area <=2.28
i_method=2
x,unew=burgers(i_method)  
N=length(x)
area=sum(unew[2:N-1])*(x[2]-x[1])   #Target area is 2.24
println("i_method="*string(i_method)*": area="*string(area))
@test area>=2.2 && area <=2.28
i_method=3
x,unew=burgers(i_method)  
N=length(x)
area=sum(unew[2:N-1])*(x[2]-x[1])   #Target area is 2.24
println("i_method="*string(i_method)*": area="*string(area))
@test area>=2.2 && area <=2.28
i_method=4
x,unew=burgers(i_method)  
N=length(x)
area=sum(unew[2:N-1])*(x[2]-x[1])   #Target area is 2.24
println("i_method="*string(i_method)*": area="*string(area))
@test area>=2.2 && area <=2.28
i_method=5
x,unew=burgers(i_method)  
N=length(x)
area=sum(unew[2:N-1])*(x[2]-x[1])   #Target area is 2.24
println("i_method="*string(i_method)*": area="*string(area))
@test area>=2.2 && area <=2.28


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
