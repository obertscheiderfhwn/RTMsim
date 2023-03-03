"""
    rtmsim_rev1(param)

RTMsim solver with the following main steps:
- Simulation initialization
- Read mesh file and prepare patches  
- Find neighbouring cells
- Assign parameters to cells
- Create local cell coordinate systems
- Calculate initial time step
- Array initialization
- Define simulation time and intermediate output times
- Boundary conditions
- (Optional initialization if `i_model=2,3,..`)
- Time evolution (for loops over all indices inside a while loop for time evolution)
    - Calculation of correction factors for cell thickness, porosity, permeability, viscosity
    - Pressure gradient calculation
    - Numerical flux function calculation
    - Update of rho, u, v, gamma and p according to conservation laws and equation of state
    - Boundary conditions
    - Prepare arrays for next time step
    - Saving of intermediate data
    - (Opional time marching etc. for i_model=2,3,...)
    - Calculation of adaptive time step 

Arguments: Data structure `param` with
- `i_model :: Int64 =1`
- `meshfilename :: String ="input.txt"`
- `tmax :: Float64 =200`
- `p_ref :: Float64 =1.01325e5`
- `rho_ref :: Float64 =1.225`
- `gamma :: Float64 =1.4`
- `mu_resin_val :: Float64 =0.06`
- `p_a_val :: Float64 =1.35e5`
- `p_init_val :: Float64 =1.0e5`  
- `t_val :: Float64 =3e-3`
- `porosity_val :: Float64 =0.7`
- `K_val :: Float64 =3e-10`
- `alpha_val :: Float64 =1.0`
- `refdir1_val :: Float64 =1.0`
- `refdir2_val :: Float64 =0.0`
- `refdir3_val :: Float64 =0.0`
- `t1_val :: Float64 =3e-3`
- `porosity1_val :: Float64 =0.7`
- `K1_val :: Float64 =3e-10`
- `alpha1_val :: Float64 =1.0`
- `refdir11_val :: Float64 =1.0`
- `refdir21_val :: Float64 =0.0`
- `refdir31_val :: Float64 =0.0`
- `t2_val :: Float64 =3e-3`
- `porosity2_val :: Float64 =0.7`
- `K2_val :: Float64 =3e-10`
- `alpha2_val :: Float64 =1.0`
- `refdir12_val :: Float64 =1.0`
- `refdir22_val :: Float64 =0.0`
- `refdir32_val :: Float64 =0.0`
- `t3_val :: Float64 =3e-3`
- `porosity3_val :: Float64 =0.7`
- `K3_val :: Float64 =3e-10`
- `alpha3_val :: Float64 =1.0`
- `refdir13_val :: Float64 =0.0`
- `refdir23_val :: Float64 =0.0`
- `refdir33_val :: Float64 =0.0`
- `t4_val :: Float64 =3e-3`
- `porosity4_val :: Float64 =0.7`
- `K4_val :: Float64 =3e-10`
- `alpha4_val :: Float64 =1.0`
- `refdir14_val :: Float64 =1.0`
- `refdir24_val :: Float64 =0.0`
- `refdir34_val :: Float64 =0.0`
- `patchtype1val :: Int64 =1`
- `patchtype2val :: Int64 =0`
- `patchtype3val :: Int64 =0`
- `patchtype4val :: Int64 =0`
- `i_restart :: Int64 =0`
- `restartfilename :: String ="results.jld2"`
- `i_interactive :: Int64 =0`
- `r_p :: Float64 =0.01`
- `n_pics :: Int64 =16`

Meaning of the variables:
- `i_model`: Identifier for physical model (Default value is 1)
- `meshfilename`: Mesh filename.
- `tmax`: Maximum simulation time.
- `p_ref rho_ref gamma mu_resin_val`: Parameters for the adiabatic equation of state and dynamic viscosity of resin used in the Darcy term.
- `p_a_val p_init_val `: Absolut pressure value for injection port and for initial cavity pressure.
- `t_val porosity_val K_val alpha_val refdir1_val refdir2_val refdir3_val`: Properties of the cells in the main preform: The vector `(refdir1_val,refdir2_val,refdir3_val)` is projected onto the cell in order to define the first principal cell direction. The second principal cell direction is perpendicular to the first one in the plane spanned by the cell nodes. The principal cell directions are used as the principal permeabilty directions. The cell properties are defined by the thickness `t_val`, the porosity `porosity_val`, the permeability `K_val` in the first principal cell direction, the permeablity `alpha_val` in the second principal direction.
- `t1_val porosity1_val K1_val alpha1_val refdir11_val refdir21_val refdir31_val` etc.: Properties for up to four additional cell regions if preform. 
- `patchtype1val patchtype2val patchtype3val patchtype4val`: These regions are used to specify the location of the pressure boundary conditions and to specify regions with different permeability, porosity and thickness properties (e.g. for different part thickness and layup or for race tracking which are regions with very high permeability typically at the boundary of the preforms). Vents need not be specified. Parameters `patchtype1val` define the patch type. Numerical values 0, 1, 2 and 3 are allowed with the following interpretation:
    - 0 .. the patch is ignored
    - 1 .. the patch represents an inlet gate, where the specified injection pressure level applies
    - 2 .. the patch specifies a preform region
    - 3 .. the patch represents a vent, where the specified initial pressure level applies
- `i_restart restartfilename`: Start with new simulation if `0` or continue previous simulation if `1` from specified file
- `i_interactive r_p`: Select the inlet ports graphically if i_interactive equal to `1` and inlet ports have specified radius
- `n_pics`: Number of intermediate output files, supposed to be a multiple of `4`

Unit tests:
- `MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]; meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf"); param=rtmsim.input_vals(1,meshfilename,200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2,0,"results.jld2",0,0.01,16); rtmsim.rtmsim_rev1(param);`
    
Addtional unit tests:
- `MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]; meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf"); param=rtmsim.input_vals(1,meshfilename,200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,0,0,0, 0,"results.jld2",0,0.01,16); rtmsim.rtmsim_rev1(param);` for starting a simulation with one pressure inlet port (sets 2, 3 and 4 are not used and consequently the preform parameters are ignored; since set 1 is a pressure inlet, also the parameters for set 1 are ignored and the only relevant parameter for the specified set is the pressure difference between injection and initial cavity pressure)
- `MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]; meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf"); param=rtmsim.input_vals(1,meshfilename,200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2, 0,"results.jld2",0,0.01,16); rtmsim.rtmsim_rev1(param);` for starting a simulation with different patches and race tracking
- `MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]; meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf"); param=rtmsim.input_vals(1,meshfilename,200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2, 1,"results.jld2",0,0.01,16); rtmsim.rtmsim_rev1(param);` for continuing the previous simulation   
"""
function rtmsim_rev1(param)
    
    if Sys.iswindows()
        meshfilename=replace(param.meshfilename,"/" => "\\")
    elseif Sys.islinux()
        meshfilename=replace(param.meshfilename,"\\" => "/")
    end  

    #----------------------------------------------------------------------
    # Simulation initialization
    #----------------------------------------------------------------------
     
    # Well defined variable types, except for strings meshfilename,restartfilename
    i_model=param.i_model
    tmax=param.tmax
    p_ref=param.p_ref;rho_ref=param.rho_ref;gamma=param.gamma;mu_resin_val=param.mu_resin_val
    p_a_val=param.p_a_val;p_init_val=param.p_init_val
    t_val=param.t_val;porosity_val=param.porosity_val;K_val=param.K_val;alpha_val=param.alpha_val;refdir1_val=param.refdir1_val;refdir2_val=param.refdir2_val;refdir3_val=param.refdir3_val
    t1_val=param.t1_val;porosity1_val=param.porosity1_val;K1_val=param.K1_val;alpha1_val=param.alpha1_val;refdir11_val=param.refdir11_val;refdir21_val=param.refdir21_val;refdir31_val=param.refdir31_val
    t2_val=param.t2_val;porosity2_val=param.porosity2_val;K2_val=param.K2_val;alpha2_val=param.alpha2_val;refdir12_val=param.refdir12_val;refdir22_val=param.refdir22_val;refdir32_val=param.refdir32_val
    t3_val=param.t3_val;porosity3_val=param.porosity3_val;K3_val=param.K3_val;alpha3_val=param.alpha3_val;refdir13_val=param.refdir13_val;refdir23_val=param.refdir23_val;refdir33_val=param.refdir33_val
    t4_val=param.t4_val;porosity4_val=param.porosity4_val;K4_val=param.K4_val;alpha4_val=param.alpha4_val;refdir14_val=param.refdir14_val;refdir24_val=param.refdir24_val;refdir34_val=param.refdir34_val
    patchtype1val=param.patchtype1val;patchtype2val=param.patchtype2val;patchtype3val=param.patchtype3val;patchtype4val=param.patchtype4val
    i_restart=param.i_restart;i_interactive=param.i_interactive;n_pics=param.n_pics;r_p=param.r_p;restartfilename=param.restartfilename

    #License statement
    println()
    println("RTMsim version 1.0")
    println()
    println("RTMsim is Julia code with GUI which simulates the mold filling in Liquid Composite Molding (LCM) manufacturing process.")
    println("Copyright (C) 2023 Christof Obertscheider / University of Applied Sciences Wiener Neustadt (FHWN)")    
    println()
    println("This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.")
    println()
    println("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.")
    println()
    println("You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.")
    println()
    println("This software is free of charge and may be used for commercial and academic purposes.  Please mention the use of this software at an appropriate place in your work.")
    println()
    println("Submit bug reports to christof.obertscheider@fhwn.ac.at ")
    println()

    #Output simulation parameter overview
    println()
    println("RTMsim started with the following parameters:")
    println("i_model=",i_model)
        if i_model!=1
            errorstring=string("Only iso-thermal RTM implemented, i.e. i_model must be =1 instead of =",string(i_model)*"\n") 
            error(errorstring)
        end
    if Sys.iswindows()
        meshfilename=replace(meshfilename,"/" => "\\")
    elseif Sys.islinux()
        meshfilename=replace(meshfilename,"\\" => "/")
    end  
    println("meshfilename=",meshfilename)
        if ~isfile(meshfilename)
            errorstring=string("File ",meshfilename," not existing"* "\n") 
            error(errorstring)
        end
    println("tmax=",string(tmax))
        if tmax<=0.0
            errorstring="tmax must be greater than zero"
            error(errorstring)
        end
        #Limit number of results time steps between n_pics_min and n_pics_max and make it multiple of 4
        n_pics_input=n_pics
        n_pics_min=Int64(4)
        n_pics_max=Int64(100)
        if mod(n_pics,4)!=0
            n_pics=round(n_pics/4)*4
        end
        if n_pics<n_pics_min
            n_pics=n_pics_min
        end
        if n_pics>n_pics_max
            n_pics=n_pics_max 
        end       
        if n_pics>n_pics_max
            n_pics=n_pics_max
        end
        n_pics=Int64(n_pics)
    if n_pics_input!=n_pics
        println("n_pics changed to n_pics=",string(n_pics))
    else
        println("n_pics=",string(n_pics),)
    end
    println("i_interactive=",string(i_interactive))   
        if i_interactive!=0 && i_interactive!=1 && i_interactive!=2
            errorstring="Wrong value for i_interactive (must be=0,1,2)"
            error(errorstring)
        end 
    if i_restart==1 
        println("i_restart,restartfilename=",string(i_restart), ",", restartfilename) 
        if i_restart!=0 && i_restart!=1 
            errorstring="Wrong value for i_restart (must be=0,1)"
            error(errorstring)
        end 
        if ~isfile(restartfilename)
            errorstring=string("File ",restartfilename," not existing"* "\n") 
            error(errorstring)
        end
    end
    println("p_ref,rho_ref,gamma,mu=", string(p_ref), ",", string(rho_ref), ",", string(gamma), ",", string(mu_resin_val))
    if p_ref<=0.0 || rho_ref<=0.0 || gamma<1.0 || mu_resin_val<=0.0
        errorstring="Wrong value for p_ref,rho_ref,gamma,mu (must be >0.0,>0.0,>1.0,>0.0)"
        error(errorstring)
    end 
    println("p_a_val,p_init_val=", string(p_a_val), ",", string(p_init_val))
        if p_a_val<=p_init_val
            errorstring="Injection pressure must be greater than initial pressure"
            error(errorstring)
        end
        if p_a_val<=0.0 || p_init_val<0.0
            errorstring="Wrong value for p_a_val,p_init_val (must be >0.0,>0.0)"
            error(errorstring)
        end 

    #Maximum number of cell neighbours
    maxnumberofneighbours=10

    #Delete old files and abort if meshfile is not existing
    if ~isfile(meshfilename)
        errorstring=string("File ",meshfilename," not existing"* "\n") 
        error(errorstring)
    end
    if i_restart==1
        cp(restartfilename,"restart.jdl2";force=true)
    end
    delete_files()     
    if i_restart==1
        cp("restart.jdl2",restartfilename;force=true)
    end
    n_out=Int64(0)

    #Assign and prepare physical parameters
    refdir_val=[refdir1_val,refdir2_val,refdir3_val]  #Vector
    u_a=0.0  
    u_b=0.0 
    u_init=0.0 
    v_a=0.0 
    v_b=0.0 
    v_init=0.0 
    p_a=p_a_val
    p_init=p_init_val
    p_b=p_a_val
    #Normalization for Delta p: p->p-p_init
        p_eps=0.001e5 #0.000e5  #
        p_a=p_a-p_init+p_eps
        # p_init=p_init-p_init+p_eps
        p_init=p_eps
        p_b=p_a-p_init+p_eps
        p_ref=p_ref  #p_ref-p_init+p_eps
    kappa=p_ref/(rho_ref^gamma)
    #Lookuptable for adiabatic law (required for stability)
        p_int1=0.0e5; rho_int1=(p_int1/kappa)^(1/gamma)
        p_int2=0.1e5; rho_int2=(p_int2/kappa)^(1/gamma)
        p_int3=0.5e5; rho_int3=(p_int3/kappa)^(1/gamma)
        p_int4=1.0e5; rho_int4=(p_int4/kappa)^(1/gamma)
        p_int5=10.0e5; rho_int5=(p_int5/kappa)^(1/gamma)
        p_int6=100.0e5; rho_int6=(p_int6/kappa)^(1/gamma)
        A=[rho_int1^2 rho_int1 1.0; rho_int3^2 rho_int3 1.0; rho_int4^2 rho_int4 1.0]
        b=[p_int1;p_int3;p_int4]
        apvals=A\b
        ap1=apvals[1];ap2=apvals[2];ap3=apvals[3]
    rho_a=(p_a/kappa)^(1.0/gamma)
    rho_b=(p_b/kappa)^(1.0/gamma)
    rho_init=(p_init/kappa)^(1.0/gamma)

    if gamma>=100  #insert here coefficients for an incompressible EOS with resin mass density as rho_ref 
                    #at p_b and 0.9*rho_ref at p_a but the EOS is for deltap and consequently normalized pressure values
        #ap1=0
        #ap2=(p_b-p_init)/(0.1*rho_ref)
        #ap3=p_b-(p_b-p_init)/0.1 
        #rho_a=p_a/ap2-ap3/ap2
        #rho_b=p_b/ap2-ap3/ap2
        #rho_init=p_init/ap2-ap3/ap2

        rho_a=rho_ref
        rho_b=rho_a
        rho_init=0.0
        p_int1=p_init; rho_int1=rho_init
        p_int2=p_init+0.9*(p_a-p_init); rho_int2=0.1*rho_a
        p_int3=p_a; rho_int3=rho_a
        A=[rho_int1^2 rho_int1 1.0; rho_int3^2 rho_int3 1.0; 2*rho_int3 1.0 0]
        b=[p_int1;p_int3;0.0]
        apvals=A\b
        ap1=apvals[1];ap2=apvals[2];ap3=apvals[3]
        #ap1*rho_new[ind]^2+ap2*rho_new[ind]+ap3


        println(string("rho_int1: ",string(rho_int1) ) )
        println(string("rho_int2: ",string(rho_int2) ) )
        println(string("rho_int3: ",string(rho_int3) ) )
        println(string("p_int1: ",string(p_int1) ) )
        println(string("p_int2: ",string(p_int2) ) )
        println(string("p_int3: ",string(p_int3) ) )
        
        println(string("ap1: ",string(ap1) ) )
        println(string("ap2: ",string(ap2) ) ) 
        println(string("ap3: ",string(ap3) ) ) 
        println(string("p_a: ",string(p_a) ) )
        println(string("p_b: ",string(p_b) ) )
        println(string("p_init: ",string(p_init) ) )
        println(string("rho_a: ",string(rho_a) ) )
        println(string("rho_b: ",string(rho_b) ) )
        println(string("rho_init: ",string(rho_init) ) )
    end
    
    T_a=295.
    T_b=295.
    T_init=295.
    gamma_a=1.0
    gamma_b=1.0    
    gamma_init=0.0
    paramset=[porosity_val,t_val,K_val,alpha_val,refdir1_val,refdir2_val,refdir3_val]
    paramset1=[porosity1_val,t1_val,K1_val,alpha1_val,refdir11_val,refdir21_val,refdir31_val]
    paramset2=[porosity2_val,t2_val,K2_val,alpha2_val,refdir12_val,refdir22_val,refdir32_val]
    paramset3=[porosity3_val,t3_val,K3_val,alpha3_val,refdir13_val,refdir23_val,refdir33_val]
    paramset4=[porosity4_val,t4_val,K4_val,alpha4_val,refdir14_val,refdir24_val,refdir34_val]


    #--------------------------------------------------------------------------
    # Read mesh file and prepare patches     
    #--------------------------------------------------------------------------
    input_mesh=rtmsim.input_args_read_mesh(meshfilename,paramset,paramset1,paramset2,paramset3,paramset4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_interactive,r_p)
    return_mesh=read_mesh(input_mesh)
    N=return_mesh.N
    cellgridid=return_mesh.cellgridid
    gridx=return_mesh.gridx
    gridy=return_mesh.gridy
    gridz=return_mesh.gridz
    cellcenterx=return_mesh.cellcenterx
    cellcentery=return_mesh.cellcentery
    cellcenterz=return_mesh.cellcenterz
    patchparameters=return_mesh.patchparameters
    patchparameters1=return_mesh.patchparameters1
    patchparameters2=return_mesh.patchparameters2
    patchparameters3=return_mesh.patchparameters3
    patchparameters4=return_mesh.patchparameters4
    patchids1=return_mesh.patchids1
    patchids2=return_mesh.patchids2
    patchids3=return_mesh.patchids3
    patchids4=return_mesh.patchids4
    inletpatchids=return_mesh.inletpatchids
    println(string("parameters for main preform: ",string(patchparameters) ) )
    if patchparameters[1]<=0.0 || patchparameters[1]>1.0 || patchparameters[2]<=0.0 || patchparameters[3]<=0 || patchparameters[4]<=0
        errorstring="Wrong value for porosity,thickness,permeability,alpha (must be between >0 and <=1,>0.0,>0.0,>0.0)"
        error(errorstring)
    end 
    if ~isempty(patchids1) 
        if patchtype1val==1        
            println("patch 1 is pressure inlet") 
        elseif patchtype1val==2        
            println(string("parameters for patch 1: ",string(patchparameters1)) )
            if patchparameters1[1]<=0.0 || patchparameters1[1]>1.0 || patchparameters1[2]<=0.0 || patchparameters1[3]<=0 || patchparameters1[4]<=0
                errorstring="Wrong value for porosity,thickness,permeability,alpha (must be between >0 and <=1,>0.0,>0.0,>0.0)"
                error(errorstring)
            end 
        elseif patchtype1val==3        
            println("patch 1 is pressure outlet") 
        end
    end
    if ~isempty(patchids2) 
        if patchtype2val==1        
            println("patch 4 is pressure inlet") 
        elseif patchtype2val==2        
            println(string("parameters for patch 2: ",string(patchparameters2)) )
            if patchparameters2[1]<=0.0 || patchparameters2[1]>1.0 || patchparameters2[2]<=0.0 || patchparameters2[3]<=0 || patchparameters2[4]<=0
                errorstring="Wrong value for porosity,thickness,permeability,alpha (must be between >0 and <=1,>0.0,>0.0,>0.0)"
                error(errorstring)
            end 
        elseif patchtype2val==3        
            println("patch 2 is pressure outlet") 
        end
    end
    if ~isempty(patchids3) 
        if patchtype3val==1        
            println("patch 3 is pressure inlet") 
        elseif patchtype3val==2        
            println(string("parameters for patch 3: ",string(patchparameters3)) )
            if patchparameters3[1]<=0.0 || patchparameters3[1]>1.0 || patchparameters3[2]<=0.0 || patchparameters3[3]<=0 || patchparameters3[4]<=0
                errorstring="Wrong value for porosity,thickness,permeability,alpha (must be between >0 and <=1,>0.0,>0.0,>0.0)"
                error(errorstring)
            end 
        elseif patchtype3val==3        
            println("patch 3 is pressure outlet") 
        end
    end
    if ~isempty(patchids4) 
        if patchtype4val==1        
            println("patch 4 is pressure inlet") 
        elseif patchtype4val==2        
            println(string("parameters for patch 4: ",string(patchparameters4)) )
            if patchparameters4[1]<=0.0 || patchparameters4[1]>1.0 || patchparameters4[2]<=0.0 || patchparameters4[3]<=0 || patchparameters4[4]<=0
                errorstring="Wrong value for porosity,thickness,permeability,alpha (must be between >0 and <=1,>0.0,>0.0,>0.0)"
                error(errorstring)
            end 
        elseif patchtype4val==3        
            println("patch 4 is pressure outlet")   
        end
    end
    if patchtype1val!=1 && patchtype2val!=1 && patchtype3val!=1 && patchtype3val!=1 && i_interactive==0 && i_restart==0
        errorstring=string("No inlet defined" * "\n") 
        error(errorstring)
    end
    if i_interactive==1 || i_interactive==2
        println("additional inlet defined interactively")   
    end

    #--------------------------------------------------------------------------
    #  Find neighbouring cells
    #--------------------------------------------------------------------------    
    input_faces=rtmsim.input_args_create_faces(cellgridid, N, maxnumberofneighbours);
    return_faces=create_faces(input_faces)
    faces=return_faces.faces
    cellneighboursarray=return_faces.cellneighboursarray
    celltype=return_faces.celltype

    #--------------------------------------------------------------------------
    #  Assign parameters to cells
    #--------------------------------------------------------------------------             
    input_parameters=rtmsim.input_args_assign_parameters(i_interactive,celltype,patchparameters,patchparameters1,patchparameters2,patchparameters3,patchparameters4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,patchids1,patchids2,patchids3,patchids4,inletpatchids,mu_resin_val,N)
    return_parameters=assign_parameters(input_parameters)
    cellthickness=return_parameters.cellthickness
    cellporosity=return_parameters.cellporosity
    cellpermeability=return_parameters.cellpermeability
    cellalpha=return_parameters.cellalpha
    celldirection=return_parameters.celldirection
    cellviscosity=return_parameters.cellviscosity
    celltype=return_parameters.celltype
    
    #--------------------------------------------------------------------------    
    #  Create local cell coordinate systems
    #--------------------------------------------------------------------------    
    input_cs=rtmsim.input_args_create_cs(N, cellgridid, gridx, gridy, gridz, cellcenterx,cellcentery,cellcenterz, faces, cellneighboursarray, celldirection, cellthickness,maxnumberofneighbours)
    return_cs=create_coordinate_systems(input_cs)
    cellvolume=return_cs.cellvolume
    cellcentertocellcenterx=return_cs.cellcentertocellcenterx
    cellcentertocellcentery=return_cs.cellcentertocellcentery
    T11=return_cs.T11
    T12=return_cs.T12
    T21=return_cs.T21
    T22=return_cs.T22
    cellfacenormalx=return_cs.cellfacenormalx
    cellfacenormaly=return_cs.cellfacenormaly
    cellfacearea=return_cs.cellfacearea

    #----------------------------------------------------------------------
    # Initial time step calculation
    #----------------------------------------------------------------------
    area=minimum(cellvolume./cellthickness)
    maxspeed=max(maximum(cellpermeability./cellviscosity),maximum(cellalpha.*cellpermeability./cellviscosity))*(p_a_val-p_init_val)/minimum(cellvolume./cellthickness)  #sqrt(area)
    betat1=1
    deltat=betat1*sqrt(area)/maxspeed
    deltat_initial=deltat


    #----------------------------------------------------------------------
    # Array initialization
    #----------------------------------------------------------------------
    rho_old=Vector{Float64}(undef, N)
    u_old=Vector{Float64}(undef, N)
    v_old=Vector{Float64}(undef, N)
    p_old=Vector{Float64}(undef, N)
    gamma_old=Vector{Float64}(undef, N)
    rho_new=Vector{Float64}(undef, N)
    u_new=Vector{Float64}(undef, N)
    v_new=Vector{Float64}(undef, N)
    p_new=Vector{Float64}(undef, N)
    gamma_new=Vector{Float64}(undef, N)
    gamma_out=Vector{Float64}(undef, N)
    for ind in 1:N
        u_old[ind]=u_init
        v_old[ind]=v_init
        rho_old[ind]=rho_init
        p_old[ind]=p_init
        gamma_old[ind]=gamma_init
    end
    for ind in 1:N
        u_new[ind]=-9e9
        v_new[ind]=-9e9
        rho_new[ind]=-9e9
        p_new[ind]=-9e9
        gamma_new[ind]=-9e9
        gamma_out[ind]=-9e9
    end        
    thickness_factor=Vector{Float64}(undef, N)
    volume_factor=Vector{Float64}(undef, N)
    face_factor=Array{Float64}(undef, N, maxnumberofneighbours)
    porosity_factor=Vector{Float64}(undef, N)
    permeability_factor=Vector{Float64}(undef, N)
    viscosity_factor=Vector{Float64}(undef, N)

    if i_restart==1
        if ~isfile(restartfilename)
            errorstring=string("File ",restartfilename," not existing"  * "\n") 
            error(errorstring)
        end
        @load restartfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out
        u_old=u_new
        v_old=v_new
        rho_old=rho_new
        p_old=p_new
        gamma_old=gamma_new
        t_restart=t
    else 
        t_restart=0
    end

    #----------------------------------------------------------------------
    # Define simulation time and intermediate output times
    #----------------------------------------------------------------------
    t_out=0
    t_progressbar=0
    t=0
    tmin=n_pics*deltat
    tmax=max(tmin,tmax)    

    #----------------------------------------------------------------------
    # Boundary conditions
    #----------------------------------------------------------------------
    for ind in 1:N
        if celltype[ind]==-1  #pressure boundary
            u_old[ind]=u_a
            v_old[ind]=v_a
            rho_old[ind]=rho_a
            p_old[ind]=p_a
            gamma_old[ind]=gamma_a
        elseif celltype[ind]==-2  #pressure outlet
            u_old[ind]=u_init
            v_old[ind]=v_init
            rho_old[ind]=rho_init
            p_old[ind]=p_init
            gamma_old[ind]=gamma_init
        end
    end

    if i_model==2
        #----------------------------------------------------------------------
        # Optional initialization if i_model=2,3,.. 
        #----------------------------------------------------------------------
        # -read text file with parameters
        # -array initialization
        # -boundary conditions
        # e.g. for vacuum infusion: Model for thickness, permeability and porosity change
        #      or temperature equation and degree of cure equation with modified viscosity
    end

    #Abort if no pressure inlet is defined, neither interactively nor as patch
    if i_restart==0
        inds1=findall(isequal(-1),celltype)
        if isempty(inds1)
            errorstring="No pressure inlet ports defined"
            error(errorsting)
        end
    end

    #----------------------------------------------------------------------
    # Time evolution
    #----------------------------------------------------------------------
    n_progressbar=20
    deltat_progressbar=tmax/n_progressbar
    p=Progress(n_progressbar)
    iter=1
    while t<=tmax           
        for ind in 1:N 
            if i_model==1
                thickness_factor[ind]=1.0  #change in cell thickness
                volume_factor[ind]=1.0  #change in cell volume do to cell thickness change
                for i_neighbour in 1:maxnumberofneighbours
                    face_factor[ind,i_neighbour]=1.0  #change is cell boundary area as average of the change in the two neighbouring cells
                end
                porosity_factor[ind]=1.0  #change in porosity
                permeability_factor[ind]=1.0  #change in permeability
                viscosity_factor[ind]=1.0  #change in viscosity
            elseif i_model==2
                #Optional initialization if i_model=2,3,.. 
            end
        end

        for ind in 1:N
            if celltype[ind]==1  || celltype[ind]==-3 
                #Pressure gradient calculation
                input_gradient=rtmsim.input_args_gradient(3,ind,p_old,cellneighboursarray,cellcentertocellcenterx,cellcentertocellcentery)
                return_gradient=numerical_gradient(input_gradient)
                dpdx=return_gradient.dpdx
                dpdy=return_gradient.dpdy
                
                #FV scheme for rho,u,v,vof conservation laws
                cellneighboursline=cellneighboursarray[ind,:]
                cellneighboursline=cellneighboursline[cellneighboursline .> 0]
                len_cellneighboursline=length(cellneighboursline)
                F_rho_num=0.0;F_rho_num_add=0.0
                F_u_num=0.0;F_u_num_add=0.0
                F_v_num=0.0;F_v_num_add=0.0
                F_gamma_num=0.0;F_gamma_num_add=0.0
                F_gamma_num1=0.0;F_gamma_num1_add=0.0
                for i_neighbour=1:len_cellneighboursline
                    i_P=ind
                    i_A=cellneighboursarray[ind,i_neighbour]      
                    rho_P=rho_old[i_P]
                    rho_A=rho_old[i_A]
                    u_P=u_old[i_P]
                    v_P=v_old[i_P]
                    uvec=[T11[ind,i_neighbour] T12[ind,i_neighbour]; T21[ind,i_neighbour] T22[ind,i_neighbour]]*[u_old[i_A];v_old[i_A]]
                    u_A=uvec[1]
                    v_A=uvec[2]
                    gamma_P=gamma_old[i_P]
                    gamma_A=gamma_old[i_A]
                    A=cellfacearea[i_P,i_neighbour]*face_factor[i_P,i_neighbour]
                    n_x=cellfacenormalx[i_P,i_neighbour]
                    n_y=cellfacenormaly[i_P,i_neighbour]
                    vars_P=[rho_P,u_P,v_P,gamma_P]
                    vars_A=[rho_A,u_A,v_A,gamma_A]
                    if i_A>0 && (celltype[i_A]==1 || celltype[i_A]==-3)  #neighbour is inner or wall cell                            
                        meshparameters=[n_x,n_y,A]
                        input_flux=rtmsim.input_args_flux(1,vars_P,vars_A,meshparameters)
                        return_flux=numerical_flux_function(input_flux)
                        F_rho_num_add=return_flux.F_rho_num_add
                        F_u_num_add=return_flux.F_u_num_add
                        F_v_num_add=return_flux.F_v_num_add
                        F_gamma_num_add=return_flux.F_gamma_num_add
                        F_gamma_num1_add=return_flux.F_gamma_num1_add

                        F_rho_num=F_rho_num+F_rho_num_add
                        F_u_num=F_u_num+F_u_num_add
                        F_v_num=F_v_num+F_v_num_add
                        F_gamma_num=F_gamma_num+F_gamma_num_add
                        F_gamma_num1=F_gamma_num1+F_gamma_num1_add  
                    end       
                    if i_A>0 && (celltype[i_A]==-1 || celltype[i_A]==-2)  #neighbour is pressure inlet or outlet
                        A=A*cellthickness[i_P]/(0.5*(cellthickness[i_P]+cellthickness[i_A]))
                        meshparameters=[n_x,n_y,A]
                        if celltype[i_A]==-2  #pressure outlet
                            n_dot_u=dot([n_x; n_y],[u_P; v_P])
                        elseif celltype[i_A]==-1  #pressure inlet
                            n_dot_u=min(0,-1/(cellviscosity[i_P]*viscosity_factor[i_P])*dot([cellpermeability[i_P]*permeability_factor[i_P] 0; 0 cellalpha[i_P]*cellpermeability[i_P]*permeability_factor[ind]]*[dpdx;dpdy],[cellfacenormalx[i_P,i_neighbour];cellfacenormaly[i_P,i_neighbour]]))  #inflow according to Darcy's law and no backflow possible
                        end
                        input_flux_boundary=rtmsim.input_args_flux_boundary(1,vars_P,vars_A,meshparameters,n_dot_u)
                        return_flux_boundary=numerical_flux_function_boundary(input_flux_boundary)
                        F_rho_num_add=return_flux_boundary.F_rho_num_add
                        F_u_num_add=return_flux_boundary.F_u_num_add
                        F_v_num_add=return_flux_boundary.F_v_num_add
                        F_gamma_num_add=return_flux_boundary.F_gamma_num_add
                        F_gamma_num1_add=return_flux_boundary.F_gamma_num1_add
                        F_rho_num=F_rho_num+F_rho_num_add
                        F_u_num=F_u_num+F_u_num_add
                        F_v_num=F_v_num+F_v_num_add
                        F_gamma_num=F_gamma_num+F_gamma_num_add
                        F_gamma_num1=F_gamma_num1+F_gamma_num1_add 
                    end      
                end

                rho_new[ind]=rho_old[ind]-deltat*F_rho_num/(cellvolume[ind]*volume_factor[ind])
                rho_new[ind]=max(rho_new[ind],0.0)
                S_u=-dpdx
                u_new[ind]=(rho_old[ind]*u_old[ind]-deltat*F_u_num/(cellvolume[ind]*volume_factor[ind])+S_u*deltat)/(rho_new[ind]+(cellviscosity[ind]*viscosity_factor[ind])/(cellpermeability[ind]*permeability_factor[ind])*deltat)
                S_v=-dpdy
                v_new[ind]=(rho_old[ind]*v_old[ind]-deltat*F_v_num/(cellvolume[ind]*volume_factor[ind])+S_v*deltat)/(rho_new[ind]+(cellviscosity[ind]*viscosity_factor[ind])/(cellalpha[ind]*cellpermeability[ind]*permeability_factor[ind])*deltat)   
                gamma_new[ind]=((cellporosity[ind]*porosity_factor[ind])*gamma_old[ind]-deltat*(F_gamma_num-gamma_old[ind]*F_gamma_num1)/(cellvolume[ind]*volume_factor[ind]))/(cellporosity[ind]*porosity_factor[ind]) 
                gamma_new[ind]=min(1,gamma_new[ind])
                gamma_new[ind]=max(0,gamma_new[ind])
                #EOS:
                if gamma>1.01
                    p_new[ind]=ap1*rho_new[ind]^2+ap2*rho_new[ind]+ap3
                else
                    p_new[ind]=kappa*rho_new[ind]^gamma
                end
            end
        end

        #boundary conditions, only for pressure boundary conditions
        for ind in 1:N
            if celltype[ind]==-1  #pressure inlet
                u_new[ind]=u_a
                v_new[ind]=v_a
                rho_new[ind]=rho_a
                p_new[ind]=p_a
            elseif celltype[ind]==-2  #pressure outlet
                u_new[ind]=u_init
                v_new[ind]=v_init
                rho_new[ind]=rho_init
                p_new[ind]=p_init
            end
            if celltype[ind]==-1  #pressure inlet
                gamma_new[ind]=gamma_a
            elseif celltype[ind]==-2  #pressure outlet
                gamma_new[ind]=gamma_init
            end
        end 

        #prepare arrays for next time step
        u_old=u_new
        v_old=v_new
        rho_old=rho_new
        p_old=p_new
        gamma_old=gamma_new

        #Progress bar with percentage during run
        prozent = (t/tmax)*100  
        if t>=t_progressbar
            #println(string(string(prozent),"%"))
            t_progressbar=t_progressbar+deltat_progressbar   
            next!(p)
        end      

        #Adaptive time stepping
        if iter>n_pics
            inds1=findall(isequal(1),celltype)
            inds2=findall(isequal(-3),celltype)
            inds=vcat(inds1,inds2)               
            weight_deltatnew=0.5  #0.1  #
            if gamma>=100
                betat2=0.1*0.1
            else
                betat2=0.1
            end
            deltat1=(1-weight_deltatnew)*deltat+weight_deltatnew* betat2*minimum( (sqrt.(cellvolume[inds]./cellthickness[inds])) ./ sqrt.(u_new[inds].^2+v_new[inds].^2) )  
            deltat=deltat1
            #deltat2=(1-weight_deltatnew)*deltat+weight_deltatnew* betat2*minimum( (sqrt.(cellvolume[inds]./cellthickness[inds])) ./ 340)  #sqrt.(gamma*p_new[inds]./rho_new[inds]) )  
            #deltat=min(deltat1,deltat2)  #minimum of convection and wave
            deltatmax=tmax/(4*n_pics) #at least four steps between writing output
            deltat=min(deltat,deltatmax)
        end

        #Save intermediate data
        if t>=t_out  || (t+deltat>tmax)
            if i_model==1
                if i_restart==1
                    t_temp=t
                    t=t+t_restart
                end
                for i in 1:N
                    gamma_out[i]=gamma_new[i] 
                end
                inds=findall(isequal(-1),celltype) #if pressure inlet cells are present
                for i in 1:length(inds)
                    gamma_out[inds[i]]=-1.0 
                end
                inds=findall(isequal(-2),celltype) #if pressure outlet cells are present, they should not be plotted in the gamma-plot because not updated
                for i in 1:length(inds)
                    gamma_out[inds[i]]=-2.0 
                end
                if t>=(tmax+t_restart)-1.5*deltat
                    t_temp1=t
                    t=(tmax+t_restart)
                end
                outputfilename=string("output_", string(n_out), ".jld2")                
                @save outputfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out
                
                #temporary output in Matlab mat-format
                #outputfilename=string("output_", string(n_out), ".mat") 
                #matwrite(outputfilename, Dict("t" => t,"rho_new" => rho_new,"u_new" => u_new,"v_new" => v_new,"p_new" => p_new,"gamma_new" => gamma_new,"gridx" => gridx,"gridy" => gridy,"gridz" => gridz,"cellgridid" => cellgridid,"N" => N,"n_out" => n_out))
                
                outputfilename=string("results.jld2")
                @save outputfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out
                if t>=(tmax+t_restart)-deltat
                    t=t_temp1
                end                    
                if i_restart==1
                    t=t_temp
                end       
            end     
            n_out=n_out+1
            t_out=t_out+tmax/n_pics            
        end
        
        
        if i_model==2
            #----------------------------------------------------------------------
            # Optional time marching etc. for i_model=2,3,.. 
            #----------------------------------------------------------------------
            # -for ind in 1:N loop with calculation of fluxes,...
            # -boundary conditions
            # -array preparation for next time step
            # -write save data to output files
        end

        iter=iter+1
        t=t+deltat 
    end
end


