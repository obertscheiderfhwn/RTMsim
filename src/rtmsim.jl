# RTMsim - A Julia module for filling simulations in Resin Transfer Moulding with the Finite Area Method
# 
# In order to use RTMsim follow the following steps:
# - Download Julia from https://julialang.org/downloads/ and add Julia to path such that can be started from command line.
# - Open Julia terminal, change to package manager with ] and 
#   add Gtk GLMakie Makie NativeFileDialog Glob LinearAlgebra JLD2 GeometryBasics Random FileIO ProgressMeter.
# - One has access to all functions through the Julia terminal. Open a Julia terminal, 
#   change to the directory with the RTMsim repository with cd("path") and 
#   call all functions after include("rtmsim.jl"). 
#   Popular functions are:
#       rtmsim.plot_mesh("..\\meshfiles\\mesh_permeameter1_foursets.bdf",1) for plotting the mesh defined in the bdf-file
#       rtmsim.plot_sets("..\\meshfiles\\mesh_permeameter1_foursets.bdf") for plotting the sets specified in the bdf-file
#       rtmsim.rtmsim_rev1(1,"..\\meshfiles\\mesh_permeameter1_foursets.bdf",200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2,0,"results.jld2",0,0.01,16) for starting a simulation with different patches and race tracking
#       rtmsim.rtmsim_rev1(1,"..\\meshfiles\\mesh_permeameter1_foursets.bdf",200, 101325,1.225,1.4,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,2,2,2,1,"results.jld2",0,0.01,16) for continuing the previous simulation
#       rtmsim.plot_mesh("..\\meshfiles\\mesh_annulusfiller1.bdf",2) for the manual selection of inlet ports
#       rtmsim.rtmsim_rev1(1,"..\\meshfiles\\mesh_annulusfiller1.bdf",200, 0.35e5,1.205,1.4,0.06, 0.35e5,0.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 0,0,0,0, 0,"results.jld2",1,0.01,16) for starting only with the interactively selected inlet ports
#       rtmsim.plot_results("results.jld2") for plotting the final filling and pressure contours
#       rtmsim.plot_overview(-1,-1) for plotting the filling contours at four equidistant time instances
#       rtmsim.plot_filling(-1,-1) for plotting the filling at different time instances selected with a slider bar
#       rtmsim.start_rtmsim("..\\inputfiles\\input.txt") for starting a simulation with the parameters specified in the text file input.txt
#       rtmsim.rtmsim_rev1(1,"..\\meshfiles\\mesh_permeameter1_meshrefinement.bdf",40, 101325,1000,200,0.06, 1.35e5,1.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-11,1,1,0,0, 3e-3,0.7,3e-9,1,1,0,0, 1,0,0,0,0,"results.jld2",0,0.01,16) for nearly incompressible fluid


module rtmsim
    using Glob, LinearAlgebra, JLD2, GeometryBasics, GLMakie, Makie, Random, FileIO, ProgressMeter
    GLMakie.activate!()    
    #using MAT  #temporary output in Matlab mat-format
            
    function start_rtmsim(inputfilename)
        # reads the text input file and calls the solver with the read parameters
        
        print("Read input file "*string(inputfilename)*"\n")
        if ~isfile(inputfilename);
            errorstring=string("File ",inputfilename," not existing"* "\n"); 
            error(errorstring);
        end
        i_model=[]; meshfilename=[]; tmax=[]; 
        p_ref=[]; rho_ref=[]; gamma=[]; mu_resin_val=[]; p_a_val=[]; p_init_val=[]; 
        t_val=[]; porosity_val=[]; K_val=[]; alpha_val=[]; refdir1_val=[]; refdir2_val=[]; refdir3_val=[]; 
        t1_val=[]; porosity1_val=[]; K1_val=[]; alpha1_val=[]; refdir11_val=[]; refdir21_val=[]; refdir31_val=[]; 
        t2_val=[]; porosity2_val=[]; K2_val=[]; alpha2_val=[]; refdir12_val=[]; refdir22_val=[]; refdir32_val=[]; 
        t3_val=[]; porosity3_val=[]; K3_val=[]; alpha3_val=[]; refdir13_val=[]; refdir23_val=[]; refdir33_val=[]; 
        t4_val=[]; porosity4_val=[]; K4_val=[]; alpha4_val=[]; refdir14_val=[]; refdir24_val=[]; refdir34_val=[];
        patchtype1val=[]; patchtype2val=[]; patchtype3val=[]; patchtype4val=[]; 
        i_restart=[]; restartfilename=[]; i_interactive=[]; r_p=[]; n_pics=[];
        open(inputfilename, "r") do fid
            i_line=1;
            while !eof(fid)
                thisline=readline(fid)
                print(string(thisline)*"\n")
                txt1=split(thisline," ")
                if i_line==1;            
                    i_model=parse(Int64,txt1[1]);
                elseif i_line==2;
                    meshfilename=txt1[1];
                elseif i_line==3;
                    tmax=parse(Float64,txt1[1]);
                elseif i_line==4;
                    p_ref=parse(Float64,txt1[1]);
                    rho_ref=parse(Float64,txt1[2]);
                    gamma=parse(Float64,txt1[3]);
                    mu_resin_val=parse(Float64,txt1[4]);
                elseif i_line==5;
                    p_a_val=parse(Float64,txt1[1]);
                    p_init_val=parse(Float64,txt1[2]);
                elseif i_line==6;
                    t_val=parse(Float64,txt1[1]);
                    porosity_val=parse(Float64,txt1[2]);
                    K_val=parse(Float64,txt1[3]);
                    alpha_val=parse(Float64,txt1[4]);
                    refdir1_val=parse(Float64,txt1[5]);
                    refdir2_val=parse(Float64,txt1[6]);
                    refdir3_val=parse(Float64,txt1[7]);
                elseif i_line==7;
                    t1_val=parse(Float64,txt1[1]);
                    porosity1_val=parse(Float64,txt1[2]);
                    K1_val=parse(Float64,txt1[3]);
                    alpha1_val=parse(Float64,txt1[4]);
                    refdir11_val=parse(Float64,txt1[5]);
                    refdir21_val=parse(Float64,txt1[6]);
                    refdir31_val=parse(Float64,txt1[7]);
                elseif i_line==8;
                    t2_val=parse(Float64,txt1[1]);
                    porosity2_val=parse(Float64,txt1[2]);
                    K2_val=parse(Float64,txt1[3]);
                    alpha2_val=parse(Float64,txt1[4]);
                    refdir12_val=parse(Float64,txt1[5]);
                    refdir22_val=parse(Float64,txt1[6]);
                    refdir32_val=parse(Float64,txt1[7]);
                elseif i_line==9;            
                    t3_val=parse(Float64,txt1[1]);
                    porosity3_val=parse(Float64,txt1[2]);
                    K3_val=parse(Float64,txt1[3]);
                    alpha3_val=parse(Float64,txt1[4]);
                    refdir13_val=parse(Float64,txt1[5]);
                    refdir23_val=parse(Float64,txt1[6]);
                    refdir33_val=parse(Float64,txt1[7]);
                elseif i_line==10;
                    t4_val=parse(Float64,txt1[1]);
                    porosity4_val=parse(Float64,txt1[2]);
                    K4_val=parse(Float64,txt1[3]);
                    alpha4_val=parse(Float64,txt1[4]);
                    refdir14_val=parse(Float64,txt1[5]);
                    refdir24_val=parse(Float64,txt1[6]);
                    refdir34_val=parse(Float64,txt1[7]);
                elseif i_line==11;
                    patchtype1val=parse(Int64,txt1[1]);
                    patchtype2val=parse(Int64,txt1[2]);
                    patchtype3val=parse(Int64,txt1[3]);
                    patchtype4val=parse(Int64,txt1[4]);
                elseif i_line==12;
                    i_restart=parse(Int64,txt1[1]);
                    restartfilename=txt1[2];
                elseif i_line==13;
                    i_interactive=parse(Int64,txt1[1]);
                    r_p= parse(Float64,txt1[2]);
                elseif i_line==14;
                    n_pics=parse(Int64,txt1[1]);
                end
                i_line=i_line+1;
                if i_line==15;break;end
            end
        end        
        print(" "*"\n")
        rtmsim_rev1(i_model,meshfilename,tmax,
            p_ref,rho_ref,gamma,mu_resin_val,
            p_a_val,p_init_val,
            t_val,porosity_val,K_val,alpha_val,refdir1_val,refdir2_val,refdir3_val,
            t1_val,porosity1_val,K1_val,alpha1_val,refdir11_val,refdir21_val,refdir31_val,
            t2_val,porosity2_val,K2_val,alpha2_val,refdir12_val,refdir22_val,refdir32_val,
            t3_val,porosity3_val,K3_val,alpha3_val,refdir13_val,refdir23_val,refdir33_val,
            t4_val,porosity4_val,K4_val,alpha4_val,refdir14_val,refdir24_val,refdir34_val,
            patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_restart,restartfilename,i_interactive,r_p,n_pics);
    end

    function rtmsim_rev1(i_model,meshfilename,tmax,
        p_ref,rho_ref,gamma,mu_resin_val,
        p_a_val,p_init_val,
        t_val,porosity_val,K_val,alpha_val,refdir1_val,refdir2_val,refdir3_val,
        t1_val,porosity1_val,K1_val,alpha1_val,refdir11_val,refdir21_val,refdir31_val,
        t2_val,porosity2_val,K2_val,alpha2_val,refdir12_val,refdir22_val,refdir32_val,
        t3_val,porosity3_val,K3_val,alpha3_val,refdir13_val,refdir23_val,refdir33_val,
        t4_val,porosity4_val,K4_val,alpha4_val,refdir14_val,refdir24_val,refdir34_val,
        patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_restart,restartfilename,i_interactive,r_p,n_pics);
        # RTMsim solver with main steps
        # - Simulation initialization
        # - Read mesh file and prepare patches  
        # - Find neighbouring cells
        # - Assign parameters to cells
        # - Create local cell coordinate systems
        # - Calculate initial time step
        # - Array initialization
        # - Define simulation time and intermediate output times
        # - Boundary conditions
        # - (Optional initialization if i_model=2,3,..)
        # - Time evolution (for loops over all indices inside a while loop for time evolution)
        #     - Calculation of correction factors for cell thickness, porosity, permeability, viscosity
        #     - Pressure gradient calculation
        #     - Numerical flux function calculation
        #     - Update of rho, u, v, gamma and p according to conservation laws and equation of state
        #     - Boundary conditions
        #     - Prepare arrays for next time step
        #     - Saving of intermediate data
        #     - (Opional time marching etc. for i_model=2,3,...)
        #     - Calculation of adaptive time step 
        #


        #----------------------------------------------------------------------
        # Simulation initialization
        #----------------------------------------------------------------------
        
        # Well defined variable types, except for strings meshfilename,restartfilename
        tmax=Float64(tmax);
        p_ref=Float64(p_ref);rho_ref=Float64(rho_ref);gamma=Float64(gamma);mu_resin_val=Float64(mu_resin_val);
        p_a_val=Float64(p_a_val);p_init_val=Float64(p_init_val);
        t_val=Float64(t_val);porosity_val=Float64(porosity_val);K_val=Float64(K_val);alpha_val=Float64(alpha_val);refdir1_val=Float64(refdir1_val);refdir2_val=Float64(refdir2_val);refdir3_val=Float64(refdir3_val);
        t1_val=Float64(t1_val);porosity1_val=Float64(porosity1_val);K1_val=Float64(K1_val);alpha1_val=Float64(alpha1_val);refdir11_val=Float64(refdir11_val);refdir21_val=Float64(refdir21_val);refdir31_val=Float64(refdir31_val);
        t2_val=Float64(t2_val);porosity2_val=Float64(porosity2_val);K2_val=Float64(K2_val);alpha2_val=Float64(alpha2_val);refdir12_val=Float64(refdir12_val);refdir22_val=Float64(refdir22_val);refdir32_val=Float64(refdir32_val);
        t3_val=Float64(t3_val);porosity3_val=Float64(porosity3_val);K3_val=Float64(K3_val);alpha3_val=Float64(alpha3_val);refdir13_val=Float64(refdir13_val);refdir23_val=Float64(refdir23_val);refdir33_val=Float64(refdir33_val);
        t4_val=Float64(t4_val);porosity4_val=Float64(porosity4_val);K4_val=Float64(K4_val);alpha4_val=Float64(alpha4_val);refdir14_val=Float64(refdir14_val);refdir24_val=Float64(refdir24_val);refdir34_val=Float64(refdir34_val);
        patchtype1val=Int64(patchtype1val);patchtype2val=Int64(patchtype2val);patchtype3val=Int64(patchtype3val);patchtype4val=Int64(patchtype4val);
        i_restart=Int64(i_restart);i_interactive=Int64(i_interactive);n_pics=Int64(n_pics);
 
        #License statement
        print("\n")
        print("RTMsim version 0.2 \n");
        print("RTMsim is Julia code with GUI which simulates the mold filling in Liquid Composite Molding (LCM) manufacturing process. \n");
        print("Copyright (C) 2022 Christof Obertscheider / University of Applied Sciences Wiener Neustadt (FHWN)  \n");    
        print("\n")
        print("This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. \n");
        print("\n")
        print("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \n");
        print("\n")
        print("You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.\n");
        print("\n")
        print("This software is free of charge and may be used for commercial and academic purposes.  Please mention the use of this software at an appropriate place in your work. \n")
        print("\n")
        print("Submit bug reports to christof.obertscheider@fhwn.ac.at \n");
        print("\n")

        #Output simulation parameter overview
        print("\n")
        print("RTMsim started with the following parameters:\n")
        print("i_model=",i_model,"\n")
            if i_model!=1;
                errorstring=string("Only iso-thermal RTM implemented, i.e. i_model must be =1 instead of =",string(i_model)*"\n"); 
                error(errorstring);
            end
        print("meshfilename=",meshfilename,"\n")
            if ~isfile(meshfilename);
                errorstring=string("File ",meshfilename," not existing"* "\n"); 
                error(errorstring);
            end
        print("tmax=",string(tmax),"\n")
            if tmax<=0.0;
                errorstring="tmax must be greater than zero"
                error(errorstring)
            end
            #Limit number of results time steps between n_pics_min and n_pics_max and make it multiple of 4
            n_pics_input=n_pics;
            n_pics_min=Int64(4);
            n_pics_max=Int64(100);
            if mod(n_pics,4)!=0;
                n_pics=round(n_pics/4)*4;
            end
            if n_pics<n_pics_min
                n_pics=n_pics_min;
            end
            if n_pics>n_pics_max
                n_pics=n_pics_max; 
            end       
            if n_pics>n_pics_max;
                n_pics=n_pics_max;
            end
            n_pics=Int64(n_pics);
        if n_pics_input!=n_pics;
            print("n_pics changed to n_pics=",string(n_pics),"\n")
        else
            print("n_pics=",string(n_pics),"\n")
        end
        print("i_interactive=",string(i_interactive),"\n");   
            if i_interactive!=0 && i_interactive!=1 && i_interactive!=2;
                errorstring="Wrong value for i_interactive (must be=0,1,2)"
                error(errorstring)
            end 
        if i_restart==1; 
            print("i_restart,restartfilename=",string(i_restart), ",", restartfilename,"\n"); 
            if i_restart!=0 && i_restart!=1 ;
                errorstring="Wrong value for i_restart (must be=0,1)"
                error(errorstring)
            end 
            if ~isfile(restartfilename);
                errorstring=string("File ",restartfilename," not existing"* "\n"); 
                error(errorstring);
            end
        end
        print("p_ref,rho_ref,gamma,mu=", string(p_ref), ",", string(rho_ref), ",", string(gamma), ",", string(mu_resin_val),"\n")
        if p_ref<=0.0 || rho_ref<=0.0 || gamma<1.0 || mu_resin_val<=0.0;
            errorstring="Wrong value for p_ref,rho_ref,gamma,mu (must be >0.0,>0.0,>1.0,>0.0)"
            error(errorstring)
        end 
        print("p_a_val,p_init_val=", string(p_a_val), ",", string(p_init_val),"\n")
            if p_a_val<=p_init_val;
                errorstring="Injection pressure must be greater than initial pressure"
                error(errorstring)
            end
            if p_a_val<=0.0 || p_init_val<0.0;
                errorstring="Wrong value for p_a_val,p_init_val (must be >0.0,>0.0)"
                error(errorstring)
            end 

        #Maximum number of cell neighbours
        maxnumberofneighbours=10;

        #Delete old files and abort if meshfile is not existing
        if ~isfile(meshfilename);
            errorstring=string("File ",meshfilename," not existing"* "\n"); 
            error(errorstring);
        end
        if i_restart==1;
            cp(restartfilename,"restart.jdl2";force=true);
        end
        delete_files();     
        if i_restart==1;
            cp("restart.jdl2",restartfilename;force=true);
        end
        n_out=Int64(0);

        #Assign and prepare physical parameters
        refdir_val=[refdir1_val,refdir2_val,refdir3_val]  #Vector
        u_a=Float64(0.0);  
        u_b=Float64(0.0); 
        u_init=Float64(0.0); 
        v_a=Float64(0.0); 
        v_b=Float64(0.0); 
        v_init=Float64(0.0); 
        p_a=p_a_val;
        p_init=p_init_val;
        p_b=p_a_val;
        #Normalization for Delta p: p->p-p_init
            p_eps=Float64(0.001e5); #Float64(0.000e5);  #
            p_a=p_a-p_init+p_eps;
            p_init=p_init-p_init+p_eps;
            p_b=p_a-p_init+p_eps;
            p_ref=p_ref;  #p_ref-p_init+p_eps;
        kappa=p_ref/(rho_ref^gamma);
        #Lookuptable for adiabatic law (required for stability)
            p_int1=Float64(0.0e5); rho_int1=(p_int1/kappa)^(1/gamma);
            p_int2=Float64(0.1e5); rho_int2=(p_int2/kappa)^(1/gamma);
            p_int3=Float64(0.5e5); rho_int3=(p_int3/kappa)^(1/gamma);
            p_int4=Float64(1.0e5); rho_int4=(p_int4/kappa)^(1/gamma);
            p_int5=Float64(10.0e5); rho_int5=(p_int5/kappa)^(1/gamma);
            p_int6=Float64(100.0e5); rho_int6=(p_int6/kappa)^(1/gamma);
            A=[rho_int1^2 rho_int1 Float64(1.0); rho_int3^2 rho_int3 Float64(1.0); rho_int4^2 rho_int4 Float64(1.0)];
            b=[p_int1;p_int3;p_int4];
            apvals=A\b;
            ap1=apvals[1];ap2=apvals[2];ap3=apvals[3];
        rho_a=(p_a/kappa)^(Float64(1)/gamma);
        rho_b=(p_b/kappa)^(Float64(1)/gamma);
        rho_init=(p_init/kappa)^(Float64(1)/gamma);

        if gamma>=100;  #insert here coefficients for an incompressible EOS with resin mass density as rho_ref 
                        #at p_b and 0.9*rho_ref at p_a but the EOS is for deltap and consequently normalized pressure values
            #ap1=0;
            #ap2=(p_b-p_init)/(0.1*rho_ref);
            #ap3=p_b-(p_b-p_init)/0.1; 
            #rho_a=p_a/ap2-ap3/ap2;
            #rho_b=p_b/ap2-ap3/ap2;
            #rho_init=p_init/ap2-ap3/ap2;

            rho_a=rho_ref;
            rho_b=rho_a;
            rho_init=0.0;
            p_int1=p_init; rho_int1=rho_init;
            p_int2=p_init+0.9*(p_a-p_init); rho_int2=0.1*rho_a;
            p_int3=p_a; rho_int3=rho_a;
            #A=[rho_int1^2 rho_int1 Float64(1.0); rho_int2^2 rho_int2 Float64(1.0); rho_int3^2 rho_int3 Float64(1.0)];
            #b=[p_int1;p_int2;p_int3];
            A=[rho_int1^2 rho_int1 Float64(1.0); rho_int3^2 rho_int3 Float64(1.0); 2*rho_int3 Float64(1.0) 0];
            b=[p_int1;p_int3;Float64(0.0)];
            apvals=A\b;
            ap1=apvals[1];ap2=apvals[2];ap3=apvals[3];
            #ap1*rho_new[ind]^2+ap2*rho_new[ind]+ap3;


            print(string("rho_int1: ",string(rho_int1) , "\n" ) );
            print(string("rho_int2: ",string(rho_int2) , "\n" ) );
            print(string("rho_int3: ",string(rho_int3) , "\n" ) );
            print(string("p_int1: ",string(p_int1) , "\n" ) );
            print(string("p_int2: ",string(p_int2) , "\n" ) );
            print(string("p_int3: ",string(p_int3) , "\n" ) );
            
            print(string("ap1: ",string(ap1) , "\n" ) );
            print(string("ap2: ",string(ap2) , "\n" ) ); 
            print(string("ap3: ",string(ap3) , "\n" ) ); 
            print(string("p_a: ",string(p_a) , "\n" ) );
            print(string("p_b: ",string(p_b) , "\n" ) );
            print(string("p_init: ",string(p_init) , "\n" ) );
            print(string("rho_a: ",string(rho_a) , "\n" ) );
            print(string("rho_b: ",string(rho_b) , "\n" ) );
            print(string("rho_init: ",string(rho_init) , "\n" ) );
        end
        
        T_a=Float64(295);
        T_b=Float64(295);
        T_init=Float64(295);
        gamma_a=Float64(1.0);
        gamma_b=Float64(1.0);    
        gamma_init=Float64(0.0);
        paramset=[porosity_val,t_val,K_val,alpha_val,refdir1_val,refdir2_val,refdir3_val];
        paramset1=[porosity1_val,t1_val,K1_val,alpha1_val,refdir11_val,refdir21_val,refdir31_val];
        paramset2=[porosity2_val,t2_val,K2_val,alpha2_val,refdir12_val,refdir22_val,refdir32_val];
        paramset3=[porosity3_val,t3_val,K3_val,alpha3_val,refdir13_val,refdir23_val,refdir33_val];
        paramset4=[porosity4_val,t4_val,K4_val,alpha4_val,refdir14_val,refdir24_val,refdir34_val];


        #--------------------------------------------------------------------------
        # Read mesh file and prepare patches     
        #--------------------------------------------------------------------------
        N,cellgridid,gridx,gridy,gridz,cellcenterx,cellcentery,cellcenterz,patchparameters,patchparameters1,patchparameters2,patchparameters3,patchparameters4,patchids1,patchids2,patchids3,patchids4,inletpatchids=
            read_mesh(meshfilename,paramset,paramset1,paramset2,paramset3,paramset4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_interactive,r_p);

        print(string("parameters for main preform: ",string(patchparameters) , "\n" ) );
        if patchparameters[1]<=0.0 || patchparameters[1]>1.0 || patchparameters[2]<=0.0 || patchparameters[3]<=0 || patchparameters[4]<=0;
            errorstring="Wrong value for porosity,thickness,permeability,alpha (must be between >0 and <=1,>0.0,>0.0,>0.0)"
            error(errorstring)
        end 
        if ~isempty(patchids1); 
            if patchtype1val==1;        
                print("patch 1 is pressure inlet \n"); 
            elseif patchtype1val==2;        
                print(string("parameters for patch 1: ",string(patchparameters1), "\n" ) );
                if patchparameters1[1]<=0.0 || patchparameters1[1]>1.0 || patchparameters1[2]<=0.0 || patchparameters1[3]<=0 || patchparameters1[4]<=0;
                    errorstring="Wrong value for porosity,thickness,permeability,alpha (must be between >0 and <=1,>0.0,>0.0,>0.0)"
                    error(errorstring)
                end 
            elseif patchtype1val==3;        
                print("patch 1 is pressure outlet \n"); 
            end
        end
        if ~isempty(patchids2); 
            if patchtype2val==1;        
                print("patch 4 is pressure inlet \n"); 
            elseif patchtype2val==2;        
                print(string("parameters for patch 2: ",string(patchparameters2), "\n" ) );
                if patchparameters2[1]<=0.0 || patchparameters2[1]>1.0 || patchparameters2[2]<=0.0 || patchparameters2[3]<=0 || patchparameters2[4]<=0;
                    errorstring="Wrong value for porosity,thickness,permeability,alpha (must be between >0 and <=1,>0.0,>0.0,>0.0)"
                    error(errorstring)
                end 
            elseif patchtype2val==3;        
                print("patch 2 is pressure outlet \n"); 
            end
        end
        if ~isempty(patchids3); 
            if patchtype3val==1;        
                print("patch 3 is pressure inlet \n"); 
            elseif patchtype3val==2;        
                print(string("parameters for patch 3: ",string(patchparameters3), "\n" ) );
                if patchparameters3[1]<=0.0 || patchparameters3[1]>1.0 || patchparameters3[2]<=0.0 || patchparameters3[3]<=0 || patchparameters3[4]<=0;
                    errorstring="Wrong value for porosity,thickness,permeability,alpha (must be between >0 and <=1,>0.0,>0.0,>0.0)"
                    error(errorstring)
                end 
            elseif patchtype3val==3;        
                print("patch 3 is pressure outlet \n"); 
            end
        end
        if ~isempty(patchids4); 
            if patchtype4val==1;        
                print("patch 4 is pressure inlet \n"); 
            elseif patchtype4val==2;        
                print(string("parameters for patch 4: ",string(patchparameters4),"\n" ) );
                if patchparameters4[1]<=0.0 || patchparameters4[1]>1.0 || patchparameters4[2]<=0.0 || patchparameters4[3]<=0 || patchparameters4[4]<=0;
                    errorstring="Wrong value for porosity,thickness,permeability,alpha (must be between >0 and <=1,>0.0,>0.0,>0.0)"
                    error(errorstring)
                end 
            elseif patchtype4val==3;        
                print("patch 4 is pressure outlet \n");   
            end
        end
        if patchtype1val!=1 && patchtype2val!=1 && patchtype3val!=1 && patchtype3val!=1 && i_interactive==0 && i_restart==0
            errorstring=string("No inlet defined" * "\n"); 
            error(errorstring);
        end
        if i_interactive==1 || i_interactive==2;
            print("additional inlet defined interactively \n");   
        end

        #--------------------------------------------------------------------------
        #  Find neighbouring cells
        #--------------------------------------------------------------------------    
        faces,cellneighboursarray,celltype = 
            create_faces(cellgridid, N, maxnumberofneighbours);


        #--------------------------------------------------------------------------
        #  Assign parameters to cells
        #--------------------------------------------------------------------------             
        cellthickness, cellporosity, cellpermeability, cellalpha, celldirection, cellviscosity, celltype = 
            assign_parameters(i_interactive,celltype,patchparameters,patchparameters1,patchparameters2,patchparameters3,patchparameters4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,patchids1,patchids2,patchids3,patchids4,inletpatchids,mu_resin_val,N);


        #--------------------------------------------------------------------------    
        #  Create local cell coordinate systems
        #--------------------------------------------------------------------------    
        cellvolume, cellcentertocellcenterx, cellcentertocellcentery, T11, T12, T21, T22, cellfacenormalx, cellfacenormaly, cellfacearea = 
            create_coordinate_systems(N, cellgridid, gridx, gridy, gridz, cellcenterx,cellcentery,cellcenterz, faces, cellneighboursarray, celldirection, cellthickness,maxnumberofneighbours);


        #----------------------------------------------------------------------
        # Initial time step calculation
        #----------------------------------------------------------------------
        area=minimum(cellvolume./cellthickness);
        maxspeed=max(maximum(cellpermeability./cellviscosity),maximum(cellalpha.*cellpermeability./cellviscosity))*(p_a_val-p_init_val)/minimum(cellvolume./cellthickness);  #sqrt(area);
        betat1=1;
        deltat=betat1*sqrt(area)/maxspeed;
        deltat_initial=deltat;


        #----------------------------------------------------------------------
        # Array initialization
        #----------------------------------------------------------------------
        rho_old=Vector{Float64}(undef, N);
        u_old=Vector{Float64}(undef, N);
        v_old=Vector{Float64}(undef, N);
        p_old=Vector{Float64}(undef, N);
        gamma_old=Vector{Float64}(undef, N);
        rho_new=Vector{Float64}(undef, N);
        u_new=Vector{Float64}(undef, N);
        v_new=Vector{Float64}(undef, N);
        p_new=Vector{Float64}(undef, N);
        gamma_new=Vector{Float64}(undef, N);
        gamma_out=Vector{Float64}(undef, N);
        for ind in 1:N;
            u_old[ind]=u_init;
            v_old[ind]=v_init;
            rho_old[ind]=rho_init;
            p_old[ind]=p_init;
            gamma_old[ind]=gamma_init;
        end
        for ind in 1:N
            u_new[ind]=-9e9;
            v_new[ind]=-9e9;
            rho_new[ind]=-9e9;
            p_new[ind]=-9e9;
            gamma_new[ind]=-9e9;
            gamma_out[ind]=-9e9;
        end        
        thickness_factor=Vector{Float64}(undef, N);
        volume_factor=Vector{Float64}(undef, N);
        face_factor=Array{Float64}(undef, N, maxnumberofneighbours);
        porosity_factor=Vector{Float64}(undef, N);
        permeability_factor=Vector{Float64}(undef, N);
        viscosity_factor=Vector{Float64}(undef, N);

        if i_restart==1;
            if ~isfile(restartfilename);
                errorstring=string("File ",restartfilename," not existing"  * "\n"); 
                error(errorstring);
            end
            @load restartfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out
            u_old=u_new;
            v_old=v_new;
            rho_old=rho_new;
            p_old=p_new;
            gamma_old=gamma_new;
            t_restart=t;
        else 
            t_restart=0;
        end

        #----------------------------------------------------------------------
        # Define simulation time and intermediate output times
        #----------------------------------------------------------------------
        t_out=0;
        t_progressbar=0;
        t=0;
        tmin=n_pics*deltat;
        tmax=max(tmin,tmax);    

        #----------------------------------------------------------------------
        # Boundary conditions
        #----------------------------------------------------------------------
        for ind in 1:N;
            if celltype[ind]==-1;  #pressure boundary
                u_old[ind]=u_a;
                v_old[ind]=v_a;
                rho_old[ind]=rho_a;
                p_old[ind]=p_a;
                gamma_old[ind]=gamma_a;
            elseif celltype[ind]==-2;  #pressure outlet
                u_old[ind]=u_init;
                v_old[ind]=v_init;
                rho_old[ind]=rho_init;
                p_old[ind]=p_init;
                gamma_old[ind]=gamma_init;
            end
        end

        if i_model==2;
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
        if i_restart==0;
            inds1=findall(isequal(-1),celltype);
            if isempty(inds1)
                errorstring="No pressure inlet ports defined";
                error(errorsting);
            end
        end

        #----------------------------------------------------------------------
        # Time evolution
        #----------------------------------------------------------------------
        n_progressbar=20;
        deltat_progressbar=tmax/n_progressbar;
        p=Progress(n_progressbar);
        iter=1;
        while t<=tmax;           
            for ind in 1:N 
                if i_model==1;
                    thickness_factor[ind]=Float64(1.0);  #change in cell thickness
                    volume_factor[ind]=Float64(1.0);  #change in cell volume do to cell thickness change
                    for i_neighbour in 1:maxnumberofneighbours
                        face_factor[ind,i_neighbour]=Float64(1.0);  #change is cell boundary area as average of the change in the two neighbouring cells
                    end
                    porosity_factor[ind]=Float64(1.0);  #change in porosity
                    permeability_factor[ind]=Float64(1.0);  #change in permeability
                    viscosity_factor[ind]=Float64(1.0);  #change in viscosity
                elseif i_model==2;
                    #Optional initialization if i_model=2,3,.. 
                end
            end

            for ind in 1:N
                if celltype[ind]==1  || celltype[ind]==-3; 
                    #Pressure gradient calculation
                    #dpdx,dpdy=numerical_gradient(1,ind,p_old,cellneighboursarray,cellcentertocellcenterx,cellcentertocellcentery);
                    dpdx,dpdy=numerical_gradient(3,ind,p_old,cellneighboursarray,cellcentertocellcenterx,cellcentertocellcentery);
                    
                    #FV scheme for rho,u,v,vof conservation laws
                    cellneighboursline=cellneighboursarray[ind,:];
                    cellneighboursline=cellneighboursline[cellneighboursline .> 0]
                    len_cellneighboursline=length(cellneighboursline)
                    F_rho_num=Float64(0.0);F_rho_num_add=Float64(0.0);
                    F_u_num=Float64(0.0);F_u_num_add=Float64(0.0);
                    F_v_num=Float64(0.0);F_v_num_add=Float64(0.0);
                    F_gamma_num=Float64(0.0);F_gamma_num_add=Float64(0.0);
                    F_gamma_num1=Float64(0.0);F_gamma_num1_add=Float64(0.0);
                    for i_neighbour=1:len_cellneighboursline;
                        i_P=ind;
                        i_A=cellneighboursarray[ind,i_neighbour];      
                        rho_P=rho_old[i_P];
                        rho_A=rho_old[i_A];
                        u_P=u_old[i_P];
                        v_P=v_old[i_P];
                        uvec=[T11[ind,i_neighbour] T12[ind,i_neighbour]; T21[ind,i_neighbour] T22[ind,i_neighbour]]*[u_old[i_A];v_old[i_A]];
                        u_A=uvec[1];
                        v_A=uvec[2];
                        gamma_P=gamma_old[i_P];
                        gamma_A=gamma_old[i_A];
                        A=cellfacearea[i_P,i_neighbour]*face_factor[i_P,i_neighbour];
                        n_x=cellfacenormalx[i_P,i_neighbour];
                        n_y=cellfacenormaly[i_P,i_neighbour];
                        vars_P=[rho_P,u_P,v_P,gamma_P];
                        vars_A=[rho_A,u_A,v_A,gamma_A];
                        if i_A>0 && (celltype[i_A]==1 || celltype[i_A]==-3);  #neighbour is inner or wall cell                            
                            meshparameters=[n_x,n_y,A];
                            F_rho_num_add,F_u_num_add,F_v_num_add,F_gamma_num_add,F_gamma_num1_add=numerical_flux_function(1,vars_P,vars_A,meshparameters);
                            F_rho_num=F_rho_num+F_rho_num_add;
                            F_u_num=F_u_num+F_u_num_add;
                            F_v_num=F_v_num+F_v_num_add;
                            F_gamma_num=F_gamma_num+F_gamma_num_add;
                            F_gamma_num1=F_gamma_num1+F_gamma_num1_add;  
                        end       
                        if i_A>0 && (celltype[i_A]==-1 || celltype[i_A]==-2);  #neighbour is pressure inlet or outlet
                            A=A*cellthickness[i_P]/(0.5*(cellthickness[i_P]+cellthickness[i_A]));
                            meshparameters=[n_x,n_y,A];
                            if celltype[i_A]==-2;  #pressure outlet
                                n_dot_u=dot([n_x; n_y],[u_P; v_P]);
                            elseif celltype[i_A]==-1;  #pressure inlet
                                n_dot_u=min(0,-1/(cellviscosity[i_P]*viscosity_factor[i_P])*dot([cellpermeability[i_P]*permeability_factor[i_P] 0; 0 cellalpha[i_P]*cellpermeability[i_P]*permeability_factor[ind]]*[dpdx;dpdy],[cellfacenormalx[i_P,i_neighbour];cellfacenormaly[i_P,i_neighbour]]));  #inflow according to Darcy's law and no backflow possible
                            end
                            F_rho_num_add,F_u_num_add,F_v_num_add,F_gamma_num_add,F_gamma_num1_add=numerical_flux_function_boundary(1,vars_P,vars_A,meshparameters,n_dot_u);
                            F_rho_num=F_rho_num+F_rho_num_add;
                            F_u_num=F_u_num+F_u_num_add;
                            F_v_num=F_v_num+F_v_num_add;
                            F_gamma_num=F_gamma_num+F_gamma_num_add;
                            F_gamma_num1=F_gamma_num1+F_gamma_num1_add; 
                        end      
                    end

                    rho_new[ind]=rho_old[ind]-deltat*F_rho_num/(cellvolume[ind]*volume_factor[ind]);
                    rho_new[ind]=max(rho_new[ind],0.0)
                    S_u=-dpdx;
                    u_new[ind]=(rho_old[ind]*u_old[ind]-deltat*F_u_num/(cellvolume[ind]*volume_factor[ind])+S_u*deltat)/(rho_new[ind]+(cellviscosity[ind]*viscosity_factor[ind])/(cellpermeability[ind]*permeability_factor[ind])*deltat);
                    S_v=-dpdy;
                    v_new[ind]=(rho_old[ind]*v_old[ind]-deltat*F_v_num/(cellvolume[ind]*volume_factor[ind])+S_v*deltat)/(rho_new[ind]+(cellviscosity[ind]*viscosity_factor[ind])/(cellalpha[ind]*cellpermeability[ind]*permeability_factor[ind])*deltat);   
                    gamma_new[ind]=((cellporosity[ind]*porosity_factor[ind])*gamma_old[ind]-deltat*(F_gamma_num-gamma_old[ind]*F_gamma_num1)/(cellvolume[ind]*volume_factor[ind]))/(cellporosity[ind]*porosity_factor[ind]); 
                    gamma_new[ind]=min(1,gamma_new[ind]);
                    gamma_new[ind]=max(0,gamma_new[ind]);
                    #EOS:
                    if gamma>1.01;
                       p_new[ind]=ap1*rho_new[ind]^2+ap2*rho_new[ind]+ap3;
                    else
                        p_new[ind]=kappa*rho_new[ind]^gamma;
                    end
                end
            end

            #boundary conditions, only for pressure boundary conditions
            for ind in 1:N;
                if celltype[ind]==-1;  #pressure inlet
                    u_new[ind]=u_a;
                    v_new[ind]=v_a;
                    rho_new[ind]=rho_a;
                    p_new[ind]=p_a;
                elseif celltype[ind]==-2;  #pressure outlet
                    u_new[ind]=u_init;
                    v_new[ind]=v_init;
                    rho_new[ind]=rho_init;
                    p_new[ind]=p_init;
                end
                if celltype[ind]==-1;  #pressure inlet
                    gamma_new[ind]=gamma_a;
                elseif celltype[ind]==-2;  #pressure outlet
                    gamma_new[ind]=gamma_init;
                end
            end 

            #prepare arrays for next time step
            u_old=u_new;
            v_old=v_new;
            rho_old=rho_new;
            p_old=p_new;
            gamma_old=gamma_new;

            #Progress bar with percentage during run
            prozent = (t/tmax)*100;  
            if t>=t_progressbar;
                #print(string(string(prozent),"%","\n"))
                t_progressbar=t_progressbar+deltat_progressbar;   
                next!(p);
            end      

            #Adaptive time stepping
            if iter>n_pics;
                inds1=findall(isequal(1),celltype);
                inds2=findall(isequal(-3),celltype);
                inds=vcat(inds1,inds2);               
                weight_deltatnew=Float64(0.5);  #0.1;  #
                if gamma>=100;
                    betat2=Float64(0.1*0.1);
                else
                    betat2=Float64(0.1);
                end
                deltat1=(1-weight_deltatnew)*deltat+weight_deltatnew* betat2*minimum( (sqrt.(cellvolume[inds]./cellthickness[inds])) ./ sqrt.(u_new[inds].^2+v_new[inds].^2) );  
                deltat=deltat1
                #deltat2=(1-weight_deltatnew)*deltat+weight_deltatnew* betat2*minimum( (sqrt.(cellvolume[inds]./cellthickness[inds])) ./ 340)  #sqrt.(gamma*p_new[inds]./rho_new[inds]) );  
                #deltat=min(deltat1,deltat2)  #minimum of convection and wave
                deltatmax=tmax/(4*n_pics); #at least four steps between writing output
                deltat=min(deltat,deltatmax);
            end

            #Save intermediate data
            if t>=t_out  || (t+deltat>tmax);
                if i_model==1;
                    if i_restart==1;
                        t_temp=t;
                        t=t+t_restart;
                    end
                    for i in 1:N
                        gamma_out[i]=gamma_new[i] 
                    end
                    inds=findall(isequal(-1),celltype); #if pressure inlet cells are present
                    for i in 1:length(inds)
                        gamma_out[inds[i]]=Float64(-1.0); 
                    end
                    inds=findall(isequal(-2),celltype); #if pressure outlet cells are present, they should not be plotted in the gamma-plot because not updated
                    for i in 1:length(inds)
                        gamma_out[inds[i]]=Float64(-2.0); 
                    end
                    if t>=(tmax+t_restart)-1.5*deltat;
                        t_temp1=t;
                        t=(tmax+t_restart);
                    end
                    outputfilename=string("output_", string(n_out), ".jld2")                
                    @save outputfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out
                    
                    #temporary output in Matlab mat-format
                    #outputfilename=string("output_", string(n_out), ".mat") 
                    #matwrite(outputfilename, Dict("t" => t,"rho_new" => rho_new,"u_new" => u_new,"v_new" => v_new,"p_new" => p_new,"gamma_new" => gamma_new,"gridx" => gridx,"gridy" => gridy,"gridz" => gridz,"cellgridid" => cellgridid,"N" => N,"n_out" => n_out))
                    
                    outputfilename=string("results.jld2")
                    @save outputfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out
                    if t>=(tmax+t_restart)-deltat;
                        t=t_temp1;
                    end                    
                    if i_restart==1;
                        t=t_temp;
                    end       
                end     
                n_out=n_out+1;
                t_out=t_out+tmax/n_pics;            
            end
            
            
            if i_model==2;
                #----------------------------------------------------------------------
                # Optional time marching etc. for i_model=2,3,.. 
                #----------------------------------------------------------------------
                # -for ind in 1:N loop with calculation of fluxes,...
                # -boundary conditions
                # -array preparation for next time step
                # -write save data to output files
            end

            iter=iter+1;
            t=t+deltat; 
        end
    end



    function numerical_gradient(i_method,ind,p_old,cellneighboursarray,cellcentertocellcenterx,cellcentertocellcentery);
        if i_method==1;
            #least square solution to determine gradient
            cellneighboursline=cellneighboursarray[ind,:];
            cellneighboursline=cellneighboursline[cellneighboursline .> 0]
            len_cellneighboursline=length(cellneighboursline)
            bvec=Vector{Float64}(undef,len_cellneighboursline);
            Amat=Array{Float64}(undef,len_cellneighboursline,2);  
            for i_neighbour in 1:len_cellneighboursline;
                i_P=ind;
                i_A=cellneighboursarray[ind,i_neighbour];  
                Amat[i_neighbour,1]=cellcentertocellcenterx[ind,i_neighbour]
                Amat[i_neighbour,2]=cellcentertocellcentery[ind,i_neighbour]
                bvec[i_neighbour]=p_old[i_A]-p_old[i_P];
            end

            if len_cellneighboursline>1;
                xvec=Amat[1:len_cellneighboursline,:]\bvec[1:len_cellneighboursline];
                dpdx=xvec[1];
                dpdy=xvec[2];        
            else
                dpdx=0;
                dpdy=0;
            end
        elseif i_method==2;
            #least square solution to determine gradient with limiter
            cellneighboursline=cellneighboursarray[ind,:];
            cellneighboursline=cellneighboursline[cellneighboursline .> 0]
            len_cellneighboursline=length(cellneighboursline)
            bvec=Vector{Float64}(undef,len_cellneighboursline);
            Amat=Array{Float64}(undef,len_cellneighboursline,2);  
            wi=Vector{Float64}(undef,len_cellneighboursline);
            for i_neighbour in 1:len_cellneighboursline;
                i_P=ind;
                i_A=cellneighboursarray[ind,i_neighbour];  
                exp_limiter=2;
                wi[i_neighbour]=1/(sqrt((cellcentertocellcenterx[ind,i_neighbour])^2+(cellcentertocellcentery[ind,i_neighbour])^2))^exp_limiter;
                Amat[i_neighbour,1]=wi[i_neighbour]*cellcentertocellcenterx[ind,i_neighbour]
                Amat[i_neighbour,2]=wi[i_neighbour]*cellcentertocellcentery[ind,i_neighbour]
                bvec[i_neighbour]=wi[i_neighbour]*(p_old[i_A]-p_old[i_P]);
            end

            if len_cellneighboursline>1
                xvec=Amat[1:len_cellneighboursline,:]\bvec[1:len_cellneighboursline];
                dpdx=xvec[1];
                dpdy=xvec[2];            
            else
                dpdx=0;
                dpdy=0;
            end
        elseif i_method==3;
            #least square solution to determine gradient - runtime optimized
            cellneighboursline=cellneighboursarray[ind,:];
            cellneighboursline=cellneighboursline[cellneighboursline .> 0]
            len_cellneighboursline=length(cellneighboursline)
            bvec=Vector{Float64}(undef,len_cellneighboursline);
            Amat=Array{Float64}(undef,len_cellneighboursline,2);  
            for i_neighbour in 1:len_cellneighboursline;
                i_P=ind;
                i_A=cellneighboursarray[ind,i_neighbour];  
                Amat[i_neighbour,1]=cellcentertocellcenterx[ind,i_neighbour]
                Amat[i_neighbour,2]=cellcentertocellcentery[ind,i_neighbour]
                bvec[i_neighbour]=p_old[i_A]-p_old[i_P];
            end
            #xvec=Amat[1:len_cellneighboursline,:]\bvec[1:len_cellneighboursline];
            #dpdx=xvec[1];
            #dpdy=xvec[2];

            if len_cellneighboursline>1
                Aplus=transpose(Amat)*Amat;
                a=Aplus[1,1]
                b=Aplus[1,2]
                c=Aplus[2,1]
                d=Aplus[2,2] 
                bvec_mod=transpose(Amat)*bvec
                inv = 1/(a * d - b * c)
                # 1 / (ad -bc) * [d -b; -c a]
                dpdx = inv * d * bvec_mod[1] - inv * b * bvec_mod[2]
                dpdy = -inv * c * bvec_mod[1] + inv * a * bvec_mod[2]
            else
                dpdx=0;
                dpdy=0;
            end

        end
        return dpdx,dpdy
    end

    function numerical_flux_function(i_method,vars_P,vars_A,meshparameters);
        if i_method==1;
            #first order upwinding
            rho_P=vars_P[1];
            u_P=vars_P[2];
            v_P=vars_P[3];
            gamma_P=vars_P[4];
            rho_A=vars_A[1];
            u_A=vars_A[2];
            v_A=vars_A[3];
            gamma_A=vars_A[4];
            n_x=meshparameters[1];
            n_y=meshparameters[2];
            A=meshparameters[3];
            n_dot_rhou=dot([n_x; n_y],0.5*(rho_P+rho_A)*[0.5*(u_P+u_A); 0.5*(v_P+v_A)]);
            phi=1;
            F_rho_num_add=n_dot_rhou*phi*A;
            if n_dot_rhou>=0;
                phi=u_P;                                
            else
                phi=u_A;
            end
            F_u_num_add=n_dot_rhou*phi*A;     
            if n_dot_rhou>=0;
                phi=v_P;  
            else
                phi=v_A;
            end
            F_v_num_add=n_dot_rhou*phi*A; 
            n_dot_u=dot([n_x; n_y],[0.5*(u_P+u_A); 0.5*(v_P+v_A)]);
            if n_dot_u>=0; 
                phi=gamma_P;  
            else
                phi=gamma_A;
            end  
            F_gamma_num_add=n_dot_u*phi*A;
            phi=1;
            F_gamma_num1_add=n_dot_u*phi*A;
        end
        return F_rho_num_add,F_u_num_add,F_v_num_add,F_gamma_num_add,F_gamma_num1_add
    end

    function numerical_flux_function_boundary(i_method,vars_P,vars_A,meshparameters,n_dot_u);
        if i_method==1;
            #first order upwinding
            rho_P=vars_P[1];
            u_P=vars_P[2];
            v_P=vars_P[3];
            gamma_P=vars_P[4];
            rho_A=vars_A[1];
            u_A=vars_A[2];
            v_A=vars_A[3];
            gamma_A=vars_A[4];
            n_x=meshparameters[1];
            n_y=meshparameters[2];
            A=meshparameters[3];        
            n_dot_rhou=n_dot_u*0.5*(rho_A+rho_P);
            phi=1;
            F_rho_num_add=n_dot_rhou*phi*A;
            if n_dot_u<=0 
                phi=u_A;                   
            else
                phi=u_P;
            end
            F_u_num_add=n_dot_rhou*phi*A;
            if n_dot_u<=0 
                phi=v_A;                   
            else
                phi=v_P;
            end
            F_v_num_add=n_dot_rhou*phi*A;
            if n_dot_u<=0 
                phi=gamma_A;                   
            else
                phi=gamma_P;
            end
            F_gamma_num_add=n_dot_u*phi*A;
            phi=1;
            F_gamma_num1_add=n_dot_u*phi*A;
        end
        return F_rho_num_add,F_u_num_add,F_v_num_add,F_gamma_num_add,F_gamma_num1_add
    end

    function delete_files();
        #delete the intermediate output files
        rm.(glob("output_*.jld2"))
    end

    function read_mesh(meshfilename,paramset,paramset1,paramset2,paramset3,paramset4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_interactive,r_p);
        #read mesh file and prepare to be used in solver:
        # - number of cells, cell ids start with 1
        # - x,y,z-coordinates of the nodes
        # - x,y,z-coordinates of the cell centers
        # - patch properties
        #
        #read other mesh files than Nastran bulk data format (bdf) based on extension and calculate the required mesh data or convert to Nastran format prepare with existing function                
        
        #read Nastran mesh
        if meshfilename[end-2:end]=="bdf"
           N,cellgridid,gridx,gridy,gridz,cellcenterx,cellcentery,cellcenterz,patchparameters,patchparameters1,patchparameters2,patchparameters3,patchparameters4,patchids1,patchids2,patchids3,patchids4,inletpatchids=
                read_nastran_mesh(meshfilename,paramset,paramset1,paramset2,paramset3,paramset4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_interactive,r_p);
        end
   
        return N,cellgridid,gridx,gridy,gridz,cellcenterx,cellcentery,cellcenterz,patchparameters,patchparameters1,patchparameters2,patchparameters3,patchparameters4,patchids1,patchids2,patchids3,patchids4,inletpatchids
    end

    function read_nastran_mesh(meshfilename,paramset,paramset1,paramset2,paramset3,paramset4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_interactive,r_p);
        #Nastran format fixed length (8 digits), GRIDS defined in global CS (local CS field empty)

        if ~isfile(meshfilename);
            errorstring=string("File ",meshfilename," not existing"* "\n"); 
            error(errorstring);
        end
        ind=Int64(1);
        gridind=Int64(1);
        setind=Int64(1);
        issetdefinition=Int64(0);
        patchorigids1=[];
        patchorigids2=[];
        patchorigids3=[];
        patchorigids4=[];
        origgridid=[];
        gridx=[];
        gridy=[];
        gridz=[];
        celloriggridid=[];
        cellgridid=Array{Int64}(undef, 0, 3);
        inletpatchids=[];
        
        open(meshfilename, "r") do fid
            line=1;
            while !eof(fid)
                thisline=readline(fid)
                if length(thisline)>=8;
                    if issetdefinition==1; 
                        if cmp( thisline[1:8],"        ")!=0;  #check if the first eight characters are empty, else issetdefinition=0;
                             issetdefinition=Int64(0);
                             setind=setind+1;
                        end
                    end
                    card=thisline[1:8];
                    if cmp(card,"GRID    ")==0
                        gridindstring=thisline[9:16];
                        origgridid=vcat(origgridid,parse(Int64,gridindstring));                        
                        txt=thisline[25:32];
                        txt=replace(txt," "=> "");txt=replace(txt,"E" => "");txt=replace(txt,"e" => "");
                        txt1=replace(txt,"-" => "e-");txt1=replace(txt1,"+" => "e+");
                        if cmp(txt1[1],'e')==0;txt2=txt1[2:end];else;txt2=txt1;end;
                        val=parse(Float64,txt2);
                        val1=val;
                        txt=thisline[33:40];
                        txt=replace(txt," "=> "");txt=replace(txt,"E" => "");txt=replace(txt,"e" => "");
                        txt1=replace(txt,"-" => "e-");txt1=replace(txt1,"+" => "e+");
                        if cmp(txt1[1],'e')==0;txt2=txt1[2:end];else;txt2=txt1;end;
                        val=parse(Float64,txt2);
                        val2=val;
                        txt=thisline[41:48];
                        txt=replace(txt," "=> "");txt=replace(txt,"E" => "");txt=replace(txt,"e" => "");
                        txt1=replace(txt,"-" => "e-");txt1=replace(txt1,"+" => "e+");
                        if cmp(txt1[1],'e')==0;txt2=txt1[2:end];else;txt2=txt1;end;
                        val=parse(Float64,txt2);
                        val3=val;
                        gridx=vcat(gridx,Float64(val1));
                        gridy=vcat(gridy,Float64(val2));
                        gridz=vcat(gridz,Float64(val3));
                        gridind=gridind+1;
                    elseif cmp(card,"CTRIA3  ")==0;        
                        celloriggridid=vcat(celloriggridid,parse(Int64,thisline[9:16]));
                        i1val=parse(Int64,thisline[25:32]);
                        i1=findfirst(isequal(i1val),origgridid);
                        i2val=parse(Int64,thisline[33:40]);
                        i2=findfirst(isequal(i2val),origgridid);
                        i3val=parse(Int64,thisline[41:48]);
                        i3=findfirst(isequal(i3val),origgridid);
                        ivec=[i1,i2,i3];
                        idel=findall(isequal(min(ivec[1],ivec[2],ivec[3])),ivec);
                        deleteat!(ivec,idel)
                        idel=findall(isequal(max(ivec[1],ivec[2])),ivec);
                        deleteat!(ivec,idel)                        
                        cellgridid=vcat(cellgridid,[min(i1,i2,i3) ivec[1] max(i1,i2,i3)]);
                        ind=ind+1;    
                    elseif cmp( card[1:3],"SET")==0 || issetdefinition==1;
                        issetdefinition=1;
                        txt1=thisline[9:end];
                        txt1=replace(txt1," "=> "");
                        txt2=split(txt1,",");
                        for i in 1:length(txt2);
                            if !isempty(txt2[i]);
                                if setind==1; 
                                    patchorigids1=vcat(patchorigids1,parse(Int64,txt2[i]))
                                elseif setind==2;
                                    patchorigids2=vcat(patchorigids2,parse(Int64,txt2[i]))
                                elseif setind==3;
                                    patchorigids3=vcat(patchorigids3,parse(Int64,txt2[i]))
                                elseif setind==4;
                                    patchorigids4=vcat(patchorigids4,parse(Int64,txt2[i]))
                                end
                            end
                        end
                    end
                end
                line+=1
            end
        end
        N=ind-1;  #total number of cells
        
        #loop to define cell center coordinates in global CS
        cellcenterx=[];
        cellcentery=[];
        cellcenterz=[];
        for ind in 1:N;
            i1=cellgridid[ind,1];
            i2=cellgridid[ind,2];
            i3=cellgridid[ind,3];
            cellcenterx=vcat(cellcenterx,(gridx[i1]+gridx[i2]+gridx[i3])/3);
            cellcentery=vcat(cellcentery,(gridy[i1]+gridy[i2]+gridy[i3])/3);
            cellcenterz=vcat(cellcenterz,(gridz[i1]+gridz[i2]+gridz[i3])/3);
        end

        if i_interactive==1;
            assign_pset(r_p,N,cellcenterx,cellcentery,cellcenterz)
            psetfilename="pset.jld2"
            if ~isfile(psetfilename);
                errorstring=string("File ",psetfilename," not existing"* "\n"); 
                error(errorstring);
            end
            @load psetfilename pset;
            inletpatchids=pset;
            if length(inletpatchids)<1;
                errorstring=string("Inlet definition empty"* "\n"); 
                error(errorstring);
            end
            patchids1=[];
            patchids2=[];
            patchids3=[];
            patchids4=[];   
            patchparameters=paramset;
            patchparameters1=[];
            patchparameters2=[];
            patchparameters3=[];
            patchparameters4=[];        
        else
            patchids1=[];
            patchids2=[];
            patchids3=[];
            patchids4=[];
            for i in 1:length(patchorigids1);
                i1=findfirst(isequal(patchorigids1[i]),celloriggridid);
                patchids1=vcat(patchids1,i1);
            end
            for i=1:length(patchorigids2);
                i1=findfirst(isequal(patchorigids2[i]),celloriggridid);
                patchids2=vcat(patchids2,i1);
            end
            for i=1:length(patchorigids3);
                i1=findfirst(isequal(patchorigids3[i]),celloriggridid);
                patchids3=vcat(patchids3,i1);
            end
            for i=1:length(patchorigids4);
                i1=findfirst(isequal(patchorigids4[i]),celloriggridid);
                patchids4=vcat(patchids4,i1);
            end
            if i_interactive==2;
                assign_pset(r_p,N,cellcenterx,cellcentery,cellcenterz)
                psetfilename="pset.jld2"
                if ~isfile(psetfilename);
                    errorstring=string("File ",psetfilename," not existing"* "\n"); 
                    error(errorstring);
                end
                @load psetfilename pset;
                inletpatchids=pset;
                if length(patchids1)<1;
                    errorstring=string("Inlet definition empty"* "\n"); 
                    error(errorstring);
                end
            end
            patchparameters=paramset;
            patchparameters1=[];
            patchparameters2=[];
            patchparameters3=[];
            patchparameters4=[];
            for i_patch in 1:4;
                if i_patch==1;
                    patchids=patchids1;
                elseif i_patch==2;
                    patchids=patchids2;
                elseif i_patch==3;
                    patchids=patchids3;
                elseif i_patch==4;
                    patchids=patchids4;
                end
                if !isempty(patchids);
                    if i_patch==1;
                        if patchtype1val==2;
                            patchparameters1=paramset1;
                        end
                    elseif i_patch==2;
                        if patchtype2val==2;
                            patchparameters2=paramset2;
                        end
                    elseif i_patch==3;
                        if patchtype3val==2; 
                            patchparameters3=paramset3;
                        end
                    elseif i_patch==4;
                        if patchtype4val==2;
                            patchparameters4=paramset4;
                        end
                    end
                end
            end
        end
           
        return N,cellgridid,gridx,gridy,gridz,cellcenterx,cellcentery,cellcenterz,patchparameters,patchparameters1,patchparameters2,patchparameters3,patchparameters4,patchids1,patchids2,patchids3,patchids4,inletpatchids
    end


    function create_faces(cellgridid, N, maxnumberofneighbours);
        # Find the set with the IDs of the neighbouring cells
        # and identify wall cells

        celltype=Vector{Int64}(undef, N);
        for i in 1:N;
            celltype[i]=1;
        end
        faces=Array{Int64}(undef, 0, 3);   #three columns: grid id1, grid id2, cell id
        i=1;
        for ind=1:N
            i1=cellgridid[ind,1];
            i2=cellgridid[ind,2];    
            faces=vcat(faces,[min(i1,i2) max(i1,i2) ind]);
            i=i+1;
            i1=cellgridid[ind,2];
            i2=cellgridid[ind,3];    
            faces=vcat(faces,[min(i1,i2) max(i1,i2) ind]);
            i=i+1;
            i1=cellgridid[ind,3];
            i2=cellgridid[ind,1];    
            faces=vcat(faces,[min(i1,i2) max(i1,i2) ind]);
            i=i+1;
        end
        facessorted=sortslices(faces,dims=1);
        vals1=unique(facessorted[:,1]);  

        # this must be generalized, currently only hard-coded number of neighbouring cells of a tria is possible
        # all considered cases had <<10 neighbouring cells 
        cellneighboursarray=Array{Int64}(undef, N, maxnumberofneighbours);
        for ind in 1:N;
            for ind_n in 1:maxnumberofneighbours;
                 cellneighboursarray[ind,ind_n]=-9;
            end
        end

        for i in 1:length(vals1);
            inds2=findall(isequal(vals1[i]), facessorted[:,1]);
            facesdetail_unsorted=facessorted[inds2,2:3];
            facesdetail=sortslices(facesdetail_unsorted,dims=1);
            for j=1:size(facesdetail,1);
                i1=facesdetail[j,2];
                inds3=findall(isequal(facesdetail[j,1]),facesdetail[:,1]);
                inds4=findall(!isequal(j),inds3);
                inds5=inds3[inds4];
                if isempty(inds5);
                    celltype[i1]=-3;  #wall
                else
                    if j==1;
                        for k in 1:length(inds5);
                            matrixrow=cellneighboursarray[i1,:];
                            indcolumn=findfirst(isequal(-9),matrixrow); 
                            if isnothing(indcolumn)
                                error("More than 10 neighbours of one tria is not supported \n");
                            else
                                cellneighboursarray[i1,indcolumn]=facesdetail[inds5[k],2];
                            end
                        end
                    else
                       for k in 1:1; 
                            matrixrow=cellneighboursarray[i1,:];
                            indcolumn=findfirst(isequal(-9),matrixrow); 
                            if isnothing(indcolumn)
                                error("More than 10 neighbours of one tria is not supported"* "\n");
                            else
                                cellneighboursarray[i1,indcolumn]=facesdetail[inds5[k],2];
                            end
                        end
                    end
                end
            end
        end

        return faces, cellneighboursarray, celltype
    end

    function assign_parameters(i_interactive,celltype,patchparameters0,patchparameters1,patchparameters2,patchparameters3,patchparameters4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,patchids1,patchids2,patchids3,patchids4,inletpatchids,mu_resin_val,N);
        #assign properties to cells

        cellthickness=Vector{Float64}(undef, N);
        cellporosity=Vector{Float64}(undef, N);
        cellpermeability=Vector{Float64}(undef, N);
        cellalpha=Vector{Float64}(undef, N);
        celldirection=Array{Float64}(undef, N,3);
        cellviscosity=Vector{Float64}(undef, N);

        if i_interactive==0 || i_interactive==2;
            if patchtype1val==1;
                for i in 1:N;
                    ind=findfirst(isequal(i),patchids1);
                    if ~isnothing(ind)
                        celltype[i]=-1;
                    end
                end                
            elseif patchtype1val==3;
                for i in 1:N;
                    ind=findfirst(isequal(i),patchids1);
                    if ~isnothing(ind)
                        celltype[i]=-2;
                    end
                end                  
            end
            if patchtype2val==1;
                for i in 1:N;
                    ind=findfirst(isequal(i),patchids2);
                    if ~isnothing(ind)
                        celltype[i]=-1;
                    end
                end
            elseif patchtype2val==3;
                for i in 1:N;
                    ind=findfirst(isequal(i),patchids2);
                    if ~isnothing(ind)
                        celltype[i]=-2;
                    end
                end
            end
            if patchtype3val==1;
                for i in 1:N;
                    ind=findfirst(isequal(i),patchids3);
                    if ~isnothing(ind)
                        celltype[i]=-1;
                    end
                end
            elseif patchtype3val==3;
                for i in 1:N;
                    ind=findfirst(isequal(i),patchids3);
                    if ~isnothing(ind)
                        celltype[i]=-2;
                    end
                end
            end
            if patchtype4val==1;
                for i in 1:N;
                    ind=findfirst(isequal(i),patchids4);
                    if ~isnothing(ind)
                        celltype[i]=-1;
                    end
                end
            elseif patchtype4val==3;
                for i in 1:N;
                    ind=findfirst(isequal(i),patchids4);
                    if ~isnothing(ind)
                        celltype[i]=-2;
                    end
                end
            end
        end
        if i_interactive==1 || i_interactive==2;
            for i in 1:N;
                ind=findfirst(isequal(i),inletpatchids);
                if ~isnothing(ind)
                    celltype[i]=-1;
                end
            end
        end        
        for ind in 1:N;
            if i_interactive==1;
                #thickness
                cellthickness[ind]=patchparameters0[2];
                #porosity
                cellporosity[ind]=patchparameters0[1]; 
                #isotropic permeability 
                cellpermeability[ind]=patchparameters0[3];            
                #alpha permeability 
                cellalpha[ind]=patchparameters0[4];
                #primary direction
                vec=[patchparameters0[5] patchparameters0[6] patchparameters0[7]];
                celldirection[ind,:]=vec/sqrt(dot(vec,vec));
                #viscosity
                cellviscosity[ind]=mu_resin_val; 
            else
                ind1=findfirst(isequal(ind),patchids1);
                ind2=findfirst(isequal(ind),patchids2);
                ind3=findfirst(isequal(ind),patchids3);
                ind4=findfirst(isequal(ind),patchids4);
                if (patchtype1val==2 && ~isnothing(ind1)) 
                    patchparameters=patchparameters1;
                elseif (patchtype2val==2 && ~isnothing(ind2)) 
                    patchparameters=patchparameters2;
                elseif (patchtype3val==2 && ~isnothing(ind3)) 
                    patchparameters=patchparameters3;
                elseif (patchtype4val==2 && ~isnothing(ind4)) 
                    patchparameters=patchparameters4;
                end
                if (patchtype1val==2 && issubset(ind,patchids1)) || (patchtype2val==2 && issubset(ind,patchids2)) || (patchtype3val==2 && issubset(ind,patchids3)) || (patchtype4val==2 && issubset(ind,patchids4));
                    #thickness
                    cellthickness[ind]=patchparameters[2];
                    #porosity
                    cellporosity[ind]=patchparameters[1]; 
                    #isotropic permeability 
                    cellpermeability[ind]=patchparameters[3];            
                    #alpha permeability 
                    cellalpha[ind]=patchparameters[4];
                    #primary direction
                    vec=[patchparameters[5] patchparameters[6] patchparameters[7]];
                    celldirection[ind,:]=vec/sqrt(dot(vec,vec));
                    #viscosity
                    cellviscosity[ind]=mu_resin_val; 
                else
                    #thickness
                    cellthickness[ind]=patchparameters0[2];
                    #porosity
                    cellporosity[ind]=patchparameters0[1]; 
                    #isotropic permeability 
                    cellpermeability[ind]=patchparameters0[3];            
                    #alpha permeability 
                    cellalpha[ind]=patchparameters0[4];
                    #primary direction
                    vec=[patchparameters0[5] patchparameters0[6] patchparameters0[7]];
                    celldirection[ind,:]=vec/sqrt(dot(vec,vec));
                    #viscosity
                    cellviscosity[ind]=mu_resin_val; 
                end
            end
        end

        return cellthickness, cellporosity, cellpermeability, cellalpha, celldirection, cellviscosity, celltype
    end

    function create_coordinate_systems(N, cellgridid, gridx, gridy, gridz, cellcenterx,cellcentery,cellcenterz, faces, cellneighboursarray, celldirection, cellthickness, maxnumberofneighbours);
        # define the local cell coordinate system and 
        # the transformation matrix from the local cell coordinate system from the neighbouring cell to the local cell coordinate system of the considered cell

        cellvolume=Vector{Float64}(undef, N);
        cellcentertocellcenterx=Array{Float64}(undef, N, maxnumberofneighbours);
        cellcentertocellcentery=Array{Float64}(undef, N, maxnumberofneighbours);
        T11=Array{Float64}(undef, N, maxnumberofneighbours);
        T12=Array{Float64}(undef, N, maxnumberofneighbours);
        T21=Array{Float64}(undef, N, maxnumberofneighbours);
        T22=Array{Float64}(undef, N, maxnumberofneighbours);
        cellfacenormalx=Array{Float64}(undef, N, maxnumberofneighbours);
        cellfacenormaly=Array{Float64}(undef, N, maxnumberofneighbours);
        cellfacearea=Array{Float64}(undef, N, maxnumberofneighbours);
        for ind in 1:N;
            cellvolume[ind]=-9;
            for ind_n in 1:maxnumberofneighbours;
                cellcentertocellcenterx[ind,ind_n]=-9.0;
                cellcentertocellcentery[ind,ind_n]=-9.0;
                T11[ind,ind_n]=-9.0;
                T12[ind,ind_n]=-9.0;
                T21[ind,ind_n]=-9.0;
                T22[ind,ind_n]=-9.0;
                cellfacenormalx[ind,ind_n]=-9.0;
                cellfacenormaly[ind,ind_n]=-9.0;
                cellfacearea[ind,ind_n]=-9.0;
            end
        end
        b1=Array{Float64}(undef, N, 3);
        b2=Array{Float64}(undef, N, 3);
        b3=Array{Float64}(undef, N, 3);
        gridxlocal=Array{Float64}(undef, N, 3);
        gridylocal=Array{Float64}(undef, N, 3);
        gridzlocal=Array{Float64}(undef, N, 3);
        theta=Vector{Float64}(undef, N);

        for ind in 1:N;
            # First, an intermediate orthonormal basis {b1, b2, b3} is created with 
            # first direction pointing from node with smallest ID to node with medium ID, 
            # second direction pointing in the orthogonal component from node with smallest ID to node with highest ID and 
            # third direction given by the cross product of the first two directions. 
            # The origin of the local coordinate system is the geometric center of the triangular cell. 
            i1=cellgridid[ind,1];
            i2=cellgridid[ind,2];
            i3=cellgridid[ind,3];  
            b1[ind,1:3]=[gridx[i2]-gridx[i1] gridy[i2]-gridy[i1] gridz[i2]-gridz[i1]];
            b1[ind,1:3]=b1[ind,1:3]/sqrt(dot(b1[ind,1:3],b1[ind,1:3]));
            a2=[gridx[i3]-gridx[i1] gridy[i3]-gridy[i1] gridz[i3]-gridz[i1]]';
            a2=a2/sqrt(dot(a2,a2));
            b2[ind,1:3]=a2-dot(b1[ind,1:3],a2)/dot(b1[ind,1:3],b1[ind,1:3])*b1[ind,1:3];
            b2[ind,1:3]=b2[ind,1:3]/sqrt(dot(b2[ind,1:3],b2[ind,1:3]));
            b3[ind,1:3]=cross(b1[ind,1:3],b2[ind,1:3]);   

            # Then the reference vector is formulated in the intermediate orthonormal basis 
            Tmat=[b1[ind,1] b2[ind,1] b3[ind,1]; b1[ind,2] b2[ind,2] b3[ind,2]; b1[ind,3] b2[ind,3] b3[ind,3]];
            xvec=celldirection[ind,:];
            bvec=Tmat\xvec;
            r1=[bvec[1] bvec[2] bvec[3]]';  #ref dir in local CS

            # In order to get the local coordinate system the basis {b1, b2, b3} is rotated by angle theta about the b3-axis.
            # Calculate the angle by which b1 must be rotated about the b3-axis to match r1 via relation rotation matrix Rz(theta)*[1;0;0]=r1, i.e. cos(theta)=r1(1) and sin(theta)=r1(2);
            theta[ind]=atan(r1[2],r1[1]);
            #Rotation of theta about nvec=b3 to get c1 and c2 
            nvec=b3[ind,:];
            xvec=b1[ind,:];
            c1=nvec*dot(nvec,xvec)+cos(theta[ind])*cross(cross(nvec,xvec),nvec)+sin(theta[ind])*cross(nvec,xvec);
            xvec=b2[ind,:];
            c2=nvec*dot(nvec,xvec)+cos(theta[ind])*cross(cross(nvec,xvec),nvec)+sin(theta[ind])*cross(nvec,xvec);
            xvec=b3[ind,:];
            c3=nvec*dot(nvec,xvec)+cos(theta[ind])*cross(cross(nvec,xvec),nvec)+sin(theta[ind])*cross(nvec,xvec);
            b1[ind,:]=c1;
            b2[ind,:]=c2;
            b3[ind,:]=c3;  
        
            #transformation of vertices into local CS
            gridxlocal[ind,1]=gridx[i1]-cellcenterx[ind];
            gridylocal[ind,1]=gridy[i1]-cellcentery[ind];
            gridzlocal[ind,1]=gridz[i1]-cellcenterz[ind];
            gridxlocal[ind,2]=gridx[i2]-cellcenterx[ind];
            gridylocal[ind,2]=gridy[i2]-cellcentery[ind];
            gridzlocal[ind,2]=gridz[i2]-cellcenterz[ind];
            gridxlocal[ind,3]=gridx[i3]-cellcenterx[ind];
            gridylocal[ind,3]=gridy[i3]-cellcentery[ind];
            gridzlocal[ind,3]=gridz[i3]-cellcenterz[ind];
            Tmat=[b1[ind,1] b2[ind,1] b3[ind,1]; b1[ind,2] b2[ind,2] b3[ind,2]; b1[ind,3] b2[ind,3] b3[ind,3]];
            xvec=[gridxlocal[ind,1] gridylocal[ind,1] gridzlocal[ind,1]]'; 
            bvec=Tmat\xvec;
            gridxlocal[ind,1]=bvec[1];gridylocal[ind,1]=bvec[2];gridzlocal[ind,1]=bvec[3];
            xvec=[gridxlocal[ind,2] gridylocal[ind,2] gridzlocal[ind,2]]'; 
            bvec=Tmat\xvec;
            gridxlocal[ind,2]=bvec[1];gridylocal[ind,2]=bvec[2];gridzlocal[ind,2]=bvec[3];
            xvec=[gridxlocal[ind,3] gridylocal[ind,3] gridzlocal[ind,3]]'; 
            bvec=Tmat\xvec;
            gridxlocal[ind,3]=bvec[1];gridylocal[ind,3]=bvec[2];gridzlocal[ind,3]=bvec[3];
        end

        cellids=[Int64(-9) Int64(-9)];
        gridids=[Int64(-9) Int64(-9)];
        x=[Float64(-9.0), Float64(-9.0), Float64(-9.0)];
        x0=[Float64(-9.0), Float64(-9.0), Float64(-9.0)];
        r0=[Float64(-9.0), Float64(-9.0), Float64(-9.0)];
        gridxlocal_neighbour=[Float64(-9.0), Float64(-9.0), Float64(-9.0)];
        gridylocal_neighbour=[Float64(-9.0), Float64(-9.0), Float64(-9.0)];
        gridzlocal_neighbour=[Float64(-9.0), Float64(-9.0), Float64(-9.0)];
        f1=[Float64(-9.0), Float64(-9.0), Float64(-9.0)];
        f2=[Float64(-9.0), Float64(-9.0), Float64(-9.0)];
        f3=[Float64(-9.0), Float64(-9.0), Float64(-9.0)];
        # In a next step the flattened geometry is created, i.e. the cell center and
        # the non-common node of the neighbouring cell is rotated about
        # the common edge to lie in the plane of the considered cell with ID ind
        for ind in 1:N;
            cellneighboursline=cellneighboursarray[ind,:];
            cellneighboursline=cellneighboursline[cellneighboursline .> 0]
            for i_neighbour in 1:length(cellneighboursline);
                # Find first the cell center of neighbouring cell in local coordinate system of cell ind
                # 1) projection of cell center P=(0,0) onto straigth line through
                #    i1 and i2 to get point Q1 and calculation of length l1 of line
                #    segment PQ1
                # 2) projection of neighbouring cell center A onto straight line
                #    through i1 and i2 to get point Q2 in global coordinate system,
                #    calculatin of length l2 of segment AQ2 and then
                #    transformation of Q2 into local coordinate system and then
                #    cellcentertocellcenterx/y(ind,1) is given by vector addition
                #    PQ1+Q1Q2+l2/l1*PQ1
            
                #for every neighbour find the two indices belonging to the boundary
                #face in between; face direction is from smaller to larger index
                #x0..local coordinates of smaller index
                #r0..vector from smaller to larger index in LCS
                inds1=findall(isequal(ind),faces[:,3]);
                inds2=findall(isequal(cellneighboursline[i_neighbour]),faces[:,3]);
                mat1=faces[inds1,:];
                mat2=faces[inds2,:];
                mat3=vcat(mat1,mat2); 
                mat4=sortslices(mat3,dims=1);
                for irow in 1:size(mat4,1)-1;
                    if mat4[irow,1]==mat4[irow+1,1] && mat4[irow,2]==mat4[irow+1,2];
                        if mat4[irow,3]==ind
                            cellids=[ind mat4[irow+1,3]];
                        else
                            cellids=[ind mat4[irow,3]];
                        end
                        gridids=[mat4[irow,1] mat4[irow,2]];
                    end
                end
                inds=[cellgridid[ind,1], cellgridid[ind,2], cellgridid[ind,3]];
                ia=findall(isequal(gridids[1]),inds);
                ib=findall(isequal(gridids[2]),inds);
                x0=[gridxlocal[ind,ia], gridylocal[ind,ia], gridzlocal[ind,ia]];
                r0=[gridxlocal[ind,ib]-gridxlocal[ind,ia], gridylocal[ind,ib]-gridylocal[ind,ia], gridzlocal[ind,ib]-gridzlocal[ind,ia]];

                #Define xvec as the vector between cell centers ind and neighbouring cell center (A) (in GCS) 
                #and transform xvec into local coordinates bvec, this gives A in LCS.
                #Find normal distance from A in LCS to the cell boundary with that cell center A in flat geometry and 
                #face normal vector can be defined.
                x=[[0.0], [0.0], [0.0]];  #P at origin of local CS
                Px=x[1];
                Py=x[2];
                Pz=x[3];
                lambda=dot(x-x0,r0)/dot(r0,r0);  
                Q1x=x0[1]+lambda*r0[1];
                Q1y=x0[2]+lambda*r0[2];
                Q1z=x0[3]+lambda*r0[3];
                vec1=[Px-Q1x, Py-Q1y, Pz-Q1z];
                l1=sqrt(dot(vec1,vec1)); 
                nvec=[(Q1x-Px), (Q1y-Py), (Q1z-Pz)];
                nvec=nvec/sqrt(dot(nvec,nvec));
                cellfacenormalx[ind,i_neighbour]=only(nvec[1]);
                cellfacenormaly[ind,i_neighbour]=only(nvec[2]); 

                Tmat=[b1[ind,1] b2[ind,1] b3[ind,1]; b1[ind,2] b2[ind,2] b3[ind,2]; b1[ind,3] b2[ind,3] b3[ind,3]];
                xvec=[cellcenterx[cellneighboursarray[ind,i_neighbour]]-cellcenterx[ind], cellcentery[cellneighboursarray[ind,i_neighbour]]-cellcentery[ind], cellcenterz[cellneighboursarray[ind,i_neighbour]]-cellcenterz[ind] ];  #A in global CS
                bvec=Tmat\xvec;
                x=[[bvec[1]], [bvec[2]], [bvec[3]]]; #A in local CS
                Ax=x[1];
                Ay=x[2];
                Az=x[3];
                lambda=dot(x-x0,r0)/dot(r0,r0);
                Q2x=x0[1]+lambda*r0[1];
                Q2y=x0[2]+lambda*r0[2];
                Q2z=x0[3]+lambda*r0[3];
                vec2=[Ax-Q2x, Ay-Q2y, Az-Q2z];
                l2=sqrt(dot(vec2,vec2));
                cellcentertocellcenterx[ind,i_neighbour]=only(Px+(Q1x-Px)+(Q2x-Q1x)+l2/l1*(Q1x-Px));
                cellcentertocellcentery[ind,i_neighbour]=only(Py+(Q1y-Py)+(Q2y-Q1y)+l2/l1*(Q1y-Py));

                vec3=[gridxlocal[ind,ib]-gridxlocal[ind,ia], gridylocal[ind,ib]-gridylocal[ind,ia], gridzlocal[ind,ib]-gridzlocal[ind,ia]];
                cellfacearea[ind,i_neighbour]=0.5*(cellthickness[cellids[1]]+cellthickness[cellids[2]])*sqrt(dot(vec3,vec3));

                #Transformation matrix for (u,v) of neighbouring cells to local coordinate system.
                #Find the two common grid points and the third non-common grid point               
                ind21=-9;  #Issues with setdiff, therefore manual implementation  #setdiff(cellgridid[cellids[2],:],gridids)
                for ind_tmp in 1:3
                    if cellgridid[cellids[2],ind_tmp]!=gridids[1] && cellgridid[cellids[2],ind_tmp]!=gridids[2]
                        ind21=cellgridid[cellids[2],ind_tmp];
                    end
                end     
                thirdgrid=only(findall(isequal(ind21),cellgridid[cellids[2],:]));
                common1grid=only(findall(isequal(gridids[1]),cellgridid[cellids[2],:]));
                common2grid=only(findall(isequal(gridids[2]),cellgridid[cellids[2],:]));   
                #construction of the third one in outside normal direction for the flat geometry
                #based on the length of the two non-common edges
                gridxlocal_neighbour[2]=only(gridxlocal[ind,ia]);  #gridxlocal(ind,common1grid);
                gridxlocal_neighbour[3]=only(gridxlocal[ind,ib]);  #gridxlocal(ind,common2grid);
                gridylocal_neighbour[2]=only(gridylocal[ind,ia]);  #gridylocal(ind,common1grid);
                gridylocal_neighbour[3]=only(gridylocal[ind,ib]);  #gridylocal(ind,common2grid);
                gridzlocal_neighbour[2]=0.0;
                gridzlocal_neighbour[3]=0.0;
                
                ind3=-9;
                for ind_tmp in 1:3
                    if cellgridid[cellids[2],ind_tmp]!=cellgridid[ind,1] && cellgridid[cellids[2],ind_tmp]!=cellgridid[ind,2] && cellgridid[cellids[2],ind_tmp]!=cellgridid[ind,3]
                        ind3=cellgridid[cellids[2],ind_tmp];
                    end
                end
                Tmat=[b1[ind,1] b2[ind,1] b3[ind,1]; b1[ind,2] b2[ind,2] b3[ind,2]; b1[ind,3] b2[ind,3] b3[ind,3]];
                xvec=[gridx[ind3]-cellcenterx[ind], gridy[ind3]-cellcentery[ind], gridz[ind3]-cellcenterz[ind]]; #A in global CS
                bvec=Tmat\xvec;
                x=[[bvec[1]], [bvec[2]], [bvec[3]]]; #A in local CS
                Ax=x[1];
                Ay=x[2];
                Az=x[3];
                lambda=dot(x-x0,r0)/dot(r0,r0);
                Q2x=x0[1]+lambda*r0[1];
                Q2y=x0[2]+lambda*r0[2];
                Q2z=x0[3]+lambda*r0[3];
                vec2=[Ax-Q2x, Ay-Q2y, Az-Q2z];
                l2=sqrt(dot(vec2,vec2));
                gridxlocal_neighbour[1]=only(Px+(Q1x-Px)+(Q2x-Q1x)+l2/l1*(Q1x-Px));
                gridylocal_neighbour[1]=only(Py+(Q1y-Py)+(Q2y-Q1y)+l2/l1*(Q1y-Py));
                gridzlocal_neighbour[1]=Float64(0);

                #Construction of LCS f1,f2,f3 according to procedure from above using the points gridxlocal_neighbour(j),gridylocal_neighbour(j)
                ivec1=[only(cellgridid[cellids[2],1]), only(cellgridid[cellids[2],2]), only(cellgridid[cellids[2],3])];                           
                min_val=min(ivec1[1],ivec1[2],ivec1[3]); 
                max_val=max(ivec1[1],ivec1[2],ivec1[3]); 
                idel1=findall(isequal(min(ivec1[1],ivec1[2],ivec1[3])),ivec1);deleteat!(ivec1,idel1);idel1=findall(isequal(max(ivec1[1],ivec1[2])),ivec1);deleteat!(ivec1,idel1);
                median_val=ivec1[1];
                if ind3==min_val; k1=1; elseif ind3==median_val; k2=1; elseif ind3==max_val; k3=1; end                
                ind4=cellgridid[cellids[2],common1grid];
                if ind4==min_val; k1=2; elseif ind4==median_val; k2=2; elseif ind4==max_val; k3=2; end             
                ind5=cellgridid[cellids[2],common2grid];
                if ind5==min_val; k1=3; elseif ind5==median_val; k2=3; elseif ind5==max_val; k3=3; end
        
                f1=[gridxlocal_neighbour[k2]-gridxlocal_neighbour[k1], gridylocal_neighbour[k2]-gridylocal_neighbour[k1], gridzlocal_neighbour[k2]-gridzlocal_neighbour[k1]];
                f1=f1/sqrt(dot(f1,f1));
                a2=[gridxlocal_neighbour[k3]-gridxlocal_neighbour[k1], gridylocal_neighbour[k3]-gridylocal_neighbour[k1], gridzlocal_neighbour[k3]-gridzlocal_neighbour[k1]];
                a2=a2/sqrt(dot(a2,a2));
                f2=a2-dot(f1,a2)/dot(f1,f1)*f1;
                f2=f2/sqrt(dot(f2,f2));
                f3=cross(f1,f2);    
        
                nvec=f3;
                xvec=f1;
                c1=nvec*dot(nvec,xvec)+cos(theta[cellneighboursarray[ind,i_neighbour]])*cross(cross(nvec,xvec),nvec)+sin(theta[cellneighboursarray[ind,i_neighbour]] )*cross(nvec,xvec);
                xvec=f2;
                c2=nvec*dot(nvec,xvec)+cos(theta[cellneighboursarray[ind,i_neighbour]])*cross(cross(nvec,xvec),nvec)+sin(theta[cellneighboursarray[ind,i_neighbour]] )*cross(nvec,xvec);
                xvec=f3;
                c3=nvec*dot(nvec,xvec)+cos(theta[cellneighboursarray[ind,i_neighbour]])*cross(cross(nvec,xvec),nvec)+sin(theta[cellneighboursarray[ind,i_neighbour]] )*cross(nvec,xvec);
                f1=c1;
                f2=c2;
                f3=c3;
                Tmat=[f1[1] f2[1] f3[1]; f1[2] f2[2] f3[2]; f1[3] f2[3] f3[3]];

                #Assign transformation matrix for the velocities in the local coordinate systems
                #(u,v)_e=T*(u,v)_f
                T11[ind,i_neighbour]=Tmat[1,1];
                T12[ind,i_neighbour]=Tmat[1,2];
                T21[ind,i_neighbour]=Tmat[2,1];
                T22[ind,i_neighbour]=Tmat[2,2];
            end

            #calculate cell volume
            vec1=[gridxlocal[ind,2]-gridxlocal[ind,1], gridylocal[ind,2]-gridylocal[ind,1], gridzlocal[ind,2]-gridzlocal[ind,1]];
            vec2=[gridxlocal[ind,3]-gridxlocal[ind,1], gridylocal[ind,3]-gridylocal[ind,1], gridzlocal[ind,3]-gridzlocal[ind,1]];
            vec3=cross(vec1,vec2);
            cellvolume[ind]=cellthickness[ind]*0.5*sqrt(dot(vec3,vec3));
        end  

        return cellvolume, cellcentertocellcenterx, cellcentertocellcentery, T11, T12, T21, T22, cellfacenormalx, cellfacenormaly, cellfacearea
    end


    function plot_mesh(meshfilename,i_mode)
        #create mesh plot with cells with i_mode==1 and
        #create mesh plots with cell center nodes with i_mode==2 for manual selection of inlet ports

        #dummy values for calling function read_mesh
        paramset=[0.5,0.3,3e-10,1.0,1.0,0.0,0.0];paramset1=paramset;paramset2=paramset;paramset3=paramset;paramset4=paramset;
        patchtype1val=-1;patchtype2val=-1;patchtype3val=-1;patchtype4val=-1;i_interactive=0;
        r_p=0.01;
        N,cellgridid,gridx,gridy,gridz,cellcenterx,cellcentery,cellcenterz,patchparameters,patchparameters1,patchparameters2,patchparameters3,patchparameters4,patchids1,patchids2,patchids3,patchids4,inletpatchids=
            read_mesh(meshfilename,paramset,paramset1,paramset2,paramset3,paramset4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_interactive,r_p);

        #for poly plot
        X=Array{Float64}(undef, 3, N);
        Y=Array{Float64}(undef, 3, N);
        Z=Array{Float64}(undef, 3, N);
        C=Array{Float32}(undef, 3, N);
        for ind in 1:N;
            X[1,ind]=gridx[cellgridid[ind,1]];
            X[2,ind]=gridx[cellgridid[ind,2]];
            X[3,ind]=gridx[cellgridid[ind,3]];
            Y[1,ind]=gridy[cellgridid[ind,1]];
            Y[2,ind]=gridy[cellgridid[ind,2]];
            Y[3,ind]=gridy[cellgridid[ind,3]];
            Z[1,ind]=gridz[cellgridid[ind,1]];
            Z[2,ind]=gridz[cellgridid[ind,2]];
            Z[3,ind]=gridz[cellgridid[ind,3]];
            C[1,ind]=1.0;
            C[2,ind]=1.0;
            C[3,ind]=1.0;
        end
        xyz = reshape([X[:] Y[:] Z[:]]', :)
        #2..for meshscatter plot
        X2=Array{Float64}(undef, 3*N);
        Y2=Array{Float64}(undef, 3*N);
        Z2=Array{Float64}(undef, 3*N);
        C2=Array{Float64}(undef, 3*N);
        for ind in 1:N;
            X2[    ind]=gridx[cellgridid[ind,1]];
            X2[  N+ind]=gridx[cellgridid[ind,2]];
            X2[2*N+ind]=gridx[cellgridid[ind,3]];
            Y2[    ind]=gridy[cellgridid[ind,1]];
            Y2[  N+ind]=gridy[cellgridid[ind,2]];
            Y2[2*N+ind]=gridy[cellgridid[ind,3]];
            Z2[    ind]=gridz[cellgridid[ind,1]];
            Z2[  N+ind]=gridz[cellgridid[ind,2]];
            Z2[2*N+ind]=gridz[cellgridid[ind,3]];
            C2[    ind]=0.0;
            C2[  N+ind]=0.0;
            C2[2*N+ind]=0.0;
        end

        #bounding box
        deltax=maximum(gridx)-minimum(gridx);
        deltay=maximum(gridy)-minimum(gridy);
        deltaz=maximum(gridz)-minimum(gridz);
        mindelta=min(deltax,deltay,deltaz);
        maxdelta=max(deltax,deltay,deltaz);
        if mindelta<maxdelta*0.001;
            eps_delta=maxdelta*0.001;
        else
            eps_delta=0;
        end 
        ax=(deltax+eps_delta)/(mindelta+eps_delta);
        ay=(deltay+eps_delta)/(mindelta+eps_delta);
        az=(deltaz+eps_delta)/(mindelta+eps_delta);
        
        if i_mode==1;
            fig = Figure(resolution=(600, 600))
            ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Mesh")
            poly!(connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C[:], strokewidth=1)
            #hidedecorations!(ax1);
            hidespines!(ax1) 
            display(fig)
        elseif i_mode==2;
            points=rand(Point3f0, length(gridx));
            for i in 1:length(gridx)
                points[i]=Point3f0(gridx[i],gridy[i],gridz[i]);
            end
            positions = Observable(points) 

            inletpos_xyz=[-9.9e9 -9.9e9 -9.9e9];
            filename="inletpostions.jld2"
            @save filename inletpos_xyz

            markersizeval=maxdelta*100;
            fig = Figure()
            ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Select inlets with p + LMB")
            p=scatter!(ax1, positions,markersize=markersizeval)
            hidedecorations!(ax1);
            hidespines!(ax1) 

            on(events(fig).mousebutton, priority = 2) do event
                if event.button == Mouse.left && event.action == Mouse.press
                    if Keyboard.p in events(fig).keyboardstate
                        plt, i = pick(fig.scene,events(fig).mouseposition[])
                        if plt == p
                            @load filename inletpos_xyz
                            t_div=100;
                            xpos=positions[][i][1];
                            ypos=positions[][i][2];
                            zpos=positions[][i][3];
                            inletpos_xyz=vcat(inletpos_xyz,[xpos ypos zpos]);
                            @save filename inletpos_xyz
                            textpos=string("(" , string(round(t_div*xpos)/t_div) , "," , string(round(t_div*ypos)/t_div) , "," , string(round(t_div*zpos)/t_div) , ")"  )
                            t1=text!(ax1,textpos,position = (xpos,ypos,zpos) ) 
                            scatter!(Point3f0(xpos,ypos,zpos),markersize=2*markersizeval,color = :black)
                            return Consume(true)
                        end
                    end
                end
                return Consume(false)
            end
            display(fig)
        end
    end 


    function plot_sets(meshfilename)
        #create a plot with the up to four cell sets defined in the mesh file

        #dummy values for calling function read_nastran_mesh
        paramset=[0.5,0.3,3e-10,1.0,1.0,0.0,0.0];paramset1=paramset;paramset2=paramset;paramset3=paramset;paramset4=paramset;
        patchtype1val=-1;patchtype2val=-1;patchtype3val=-1;patchtype4val=-1;i_interactive=0;
        r_p=0.01;
        N,cellgridid,gridx,gridy,gridz,cellcenterx,cellcentery,cellcenterz,patchparameters,patchparameters1,patchparameters2,patchparameters3,patchparameters4,patchids1,patchids2,patchids3,patchids4,inletpatchids=
            read_mesh(meshfilename,paramset,paramset1,paramset2,paramset3,paramset4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_interactive,r_p);

        if isempty(patchids1);
            n_patch=0;
            errorstring=string("No sets defined"* "\n"); 
            error(errorstring);
        else
            if isempty(patchids2);
                n_patch=1;
            else
                if isempty(patchids3);
                    n_patch=2;
                else
                    if isempty(patchids4);
                        n_patch=3;
                    else
                        n_patch=4;
                    end
                end
            end
        end
      
        #for poly plot
        X=Array{Float64}(undef, 3, N);
        Y=Array{Float64}(undef, 3, N);
        Z=Array{Float64}(undef, 3, N);
        C=Array{Float32}(undef, 3, N);
        C_patch1=Array{Float32}(undef, 3, N);
        C_patch2=Array{Float32}(undef, 3, N);
        C_patch3=Array{Float32}(undef, 3, N);
        C_patch4=Array{Float32}(undef, 3, N);
        for ind in 1:N;
            X[1,ind]=gridx[cellgridid[ind,1]];
            X[2,ind]=gridx[cellgridid[ind,2]];
            X[3,ind]=gridx[cellgridid[ind,3]];
            Y[1,ind]=gridy[cellgridid[ind,1]];
            Y[2,ind]=gridy[cellgridid[ind,2]];
            Y[3,ind]=gridy[cellgridid[ind,3]];
            Z[1,ind]=gridz[cellgridid[ind,1]];
            Z[2,ind]=gridz[cellgridid[ind,2]];
            Z[3,ind]=gridz[cellgridid[ind,3]];
            C[1,ind]=1.0;
            C[2,ind]=1.0;
            C[3,ind]=1.0;
            if issubset(ind, patchids1)
                C_patch1[1,ind]=1.0;
                C_patch1[2,ind]=1.0;
                C_patch1[3,ind]=1.0;
            else
                C_patch1[1,ind]=0.0;
                C_patch1[2,ind]=0.0;
                C_patch1[3,ind]=0.0;
            end   
            if issubset(ind, patchids2)
                C_patch2[1,ind]=1.0;
                C_patch2[2,ind]=1.0;
                C_patch2[3,ind]=1.0;
            else
                C_patch2[1,ind]=0.0;
                C_patch2[2,ind]=0.0;
                C_patch2[3,ind]=0.0;
            end         
            if issubset(ind, patchids3)
                C_patch3[1,ind]=1.0;
                C_patch3[2,ind]=1.0;
                C_patch3[3,ind]=1.0;
            else
                C_patch3[1,ind]=0.0;
                C_patch3[2,ind]=0.0;
                C_patch3[3,ind]=0.0;
            end  
            if issubset(ind, patchids4)
                C_patch4[1,ind]=1.0;
                C_patch4[2,ind]=1.0;
                C_patch4[3,ind]=1.0;
            else
                C_patch4[1,ind]=0.0;
                C_patch4[2,ind]=0.0;
                C_patch4[3,ind]=0.0;
            end        
        end
        xyz = reshape([X[:] Y[:] Z[:]]', :)

        #2..for meshscatter plot
        X2=Array{Float64}(undef, 3*N);
        Y2=Array{Float64}(undef, 3*N);
        Z2=Array{Float64}(undef, 3*N);
        C2=Array{Float64}(undef, 3*N);
        for ind in 1:N;
            X2[    ind]=gridx[cellgridid[ind,1]];
            X2[  N+ind]=gridx[cellgridid[ind,2]];
            X2[2*N+ind]=gridx[cellgridid[ind,3]];
            Y2[    ind]=gridy[cellgridid[ind,1]];
            Y2[  N+ind]=gridy[cellgridid[ind,2]];
            Y2[2*N+ind]=gridy[cellgridid[ind,3]];
            Z2[    ind]=gridz[cellgridid[ind,1]];
            Z2[  N+ind]=gridz[cellgridid[ind,2]];
            Z2[2*N+ind]=gridz[cellgridid[ind,3]];
            C2[    ind]=0.0;
            C2[  N+ind]=0.0;
            C2[2*N+ind]=0.0;
        end

        resolution_val=300;
        if n_patch==1;
            fig = Figure(resolution=(1*resolution_val, resolution_val))
        elseif n_patch==2;
            fig = Figure(resolution=(2*resolution_val, 2*resolution_val))
        elseif n_patch==3;
            fig = Figure(resolution=(3*resolution_val, resolution_val))
        elseif n_patch==4;
            fig = Figure(resolution=(4*resolution_val, resolution_val))
        end

        #bounding box
        deltax=maximum(gridx)-minimum(gridx);
        deltay=maximum(gridy)-minimum(gridy);
        deltaz=maximum(gridz)-minimum(gridz);
        mindelta=min(deltax,deltay,deltaz);
        maxdelta=max(deltax,deltay,deltaz); 
        if mindelta<maxdelta*0.001;
            eps_delta=maxdelta*0.001;
        else
            eps_delta=0;
        end 
        ax=(deltax+eps_delta)/(mindelta+eps_delta);
        ay=(deltay+eps_delta)/(mindelta+eps_delta);
        az=(deltaz+eps_delta)/(mindelta+eps_delta);
        ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Set 1")
        poly!(connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_patch1[:], strokewidth=1)
        hidedecorations!(ax1);hidespines!(ax1) 
        if n_patch>=2;
            ax2 = Axis3(fig[1, 2]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Set 2")
            poly!(connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_patch2[:], strokewidth=1)
            hidedecorations!(ax2);hidespines!(ax2) 
        end
        if n_patch>=3;
            ax3 = Axis3(fig[1, 3]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Set 3")
            poly!(connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_patch3[:], strokewidth=1)
            hidedecorations!(ax3);hidespines!(ax3) 
        end
        if n_patch>=4;
            ax4 = Axis3(fig[1, 4]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Set 4")
            poly!(connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_patch4[:], strokewidth=1)
            hidedecorations!(ax4);hidespines!(ax4) 
        end
        display(fig)
    end

    function assign_pset(r_p,N,cellcenterx,cellcentery,cellcenterz)
        #create the cell set from the manually selected inlet port nodes

        filename="inletpostions.jld2"
        @load filename inletpos_xyz
        n_p=size(inletpos_xyz,1)-1;
        patchpids=[];
        i=1;
        for i_p in 2:n_p+1;
            r_p_temp=r_p;
            i_add=1;
            for ind in 1:N
                vec1=[cellcenterx[ind]-inletpos_xyz[i_p,1],cellcentery[ind]-inletpos_xyz[i_p,2],cellcenterz[ind]-inletpos_xyz[i_p,3]]
                if sqrt(dot(vec1,vec1))<=r_p_temp;
                   patchpids=vcat(patchpids,ind);
                   i=i+1;
                   i_add=0;
                end
            end
            while i_add==1;
                r_p_temp=1.1*r_p_temp;
                i_firstcell=0;
                for ind in 1:N
                    vec1=[cellcenterx[ind]-inletpos_xyz[i_p,1],cellcentery[ind]-inletpos_xyz[i_p,2],cellcenterz[ind]-inletpos_xyz[i_p,3]]
                    if sqrt(dot(vec1,vec1))<=r_p_temp  && i_firstcell==0;
                        patchpids=vcat(patchpids,ind);
                        i_firstcell=1;
                        i=i+1;
                        i_add=0;
                    end
                end
            end
        end
        pset=patchpids;
        psetfilename="pset.jld2"
        @save psetfilename pset;
    end

    function plot_results(resultsfilename)
        #create contour plots of the filling factor and the pressure after loading a results file
        #default call: rtmsim.plot_results("results.jld2")

        if ~isfile(resultsfilename);
            errorstring=string("File ",resultsfilename," not existing"* "\n"); 
            error(errorstring);
        end
        t_digits=2; 
        t_div=10^2;
        @load resultsfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out

        gamma_plot=Vector{Float64}(undef, N);
        deltap=maximum(p_new)-minimum(p_new);     
        for ind=1:N;
            if gamma_out[ind]>0.8;
                gamma_plot[ind]=1;
            else
                gamma_plot[ind]=0;
            end
        end
        deltagamma=maximum(gamma_plot)-minimum(gamma_plot);

        #for poly plot
        inds0=findall(gamma_out.>-0.5);
        N0=length(inds0);
        X=Array{Float64}(undef, 3, N0);
        Y=Array{Float64}(undef, 3, N0);
        Z=Array{Float64}(undef, 3, N0);
        C_p=Array{Float32}(undef, 3, N0);        
        C_gamma=Array{Float32}(undef, 3, N0);
        inds1=findall(gamma_out.<-0.5);
        N1=length(inds1);
        X1=Array{Float64}(undef, 3, N1);
        Y1=Array{Float64}(undef, 3, N1);
        Z1=Array{Float64}(undef, 3, N1);  
        C1_gamma=Array{Float32}(undef, 3, N1);
        C1_p=Array{Float32}(undef, 3, N1);
        for i in 1:N0;
            ind=inds0[i];
            X[1,i]=gridx[cellgridid[ind,1]];
            X[2,i]=gridx[cellgridid[ind,2]];
            X[3,i]=gridx[cellgridid[ind,3]];
            Y[1,i]=gridy[cellgridid[ind,1]];
            Y[2,i]=gridy[cellgridid[ind,2]];
            Y[3,i]=gridy[cellgridid[ind,3]];
            Z[1,i]=gridz[cellgridid[ind,1]];
            Z[2,i]=gridz[cellgridid[ind,2]];
            Z[3,i]=gridz[cellgridid[ind,3]];
            C_gamma[1,i]=gamma_plot[ind]/deltagamma;
            C_gamma[2,i]=gamma_plot[ind]/deltagamma;
            C_gamma[3,i]=gamma_plot[ind]/deltagamma;
            C_p[1,i]=p_new[ind]/deltap;
            C_p[2,i]=p_new[ind]/deltap;
            C_p[3,i]=p_new[ind]/deltap;
        end
        xyz = reshape([X[:] Y[:] Z[:]]', :)        
        for i in 1:N1
            ind=inds1[i];
            X1[1,i]=gridx[cellgridid[ind,1]];
            X1[2,i]=gridx[cellgridid[ind,2]];
            X1[3,i]=gridx[cellgridid[ind,3]];
            Y1[1,i]=gridy[cellgridid[ind,1]];
            Y1[2,i]=gridy[cellgridid[ind,2]];
            Y1[3,i]=gridy[cellgridid[ind,3]];
            Z1[1,i]=gridz[cellgridid[ind,1]];
            Z1[2,i]=gridz[cellgridid[ind,2]];
            Z1[3,i]=gridz[cellgridid[ind,3]];
            C1_gamma[1,i]=0.5;
            C1_gamma[2,i]=0.5;
            C1_gamma[3,i]=0.5;
            C1_p[1,i]=p_new[ind]/deltap;
            C1_p[2,i]=p_new[ind]/deltap;
            C1_p[3,i]=p_new[ind]/deltap;
        end
        xyz1 = reshape([X1[:] Y1[:] Z1[:]]', :)

        #bounding box
        deltax=maximum(gridx)-minimum(gridx);
        deltay=maximum(gridy)-minimum(gridy);
        deltaz=maximum(gridz)-minimum(gridz);
        mindelta=min(deltax,deltay,deltaz);
        maxdelta=max(deltax,deltay,deltaz);
        if mindelta<maxdelta*0.001;
            eps_delta=maxdelta*0.001;
        else
            eps_delta=0;
        end 
        ax=(deltax+eps_delta)/(mindelta+eps_delta);
        ay=(deltay+eps_delta)/(mindelta+eps_delta);
        az=(deltaz+eps_delta)/(mindelta+eps_delta);

        resolution_val=600;
        fig = Figure(resolution=(2*resolution_val, resolution_val))
        ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
        poly!(connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
        if N1>0; 
            poly!(connect(xyz1, Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
        end
        hidedecorations!(ax1);
        hidespines!(ax1) 
        ax2 = Axis3(fig[1, 2]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Pressure at t=", string(round(t_div*t)/t_div) ,"s"))
        poly!(connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_p[:], strokewidth=1, colorrange=(0,1))
        if N1>0; 
            poly!(connect(xyz1, Point{3}), connect(1:length(X1), TriangleFace); color=C1_p[:], strokewidth=1, colorrange=(0,1))
        end
        Colorbar(fig[1, 3], limits = (0, deltap), colormap = :viridis,  vertical=true, height=Relative(0.5));  
        #Colorbar(fig[2, 2], limits = (0, deltap), colormap = :viridis,  vertical=false, width=Relative(0.5));
        hidedecorations!(ax2);
        hidespines!(ax2) 
        display(fig)
    end 


    function plot_overview(n_out,n_pics)
        #create filling factor contour plots
        #n_out ist the index of the last output file, if n_out==-1 the output file with the highest index is chosen
        #consider the last n_pics for creating the contour plots at four equidistant time intervals, if n_pics==-1 all available output files are considered
        #default call is plot_overview(-1,-1)

        val=0;
		n_out_start=-1;
        if n_out==-1;
            vec1=glob("output_*.jld2");
            for i=1:length(vec1);
                vec2=split(vec1[i],".")
                vec3=split(vec2[1],"_")
                val=max(val,parse(Int64,vec3[2]))
				if i==1;
				    n_out_start=parse(Int64,vec3[2]);
				end
            end            
            n_out=val;
        end
		if n_pics==-1;
		    n_pics=(n_out-n_out_start);
		end
        if mod(n_pics,4)!=0;
            errorstring=string("n_pics must be multiple of four"* "\n"); 
            error(errorstring);
        end
        t_digits=2; 
        t_div=10^2;

        resolution_val=300;
        fig = Figure(resolution=(4*resolution_val, resolution_val))
        i_out=n_out-3*Int64(n_pics/4);
        for i_plot in 1:4;          
            outputfilename=string("output_",string(i_out),".jld2");
            if ~isfile(outputfilename);
                errorstring=string("File ",outputfilename," not existing"* "\n"); 
                error(errorstring);
            else
                loadfilename="results_temp.jld2"
                cp(outputfilename,loadfilename;force=true);
                @load loadfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out
            end

            gamma_plot=Vector{Float64}(undef, N);
            deltap=maximum(p_new)-minimum(p_new);      
            for ind=1:N;
                if gamma_out[ind]>0.8;
                    gamma_plot[ind]=1;
                else
                    gamma_plot[ind]=0;
                end
            end
            deltagamma=maximum(gamma_plot)-minimum(gamma_plot);

            #for poly plot
            inds0=findall(gamma_out.>-0.5);
            N0=length(inds0);
            X=Array{Float64}(undef, 3, N0);
            Y=Array{Float64}(undef, 3, N0);
            Z=Array{Float64}(undef, 3, N0);
            C_p=Array{Float32}(undef, 3, N0);        
            C_gamma=Array{Float32}(undef, 3, N0);
            inds1=findall(gamma_out.<-0.5);
            N1=length(inds1);
            X1=Array{Float64}(undef, 3, N1);
            Y1=Array{Float64}(undef, 3, N1);
            Z1=Array{Float64}(undef, 3, N1);  
            C1_p=Array{Float32}(undef, 3, N1);
            C1_gamma=Array{Float32}(undef, 3, N1);
            for i in 1:N0;
                ind=inds0[i];
                X[1,i]=gridx[cellgridid[ind,1]];
                X[2,i]=gridx[cellgridid[ind,2]];
                X[3,i]=gridx[cellgridid[ind,3]];
                Y[1,i]=gridy[cellgridid[ind,1]];
                Y[2,i]=gridy[cellgridid[ind,2]];
                Y[3,i]=gridy[cellgridid[ind,3]];
                Z[1,i]=gridz[cellgridid[ind,1]];
                Z[2,i]=gridz[cellgridid[ind,2]];
                Z[3,i]=gridz[cellgridid[ind,3]];
                C_gamma[1,i]=gamma_plot[ind]/deltagamma;
                C_gamma[2,i]=gamma_plot[ind]/deltagamma;
                C_gamma[3,i]=gamma_plot[ind]/deltagamma;
                C_p[1,i]=p_new[ind]/deltap;
                C_p[2,i]=p_new[ind]/deltap;
                C_p[3,i]=p_new[ind]/deltap;
            end
            xyz = reshape([X[:] Y[:] Z[:]]', :)
            for i in 1:N1
                ind=inds1[i];
                X1[1,i]=gridx[cellgridid[ind,1]];
                X1[2,i]=gridx[cellgridid[ind,2]];
                X1[3,i]=gridx[cellgridid[ind,3]];
                Y1[1,i]=gridy[cellgridid[ind,1]];
                Y1[2,i]=gridy[cellgridid[ind,2]];
                Y1[3,i]=gridy[cellgridid[ind,3]];
                Z1[1,i]=gridz[cellgridid[ind,1]];
                Z1[2,i]=gridz[cellgridid[ind,2]];
                Z1[3,i]=gridz[cellgridid[ind,3]];
                C1_gamma[1,i]=0.5;
                C1_gamma[2,i]=0.5;
                C1_gamma[3,i]=0.5;
                C1_p[1,i]=p_new[ind]/deltap;
                C1_p[2,i]=p_new[ind]/deltap;
                C1_p[3,i]=p_new[ind]/deltap;
            end
            xyz1 = reshape([X1[:] Y1[:] Z1[:]]', :)

            #bounding box
            deltax=maximum(gridx)-minimum(gridx);
            deltay=maximum(gridy)-minimum(gridy);
            deltaz=maximum(gridz)-minimum(gridz);
            mindelta=min(deltax,deltay,deltaz);
            maxdelta=max(deltax,deltay,deltaz);
            if mindelta<maxdelta*0.001;
                eps_delta=maxdelta*0.001;
            else
                eps_delta=0;
            end 
            ax=(deltax+eps_delta)/(mindelta+eps_delta);
            ay=(deltay+eps_delta)/(mindelta+eps_delta);
            az=(deltaz+eps_delta)/(mindelta+eps_delta);
            if i_plot==1;
                ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
                poly!(connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
                if N1>0;
                    poly!(connect(xyz1, Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
                end
                hidedecorations!(ax1);
                hidespines!(ax1) 
            elseif i_plot==2
                ax2 = Axis3(fig[1, 2]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
                poly!(connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
                if N1>0;
                    poly!(connect(xyz1, Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
                end
                hidedecorations!(ax2);
                hidespines!(ax2) 
            elseif i_plot==3
                ax3 = Axis3(fig[1, 3]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
                poly!(connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
                if N1>0;
                    poly!(connect(xyz1, Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
                end
                hidedecorations!(ax3);
                hidespines!(ax3) 
            elseif i_plot==4
                ax4 = Axis3(fig[1, 4]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
                poly!(connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
                if N1>0;
                    poly!(connect(xyz1, Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
                end
                hidedecorations!(ax4);
                hidespines!(ax4) 
            end
            i_out=i_out+Int64(n_pics/4);
        end
        display(fig)
    end 


    function plot_filling(n_out,n_pics)
        #create a window showing the filling factor contour plot at a selected time instance. Selection is with slider bar.
        #n_out ist the index of the last output file, if n_out==-1 the output file with the highest index is chosen
        #consider the last n_pics for creating the contour plots at four equidistant time intervals, if n_pics==-1 all available output files are considered
        #default call is plot_filling(-1,-1)

        val=0;
		n_out_start=-1;
        if n_out==-1;
            vec1=glob("output_*.jld2");
            for i=1:length(vec1);
                vec2=split(vec1[i],".")
                vec3=split(vec2[1],"_")
                val=max(val,parse(Int64,vec3[2]))
				if i==1;
				    n_out_start=parse(Int64,vec3[2]);
				end
            end            
            n_out=val;
        end
		if n_pics==-1;
		    n_pics=(n_out-n_out_start);
		end
        #plot the last n_pics pictures
        if n_pics<4;
            errorstring=string("Makes no sense for n_pics<4"* "\n"); 
            error(errorstring);
        end
        t_digits=2; 
        t_div=10^2;
        
        time_vector=[];
        output_array=[];
        inds=[];
        inds0=[];
        inds1=[];
        N=Int64(0);
        N0=Int64(0);
        N1=Int64(0);
        ax=Float64(0.0);
        ay=Float64(0.0);
        az=Float64(0.0);
        t=Float64(0);
        xyz=Vector{Float64};
        xyz1=Vector{Float64};
        X=Vector{Float64};
        Y=Vector{Float64};
        Z=Vector{Float64};
        C=Vector{Float64};
        X1=Vector{Float64};
        Y1=Vector{Float64};
        Z1=Vector{Float64};
        C1_gamma=Vector{Float64};        
        i_out=n_out-n_pics;
        i_firstfile=1;
        for i_plot in 1:n_pics+1;          
            outputfilename=string("output_",string(i_out),".jld2");
            if ~isfile(outputfilename);
                errorstring=string("File ",outputfilename," not existing"* "\n"); 
                error(errorstring);
            else
                loadfilename="results_temp.jld2"
                cp(outputfilename,loadfilename;force=true);
                @load loadfilename t rho_new u_new v_new p_new gamma_new gamma_out gridx gridy gridz cellgridid N n_out
                #print(string(i_plot)*" \n")
                #print(string(i_out)*" \n")
                if i_firstfile==1;
                    i_firstfile=0;
                    #for poly plot
                    inds0=findall(gamma_out.>-0.5);
                    N0=length(inds0);
                    X=Array{Float64}(undef, 3, N0);
                    Y=Array{Float64}(undef, 3, N0);
                    Z=Array{Float64}(undef, 3, N0);
                    C_p=Array{Float32}(undef, 3, N0);        
                    C_gamma=Array{Float32}(undef, 3, N0);
                    inds1=findall(gamma_out.<-0.5);
                    N1=length(inds1);
                    X1=Array{Float64}(undef, 3, N1);
                    Y1=Array{Float64}(undef, 3, N1);
                    Z1=Array{Float64}(undef, 3, N1);  
                    C1_gamma=Array{Float32}(undef, 3, N1);
                    for i in 1:N0;
                        ind=inds0[i];
                        X[1,i]=gridx[cellgridid[ind,1]];
                        X[2,i]=gridx[cellgridid[ind,2]];
                        X[3,i]=gridx[cellgridid[ind,3]];
                        Y[1,i]=gridy[cellgridid[ind,1]];
                        Y[2,i]=gridy[cellgridid[ind,2]];
                        Y[3,i]=gridy[cellgridid[ind,3]];
                        Z[1,i]=gridz[cellgridid[ind,1]];
                        Z[2,i]=gridz[cellgridid[ind,2]];
                        Z[3,i]=gridz[cellgridid[ind,3]];
                    end
                    xyz = reshape([X[:] Y[:] Z[:]]', :)
                    for i in 1:N1
                        ind=inds1[i];
                        X1[1,i]=gridx[cellgridid[ind,1]];
                        X1[2,i]=gridx[cellgridid[ind,2]];
                        X1[3,i]=gridx[cellgridid[ind,3]];
                        Y1[1,i]=gridy[cellgridid[ind,1]];
                        Y1[2,i]=gridy[cellgridid[ind,2]];
                        Y1[3,i]=gridy[cellgridid[ind,3]];
                        Z1[1,i]=gridz[cellgridid[ind,1]];
                        Z1[2,i]=gridz[cellgridid[ind,2]];
                        Z1[3,i]=gridz[cellgridid[ind,3]];
                        C1_gamma[1,i]=0.5;
                        C1_gamma[2,i]=0.5;
                        C1_gamma[3,i]=0.5;
                    end
                    xyz1 = reshape([X1[:] Y1[:] Z1[:]]', :)

                    #bounding box
                    deltax=maximum(gridx)-minimum(gridx);
                    deltay=maximum(gridy)-minimum(gridy);
                    deltaz=maximum(gridz)-minimum(gridz);
                    mindelta=min(deltax,deltay,deltaz);
                    maxdelta=max(deltax,deltay,deltaz);
                    if mindelta<maxdelta*0.001;
                        eps_delta=maxdelta*0.001;
                    else
                        eps_delta=0;
                    end 
                    ax=(deltax+eps_delta)/(mindelta+eps_delta);
                    ay=(deltay+eps_delta)/(mindelta+eps_delta);
                    az=(deltaz+eps_delta)/(mindelta+eps_delta);
                    time_vector=t;
                    output_array=gamma_out; 
                    N_val=N;

                else
                    time_vector=vcat(time_vector,t);
                    output_array=hcat(output_array,gamma_out);
                end
            end
            i_out=i_out+1;
        end

        gamma_plot=output_array[:,end]  
        for ind=1:N;
            if gamma_plot[ind]>0.8;
                gamma_plot[ind]=1;
            else
                gamma_plot[ind]=0;
            end
        end
        deltagamma=maximum(gamma_plot)-minimum(gamma_plot);
        
        C_gamma=Array{Float32}(undef, 3, N0);
        for i in 1:N0;
            ind=inds0[i];
            C_gamma[1,i]=gamma_plot[ind]/deltagamma;
            C_gamma[2,i]=gamma_plot[ind]/deltagamma;
            C_gamma[3,i]=gamma_plot[ind]/deltagamma;
        end

        resolution_val=600;
        fig = Figure(resolution=(resolution_val, resolution_val))   
        ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title=string("Filling factor at t=", string(round(t_div*t)/t_div) ,"s"))
        #p1=poly!(ax1,connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
        #if N1>0;
        #    p2=poly!(ax1,connect(xyz1, Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
        #end
        hidedecorations!(ax1);
        hidespines!(ax1) 
        #sl_t = Slider(fig[2, 1], range = time_vector[1]:  (time_vector[end]-time_vector[1])/n_pics :time_vector[end], startvalue =  time_vector[end] );
        sl_t = Slider(fig[2, 1], range = time_vector[1]:  (time_vector[end]-time_vector[1])/n_pics :time_vector[end], startvalue =  time_vector[1] );
        point = lift(sl_t.value) do x           
            if x<0.5*(time_vector[end]+time_vector[1]);
                gamma_plot=output_array[:,1]
            else
                gamma_plot=output_array[:,end]
            end
            time_val=x
            tind=1
            tind1=1
            tind2=2
            for i in 1:length(time_vector)-1;
                if x>=0.5*(time_vector[i]+time_vector[i+1])
                    tind=i+1;
                end
            end            
            gamma_plot=output_array[:,tind]
            time_val=time_vector[tind]

            for ind=1:N;
                if gamma_plot[ind]>0.8;
                    gamma_plot[ind]=1;
                else
                    gamma_plot[ind]=0;
                end
            end
            deltagamma=maximum(gamma_plot)-minimum(gamma_plot);
            for i in 1:N0;
                ind=inds0[i];
                C_gamma[1,i]=gamma_plot[ind]/deltagamma;
                C_gamma[2,i]=gamma_plot[ind]/deltagamma;
                C_gamma[3,i]=gamma_plot[ind]/deltagamma;
            end
            #empty!(ax1.scene)
            p1=poly!(ax1,connect(xyz, Point{3}), connect(1:length(X), TriangleFace); color=C_gamma[:], strokewidth=1, colorrange=(0,1))
            if N1>0;
                p2=poly!(ax1,connect(xyz1, Point{3}), connect(1:length(X1), TriangleFace); color=C1_gamma[:], strokewidth=1, colorrange=(0,1),colormap = (:bone))
            end
            ax1.title=string("Filling factor at t=", string(round(t_div*time_val)/t_div) ,"s")
            hidedecorations!(ax1);
            hidespines!(ax1) 
            display(fig)
        end
    end 

end
