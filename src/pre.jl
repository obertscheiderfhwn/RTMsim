"""
    start_rtmsim(inputfilename)
    
Reads the text input file and calls the solver with the read parameters. 

Arguments:
- inputfilename :: String

The complete set of input parameters can be accessed in the input file. The following paragraph shows an example for such an input file:
```
    1    #i_model 
    meshfiles/mesh_permeameter1_foursets.bdf    #meshfilename 
    200    #tmax 
    1.01325e5 1.225 1.4 0.06    #p_ref rho_ref gamma mu_resin_val 
    1.35e5 1.0e5    #p_a_val p_init_val 
    3e-3 0.7 3e-10 1 1 0 0    #t_val porosity_val K_val alpha_val refdir1_val refdir2_val refdir3_val 
    3e-3 0.7 3e-10 1 1 0 0    #t1_val porosity1_val K1_val alpha1_val refdir11_val refdir21_val refdir31_val 
    3e-3 0.7 3e-10 1 1 0 0    #t2_val porosity2_val K2_val alpha2_val refdir12_val refdir22_val refdir32_val
    3e-3 0.7 3e-10 1 1 0 0    #t3_val porosity3_val K3_val alpha3_val refdir13_val refdir23_val refdir33_val
    3e-3 0.7 3e-10 1 1 0 0    #t4_val porosity4_val K4_val alpha4_val refdir14_val refdir24_val refdir34_val 
    1 0 0 0    #patchtype1val patchtype2val patchtype3val patchtype4val 
    0 results.jld2    #i_restart restartfilename
    0 0.01    #i_interactive r_p
    16    #n_pics
```

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
Entries are separated by one blank.

Unit test:
- `MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]; inputfilename=joinpath(MODULE_ROOT,"inputfiles","input.txt"); rtmsim.start_rtmsim(inputfilename);`
"""
function start_rtmsim(inputfilename)     
    if Sys.iswindows()
        inputfilename=replace(inputfilename,"/" => "\\")
    elseif Sys.islinux()
        inputfilename=replace(inputfilename,"\\" => "/")
    end        
    println("Read input file "*string(inputfilename))
    if ~isfile(inputfilename)
        errorstring=string("File ",inputfilename," not existing"* "\n") 
        error(errorstring)
    end


    A_input=readdlm(inputfilename,' ')
    param=rtmsim.input_vals(A_input[1,1],
                            A_input[2,1],
                            A_input[3,1],
                            A_input[4,1],A_input[4,2],A_input[4,3],A_input[4,4],
                            A_input[5,1],A_input[5,2],
                            A_input[6,1],A_input[6,2],A_input[6,3],A_input[6,4],A_input[6,5],A_input[6,6],A_input[6,7],
                            A_input[7,1],A_input[7,2],A_input[7,3],A_input[7,4],A_input[7,5],A_input[7,6],A_input[7,7],
                            A_input[8,1],A_input[8,2],A_input[8,3],A_input[8,4],A_input[8,5],A_input[8,6],A_input[8,7],
                            A_input[9,1],A_input[9,2],A_input[9,3],A_input[9,4],A_input[9,5],A_input[9,6],A_input[9,7],
                            A_input[10,1],A_input[10,2],A_input[10,3],A_input[10,4],A_input[10,5],A_input[10,6],A_input[10,7],
                            A_input[11,1],A_input[11,2],A_input[11,3],A_input[11,4],
                            A_input[12,1],A_input[12,2],
                            A_input[13,1],A_input[13,2],
                            A_input[14,1])
    #println("param=",param)

    println(" ")
    rtmsim_rev1(param)
end


"""
    function read_mesh(input_struct)

Read mesh file and prepare to be used in solver:
- number of cells, cell ids start with 1
- x,y,z-coordinates of the nodes
- x,y,z-coordinates of the cell centers
- patch properties
Read other mesh files than Nastran bulk data format (bdf) based on extension and calculate the required mesh data or convert to Nastran format prepare with existing function                
"""
function read_mesh(input_struct)
    meshfilename=input_struct.meshfilename
    #read Nastran mesh
    if meshfilename[end-2:end]=="bdf"
        return_struct=read_nastran_mesh(input_struct)      
    end
    return return_struct
end


"""
    function read_nastran_mesh(input_struct)
        
Read file in Nastran format with fixed length (8 digits), nodes (`GRIDS`) defined in global coordinate system.

Arguments:
- meshfilename :: String
- paramset, paramset1, paramset2, paramset3, paramset3 :: Vector{Float}
- patchtype1val,patchtype1val1,patchtype1val2,patchtype1val3,patchtype1val4 :: Int
- i_interactive :: Int
- r_p :: Float

Unit test:
- `MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]; meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf"); paramset=[0.5,0.3,3e-10,1.0,1.0,0.0,0.0];paramset1=paramset;paramset2=paramset;paramset3=paramset;paramset4=paramset;patchtype1val=-1;patchtype2val=-1;patchtype3val=-1;patchtype4val=-1;i_interactive=0;r_p=0.01; input_struct=rtmsim.input_args_read_mesh(meshfilename,paramset,paramset1,paramset2,paramset3,paramset4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_interactive,r_p); return_struct=rtmsim.read_mesh(input_struct);`
"""
function read_nastran_mesh(input_struct)
    meshfilename=input_struct.meshfilename
    paramset=input_struct.paramset
    paramset1=input_struct.paramset1
    paramset2=input_struct.paramset2
    paramset3=input_struct.paramset3
    paramset4=input_struct.paramset4
    patchtype1val=input_struct.patchtype1val
    patchtype2val=input_struct.patchtype2val
    patchtype3val=input_struct.patchtype3val
    patchtype4val=input_struct.patchtype4val
    i_interactive=input_struct.i_interactive
    r_p=input_struct.r_p
    if ~isfile(meshfilename)
        errorstring=string("File ",meshfilename," not existing"* "\n") 
        error(errorstring)
    end
    ind=Int64(1)
    gridind=Int64(1)
    setind=Int64(1)
    issetdefinition=Int64(0)
    patchorigids1=[]
    patchorigids2=[]
    patchorigids3=[]
    patchorigids4=[]
    origgridid=[]
    gridx=[]
    gridy=[]
    gridz=[]
    celloriggridid=[]
    cellgridid=Array{Int64}(undef, 0, 3)
    inletpatchids=[]

    open(meshfilename, "r") do fid
        line=1
        while !eof(fid)
            thisline=readline(fid)
            if length(thisline)>=8
                if issetdefinition==1 
                    if cmp( thisline[1:8],"        ")!=0  #check if the first eight characters are empty, else issetdefinition=0
                        issetdefinition=Int64(0)
                        setind=setind+1
                    end
                end
                card=thisline[1:8]
                if cmp(card,"GRID    ")==0
                    gridindstring=thisline[9:16]
                    origgridid=vcat(origgridid,parse(Int64,gridindstring))                        
                    txt=thisline[25:32]
                    txt=replace(txt," "=> "");txt=replace(txt,"E" => "");txt=replace(txt,"e" => "")
                    txt1=replace(txt,"-" => "e-");txt1=replace(txt1,"+" => "e+")
                    if cmp(txt1[1],'e')==0;txt2=txt1[2:end];else;txt2=txt1;end
                    val=parse(Float64,txt2)
                    val1=val
                    txt=thisline[33:40]
                    txt=replace(txt," "=> "");txt=replace(txt,"E" => "");txt=replace(txt,"e" => "")
                    txt1=replace(txt,"-" => "e-");txt1=replace(txt1,"+" => "e+")
                    if cmp(txt1[1],'e')==0;txt2=txt1[2:end];else;txt2=txt1;end
                    val=parse(Float64,txt2)
                    val2=val
                    txt=thisline[41:48]
                    txt=replace(txt," "=> "");txt=replace(txt,"E" => "");txt=replace(txt,"e" => "")
                    txt1=replace(txt,"-" => "e-");txt1=replace(txt1,"+" => "e+")
                    if cmp(txt1[1],'e')==0;txt2=txt1[2:end];else;txt2=txt1;end
                    val=parse(Float64,txt2)
                    val3=val
                    gridx=vcat(gridx,val1)
                    gridy=vcat(gridy,val2)
                    gridz=vcat(gridz,val3)
                    gridind=gridind+1
                elseif cmp(card,"CTRIA3  ")==0        
                    celloriggridid=vcat(celloriggridid,parse(Int64,thisline[9:16]))
                    i1val=parse(Int64,thisline[25:32])
                    i1=findfirst(isequal(i1val),origgridid)
                    i2val=parse(Int64,thisline[33:40])
                    i2=findfirst(isequal(i2val),origgridid)
                    i3val=parse(Int64,thisline[41:48])
                    i3=findfirst(isequal(i3val),origgridid)
                    ivec=[i1,i2,i3]
                    idel=findall(isequal(min(ivec[1],ivec[2],ivec[3])),ivec)
                    deleteat!(ivec,idel)
                    idel=findall(isequal(max(ivec[1],ivec[2])),ivec)
                    deleteat!(ivec,idel)                        
                    cellgridid=vcat(cellgridid,[min(i1,i2,i3) ivec[1] max(i1,i2,i3)])
                    ind=ind+1    
                elseif cmp( card[1:3],"SET")==0 || issetdefinition==1
                    issetdefinition=1
                    txt1=thisline[9:end]
                    txt1=replace(txt1," "=> "")
                    txt2=split(txt1,",")
                    for i in 1:length(txt2)
                        if !isempty(txt2[i])
                            if setind==1 
                                patchorigids1=vcat(patchorigids1,parse(Int64,txt2[i]))
                            elseif setind==2
                                patchorigids2=vcat(patchorigids2,parse(Int64,txt2[i]))
                            elseif setind==3
                                patchorigids3=vcat(patchorigids3,parse(Int64,txt2[i]))
                            elseif setind==4
                                patchorigids4=vcat(patchorigids4,parse(Int64,txt2[i]))
                            end
                        end
                    end
                end
            end
            line+=1
        end
    end
    N=ind-1  #total number of cells

    #loop to define cell center coordinates in global CS
    cellcenterx=[]
    cellcentery=[]
    cellcenterz=[]
    for ind in 1:N
        i1=cellgridid[ind,1]
        i2=cellgridid[ind,2]
        i3=cellgridid[ind,3]
        cellcenterx=vcat(cellcenterx,(gridx[i1]+gridx[i2]+gridx[i3])/3)
        cellcentery=vcat(cellcentery,(gridy[i1]+gridy[i2]+gridy[i3])/3)
        cellcenterz=vcat(cellcenterz,(gridz[i1]+gridz[i2]+gridz[i3])/3)
    end

    if i_interactive==1
        input_pset=rtmsim.input_args_assign_pset(r_p,N,cellcenterx,cellcentery,cellcenterz)
        assign_pset(input_pset)
        psetfilename="pset.jld2"
        if ~isfile(psetfilename)
            errorstring=string("File ",psetfilename," not existing"* "\n") 
            error(errorstring)
        end
        @load psetfilename pset
        inletpatchids=pset
        if length(inletpatchids)<1
            errorstring=string("Inlet definition empty"* "\n") 
            error(errorstring)
        end
        patchids1=[]
        patchids2=[]
        patchids3=[]
        patchids4=[]   
        patchparameters=paramset
        patchparameters1=[]
        patchparameters2=[]
        patchparameters3=[]
        patchparameters4=[]        
    else
        patchids1=[]
        patchids2=[]
        patchids3=[]
        patchids4=[]
        for i in 1:length(patchorigids1)
            i1=findfirst(isequal(patchorigids1[i]),celloriggridid)
            patchids1=vcat(patchids1,i1)
        end
        for i=1:length(patchorigids2)
            i1=findfirst(isequal(patchorigids2[i]),celloriggridid)
            patchids2=vcat(patchids2,i1)
        end
        for i=1:length(patchorigids3)
            i1=findfirst(isequal(patchorigids3[i]),celloriggridid)
            patchids3=vcat(patchids3,i1)
        end
        for i=1:length(patchorigids4)
            i1=findfirst(isequal(patchorigids4[i]),celloriggridid)
            patchids4=vcat(patchids4,i1)
        end
        if i_interactive==2
            input_pset=rtmsim.input_args_assign_pset(r_p,N,cellcenterx,cellcentery,cellcenterz)
            assign_pset(input_pset)
            psetfilename="pset.jld2"
            if ~isfile(psetfilename)
                errorstring=string("File ",psetfilename," not existing"* "\n") 
                error(errorstring)
            end
            @load psetfilename pset
            inletpatchids=pset
            if length(patchids1)<1
                errorstring=string("Inlet definition empty"* "\n") 
                error(errorstring)
            end
        end
        patchparameters=paramset
        patchparameters1=[]
        patchparameters2=[]
        patchparameters3=[]
        patchparameters4=[]
        for i_patch in 1:4
            if i_patch==1
                patchids=patchids1
            elseif i_patch==2
                patchids=patchids2
            elseif i_patch==3
                patchids=patchids3
            elseif i_patch==4
                patchids=patchids4
            end
            if !isempty(patchids)
                if i_patch==1
                    if patchtype1val==2
                        patchparameters1=paramset1
                    end
                elseif i_patch==2
                    if patchtype2val==2
                        patchparameters2=paramset2
                    end
                elseif i_patch==3
                    if patchtype3val==2 
                        patchparameters3=paramset3
                    end
                elseif i_patch==4
                    if patchtype4val==2
                        patchparameters4=paramset4
                    end
                end
            end
        end
    end
    
    return_struct=rtmsim.return_args_read_mesh(N,cellgridid,gridx,gridy,gridz,cellcenterx,cellcentery,cellcenterz,patchparameters,patchparameters1,patchparameters2,patchparameters3,patchparameters4,patchids1,patchids2,patchids3,patchids4,inletpatchids)
    return return_struct
end


"""
    function plot_mesh(meshfilename,i_mode)

Create mesh plot with cells with `i_mode==1` and create mesh plots with cell center nodes with `i_mode==2` for manual selection of inlet ports.

Arguments:
- meshfilename :: String
- i_mode :: Int

Unit test:
- `MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]; meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf"); rtmsim.plot_mesh(meshfilename,1);`

Additional unit tests:
- `MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]; meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf"); rtmsim.plot_mesh(meshfilename,2);` for the manual selection of inlet ports with left mouse button click while key p is pressed 
- `MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]; meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf"); param=rtmsim.input_vals(1,meshfilename,200, 0.35e5,1.205,1.4,0.06, 0.35e5,0.00e5, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 3e-3,0.7,3e-10,1,1,0,0, 0,0,0,0, 0,"results.jld2",1,0.01,16); rtmsim.rtmsim_rev1(param);` for starting only with the interactively selected inlet ports
"""
function plot_mesh(meshfilename,i_mode)
    if Sys.iswindows()
        meshfilename=replace(meshfilename,"/" => "\\")
    elseif Sys.islinux()
        meshfilename=replace(meshfilename,"\\" => "/")
    end      
    #dummy values for calling function read_mesh
    paramset=[0.5,0.3,3e-10,1.0,1.0,0.0,0.0];paramset1=paramset;paramset2=paramset;paramset3=paramset;paramset4=paramset
    patchtype1val=-1;patchtype2val=-1;patchtype3val=-1;patchtype4val=-1;i_interactive=0
    r_p=0.01
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

    #for poly plot
    X=Array{Float64}(undef, 3, N)
    Y=Array{Float64}(undef, 3, N)
    Z=Array{Float64}(undef, 3, N)
    C=Array{Float32}(undef, 3, N)
    for ind in 1:N
        X[1,ind]=gridx[cellgridid[ind,1]]
        X[2,ind]=gridx[cellgridid[ind,2]]
        X[3,ind]=gridx[cellgridid[ind,3]]
        Y[1,ind]=gridy[cellgridid[ind,1]]
        Y[2,ind]=gridy[cellgridid[ind,2]]
        Y[3,ind]=gridy[cellgridid[ind,3]]
        Z[1,ind]=gridz[cellgridid[ind,1]]
        Z[2,ind]=gridz[cellgridid[ind,2]]
        Z[3,ind]=gridz[cellgridid[ind,3]]
        C[1,ind]=1.0
        C[2,ind]=1.0
        C[3,ind]=1.0
    end
    xyz = reshape([X[:] Y[:] Z[:]]', :)
    #2..for meshscatter plot
    X2=Array{Float64}(undef, 3*N)
    Y2=Array{Float64}(undef, 3*N)
    Z2=Array{Float64}(undef, 3*N)
    C2=Array{Float64}(undef, 3*N)
    for ind in 1:N
        X2[    ind]=gridx[cellgridid[ind,1]]
        X2[  N+ind]=gridx[cellgridid[ind,2]]
        X2[2*N+ind]=gridx[cellgridid[ind,3]]
        Y2[    ind]=gridy[cellgridid[ind,1]]
        Y2[  N+ind]=gridy[cellgridid[ind,2]]
        Y2[2*N+ind]=gridy[cellgridid[ind,3]]
        Z2[    ind]=gridz[cellgridid[ind,1]]
        Z2[  N+ind]=gridz[cellgridid[ind,2]]
        Z2[2*N+ind]=gridz[cellgridid[ind,3]]
        C2[    ind]=0.0
        C2[  N+ind]=0.0
        C2[2*N+ind]=0.0
    end

    #bounding box
    deltax=maximum(gridx)-minimum(gridx)
    deltay=maximum(gridy)-minimum(gridy)
    deltaz=maximum(gridz)-minimum(gridz)
    mindelta=min(deltax,deltay,deltaz)
    maxdelta=max(deltax,deltay,deltaz)
    if mindelta<maxdelta*0.001
        eps_delta=maxdelta*0.001
    else
        eps_delta=0
    end 
    ax=(deltax+eps_delta)/(mindelta+eps_delta)
    ay=(deltay+eps_delta)/(mindelta+eps_delta)
    az=(deltaz+eps_delta)/(mindelta+eps_delta)

    if i_mode==1
        fig = Figure(resolution=(600, 600))
        ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Mesh")
        poly!(connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C[:], strokewidth=1)
        #hidedecorations!(ax1)
        hidespines!(ax1) 
        display(fig)
    elseif i_mode==2
        points=rand(Point3f0, length(gridx))
        for i in 1:length(gridx)
            points[i]=Point3f0(gridx[i],gridy[i],gridz[i])
        end
        positions = Observable(points) 

        inletpos_xyz=[-9.9e9 -9.9e9 -9.9e9]
        filename="inletpostions.jld2"
        @save filename inletpos_xyz

        markersizeval=maxdelta*100
        fig = Figure()
        ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Select inlets with p + LMB")
        p=scatter!(ax1, positions,markersize=markersizeval)
        hidedecorations!(ax1)
        hidespines!(ax1) 

        on(events(fig).mousebutton, priority = 2) do event
            if event.button == Mouse.left && event.action == Mouse.press
                if Keyboard.p in events(fig).keyboardstate
                    plt, i = pick(fig.scene,events(fig).mouseposition[])
                    if plt == p
                        @load filename inletpos_xyz
                        t_div=100
                        xpos=positions[][i][1]
                        ypos=positions[][i][2]
                        zpos=positions[][i][3]
                        inletpos_xyz=vcat(inletpos_xyz,[xpos ypos zpos])
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


"""
    function plot_sets(meshfilename)
    
Create a plot with the up to four cell sets defined in the mesh file.

Unit test:
- `MODULE_ROOT=splitdir(splitdir(pathof(rtmsim))[1])[1]; meshfilename=joinpath(MODULE_ROOT,"meshfiles","mesh_permeameter1_foursets.bdf"); rtmsim.plot_sets(meshfilename);`
"""
function plot_sets(meshfilename)
    #dummy values for calling function read_nastran_mesh
    paramset=[0.5,0.3,3e-10,1.0,1.0,0.0,0.0];paramset1=paramset;paramset2=paramset;paramset3=paramset;paramset4=paramset
    patchtype1val=-1;patchtype2val=-1;patchtype3val=-1;patchtype4val=-1;i_interactive=0
    r_p=0.01
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

    if isempty(patchids1)
        n_patch=0
        errorstring=string("No sets defined"* "\n") 
        error(errorstring)
    else
        if isempty(patchids2)
            n_patch=1
        else
            if isempty(patchids3)
                n_patch=2
            else
                if isempty(patchids4)
                    n_patch=3
                else
                    n_patch=4
                end
            end
        end
    end

    #for poly plot
    X=Array{Float64}(undef, 3, N)
    Y=Array{Float64}(undef, 3, N)
    Z=Array{Float64}(undef, 3, N)
    C=Array{Float32}(undef, 3, N)
    C_patch1=Array{Float32}(undef, 3, N)
    C_patch2=Array{Float32}(undef, 3, N)
    C_patch3=Array{Float32}(undef, 3, N)
    C_patch4=Array{Float32}(undef, 3, N)
    for ind in 1:N
        X[1,ind]=gridx[cellgridid[ind,1]]
        X[2,ind]=gridx[cellgridid[ind,2]]
        X[3,ind]=gridx[cellgridid[ind,3]]
        Y[1,ind]=gridy[cellgridid[ind,1]]
        Y[2,ind]=gridy[cellgridid[ind,2]]
        Y[3,ind]=gridy[cellgridid[ind,3]]
        Z[1,ind]=gridz[cellgridid[ind,1]]
        Z[2,ind]=gridz[cellgridid[ind,2]]
        Z[3,ind]=gridz[cellgridid[ind,3]]
        C[1,ind]=1.0
        C[2,ind]=1.0
        C[3,ind]=1.0
        if issubset(ind, patchids1)
            C_patch1[1,ind]=1.0
            C_patch1[2,ind]=1.0
            C_patch1[3,ind]=1.0
        else
            C_patch1[1,ind]=0.0
            C_patch1[2,ind]=0.0
            C_patch1[3,ind]=0.0
        end   
        if issubset(ind, patchids2)
            C_patch2[1,ind]=1.0
            C_patch2[2,ind]=1.0
            C_patch2[3,ind]=1.0
        else
            C_patch2[1,ind]=0.0
            C_patch2[2,ind]=0.0
            C_patch2[3,ind]=0.0
        end         
        if issubset(ind, patchids3)
            C_patch3[1,ind]=1.0
            C_patch3[2,ind]=1.0
            C_patch3[3,ind]=1.0
        else
            C_patch3[1,ind]=0.0
            C_patch3[2,ind]=0.0
            C_patch3[3,ind]=0.0
        end  
        if issubset(ind, patchids4)
            C_patch4[1,ind]=1.0
            C_patch4[2,ind]=1.0
            C_patch4[3,ind]=1.0
        else
            C_patch4[1,ind]=0.0
            C_patch4[2,ind]=0.0
            C_patch4[3,ind]=0.0
        end        
    end
    xyz = reshape([X[:] Y[:] Z[:]]', :)

    #2..for meshscatter plot
    X2=Array{Float64}(undef, 3*N)
    Y2=Array{Float64}(undef, 3*N)
    Z2=Array{Float64}(undef, 3*N)
    C2=Array{Float64}(undef, 3*N)
    for ind in 1:N
        X2[    ind]=gridx[cellgridid[ind,1]]
        X2[  N+ind]=gridx[cellgridid[ind,2]]
        X2[2*N+ind]=gridx[cellgridid[ind,3]]
        Y2[    ind]=gridy[cellgridid[ind,1]]
        Y2[  N+ind]=gridy[cellgridid[ind,2]]
        Y2[2*N+ind]=gridy[cellgridid[ind,3]]
        Z2[    ind]=gridz[cellgridid[ind,1]]
        Z2[  N+ind]=gridz[cellgridid[ind,2]]
        Z2[2*N+ind]=gridz[cellgridid[ind,3]]
        C2[    ind]=0.0
        C2[  N+ind]=0.0
        C2[2*N+ind]=0.0
    end

    resolution_val=300
    if n_patch==1
        fig = Figure(resolution=(1*resolution_val, resolution_val))
    elseif n_patch==2
        fig = Figure(resolution=(2*resolution_val, 2*resolution_val))
    elseif n_patch==3
        fig = Figure(resolution=(3*resolution_val, resolution_val))
    elseif n_patch==4
        fig = Figure(resolution=(4*resolution_val, resolution_val))
    end

    #bounding box
    deltax=maximum(gridx)-minimum(gridx)
    deltay=maximum(gridy)-minimum(gridy)
    deltaz=maximum(gridz)-minimum(gridz)
    mindelta=min(deltax,deltay,deltaz)
    maxdelta=max(deltax,deltay,deltaz) 
    if mindelta<maxdelta*0.001
        eps_delta=maxdelta*0.001
    else
        eps_delta=0
    end 
    ax=(deltax+eps_delta)/(mindelta+eps_delta)
    ay=(deltay+eps_delta)/(mindelta+eps_delta)
    az=(deltaz+eps_delta)/(mindelta+eps_delta)
    ax1 = Axis3(fig[1, 1]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Set 1")
    poly!(connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_patch1[:], strokewidth=1)
    hidedecorations!(ax1);hidespines!(ax1) 
    if n_patch>=2
        ax2 = Axis3(fig[1, 2]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Set 2")
        poly!(connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_patch2[:], strokewidth=1)
        hidedecorations!(ax2);hidespines!(ax2) 
    end
    if n_patch>=3
        ax3 = Axis3(fig[1, 3]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Set 3")
        poly!(connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_patch3[:], strokewidth=1)
        hidedecorations!(ax3);hidespines!(ax3) 
    end
    if n_patch>=4
        ax4 = Axis3(fig[1, 4]; aspect=(ax,ay,az), perspectiveness=0.5,viewmode = :fitzoom,title="Set 4")
        poly!(connect(xyz, Makie.Point{3}), connect(1:length(X), TriangleFace); color=C_patch4[:], strokewidth=1)
        hidedecorations!(ax4);hidespines!(ax4) 
    end
    display(fig)
end


"""
    function assign_pset(r_p,N,cellcenterx,cellcentery,cellcenterz)

Create the cell set from the manually selected inlet port nodes
"""
function assign_pset(input_struct)
    r_p=input_struct.r_p
    N=input_struct.N
    cellcenterx=input_struct.cellcenterx
    cellcentery=input_struct.cellcentery
    cellcenterz=input_struct.cellcenterz
    filename="inletpostions.jld2"
    @load filename inletpos_xyz
    n_p=size(inletpos_xyz,1)-1
    patchpids=[]
    i=1
    for i_p in 2:n_p+1
        r_p_temp=r_p
        i_add=1
        for ind in 1:N
            vec1=[cellcenterx[ind]-inletpos_xyz[i_p,1],cellcentery[ind]-inletpos_xyz[i_p,2],cellcenterz[ind]-inletpos_xyz[i_p,3]]
            if sqrt(dot(vec1,vec1))<=r_p_temp
            patchpids=vcat(patchpids,ind)
            i=i+1
            i_add=0
            end
        end
        while i_add==1
            r_p_temp=1.1*r_p_temp
            i_firstcell=0
            for ind in 1:N
                vec1=[cellcenterx[ind]-inletpos_xyz[i_p,1],cellcentery[ind]-inletpos_xyz[i_p,2],cellcenterz[ind]-inletpos_xyz[i_p,3]]
                if sqrt(dot(vec1,vec1))<=r_p_temp  && i_firstcell==0
                    patchpids=vcat(patchpids,ind)
                    i_firstcell=1
                    i=i+1
                    i_add=0
                end
            end
        end
    end
    pset=patchpids
    psetfilename="pset.jld2"
    @save psetfilename pset
end
