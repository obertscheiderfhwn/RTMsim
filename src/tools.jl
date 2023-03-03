"""
    function create_faces(input_struct)

Find the set with the IDs of the neighbouring cells and identify wall cells.

Arguments: Data structure with
- `cellgridid :: Array{Int,2}`
- `N :: Int`
- `maxnumberofneighbours :: Int`

Meaning of the arguments:
- `cellgridid`: The i-th line contains the three IDs of the nodes which form the cell
- `N`: Number of cells.
- `maxnumberofneighbours`: Number of columns of array `cellneighboursarray`. Default value is `10`. If more cell neighbours in the mesh, an error occurs and this value must be increased. 

Return: Data structure with
- `faces :: Array{Int,2}`
- `cellneighboursarray :: Array{Int,2}`
- `celltype :: Vector{Int}`
"""
function create_faces(input_struct)
    cellgridid=input_struct.cellgridid
    N=input_struct.N
    maxnumberofneighbours=input_struct.maxnumberofneighbours
    celltype=Vector{Int64}(undef, N)
    for i in 1:N
        celltype[i]=1
    end
    faces=Array{Int64}(undef, 0, 3)   #three columns: grid id1, grid id2, cell id
    i=1
    for ind=1:N
        i1=cellgridid[ind,1]
        i2=cellgridid[ind,2]    
        faces=vcat(faces,[min(i1,i2) max(i1,i2) ind])
        i=i+1
        i1=cellgridid[ind,2]
        i2=cellgridid[ind,3]    
        faces=vcat(faces,[min(i1,i2) max(i1,i2) ind])
        i=i+1
        i1=cellgridid[ind,3]
        i2=cellgridid[ind,1]    
        faces=vcat(faces,[min(i1,i2) max(i1,i2) ind])
        i=i+1
    end
    facessorted=sortslices(faces,dims=1)
    vals1=unique(facessorted[:,1])  

    # this must be generalized, currently only hard-coded number of neighbouring cells of a tria is possible
    # all considered cases had <<10 neighbouring cells 
    cellneighboursarray=Array{Int64}(undef, N, maxnumberofneighbours)
    for ind in 1:N
        for ind_n in 1:maxnumberofneighbours
            cellneighboursarray[ind,ind_n]=-9
        end
    end

    for i in 1:length(vals1)
        inds2=findall(isequal(vals1[i]), facessorted[:,1])
        facesdetail_unsorted=facessorted[inds2,2:3]
        facesdetail=sortslices(facesdetail_unsorted,dims=1)
        for j=1:size(facesdetail,1)
            i1=facesdetail[j,2]
            inds3=findall(isequal(facesdetail[j,1]),facesdetail[:,1])
            inds4=findall(!isequal(j),inds3)
            inds5=inds3[inds4]
            if isempty(inds5)
                celltype[i1]=-3  #wall
            else
                if j==1
                    for k in 1:length(inds5)
                        matrixrow=cellneighboursarray[i1,:]
                        indcolumn=findfirst(isequal(-9),matrixrow) 
                        if isnothing(indcolumn)
                            error("More than 10 neighbours of one tria is not supported \n")
                        else
                            cellneighboursarray[i1,indcolumn]=facesdetail[inds5[k],2]
                        end
                    end
                else
                for k in 1:1 
                        matrixrow=cellneighboursarray[i1,:]
                        indcolumn=findfirst(isequal(-9),matrixrow) 
                        if isnothing(indcolumn)
                            error("More than 10 neighbours of one tria is not supported"* "\n")
                        else
                            cellneighboursarray[i1,indcolumn]=facesdetail[inds5[k],2]
                        end
                    end
                end
            end
        end
    end

    return_struct=rtmsim.return_args_create_faces(faces, cellneighboursarray, celltype)
    return return_struct
end


"""
function delete_files()

Deletes the intermediate jld2 output files.
"""
function delete_files()
    #delete the intermediate output files
    rm.(glob("output_*.jld2"))
end


"""
    function assign_parameters(input_struct)

Assign properties to cells.

Arguments: Data structure with
- `i_interactive :: Int`
- `celltype :: Vector{Int}`
- `patchparameters0 :: Vector{Float}`
- `patchparameters1 :: Vector{Float}`
- `patchparameters2 :: Vector{Float}`
- `patchparameters3 :: Vector{Float}`
- `patchparameters4 :: Vector{Float}`
- `patchtype1val :: Int`
- `patchtype2val :: Int`
- `patchtype3val :: Int`
- `patchtype4val :: Int`
- `patchids1 :: Vector{Int}`
- `patchids2 :: Vector{Int}`
- `patchids3 :: Vector{Int}`
- `patchids4 :: Vector{Int}`
- `inletpatchids :: Vector{Int}`
- `mu_resin_val :: Float`
- `N :: Int`

Meaning of the arguments:
- `i_interactive`: 1..if pressure inlet cells are selected manually, 0..else
- `celltype[i]`: Describes cell type. -1..pressure inlet, -2..pressure outlet, -3..cell with wall boundary, 1..interior cell
- `patchparameters0`,`patchparameters1`,`patchparameters2`,`patchparameters3`,`patchparameters4` :: 7-element vector with parameters for the main preform and the four sets if used as patch. The seven elements are cellporosity,cellthickness,cellpermeability,cellalpha,celldirection[1],celldirection[2],celldirection[3].
- `patchtype1val`,`patchtype2val`,`patchtype3val`,`patchtype4val` :: Type of patch. 1..
- `patchids1`,`patchids2`,`patchids3`,`patchids4` :: Type of set. 0..ignored, 1..pressure inlet, 2..patch, 3..pressure outlet
- `inletpatchids` :: Vector with the IDs of the cells which are inlet cells.
- `mu_resin_val` :: Kinematic viscosity value.
- `N`: Number of cells.

Return: Data structure with
- `cellthickness :: Vector{Float64}`
- `cellporosity :: Vector{Float64}`
- `cellpermeability :: Vector{Float64}`
- `cellalpha :: Vector{Float64}`
- `celldirection :: Array{Float64,2}`
- `cellviscosity :: Vector{Float64}`
- `celltype :: Vector{Int64}`
"""
function assign_parameters(input_struct)
    i_interactive=input_struct.i_interactive
    celltype=input_struct.celltype
    patchparameters0=input_struct.patchparameters0
    patchparameters1=input_struct.patchparameters1
    patchparameters2=input_struct.patchparameters2
    patchparameters3=input_struct.patchparameters3
    patchparameters4=input_struct.patchparameters4
    patchtype1val=input_struct.patchtype1val
    patchtype2val=input_struct.patchtype2val
    patchtype3val=input_struct.patchtype3val
    patchtype4val=input_struct.patchtype4val
    patchids1=input_struct.patchids1
    patchids2=input_struct.patchids2
    patchids3=input_struct.patchids3
    patchids4=input_struct.patchids4
    inletpatchids=input_struct.inletpatchids
    mu_resin_val=input_struct.mu_resin_val
    N=input_struct.N

    cellthickness=Vector{Float64}(undef, N)
    cellporosity=Vector{Float64}(undef, N)
    cellpermeability=Vector{Float64}(undef, N)
    cellalpha=Vector{Float64}(undef, N)
    celldirection=Array{Float64}(undef, N,3)
    cellviscosity=Vector{Float64}(undef, N)

    if i_interactive==0 || i_interactive==2
        if patchtype1val==1
            for i in 1:N
                ind=findfirst(isequal(i),patchids1)
                if ~isnothing(ind)
                    celltype[i]=-1
                end
            end                
        elseif patchtype1val==3
            for i in 1:N
                ind=findfirst(isequal(i),patchids1)
                if ~isnothing(ind)
                    celltype[i]=-2
                end
            end                  
        end
        if patchtype2val==1
            for i in 1:N
                ind=findfirst(isequal(i),patchids2)
                if ~isnothing(ind)
                    celltype[i]=-1
                end
            end
        elseif patchtype2val==3
            for i in 1:N
                ind=findfirst(isequal(i),patchids2)
                if ~isnothing(ind)
                    celltype[i]=-2
                end
            end
        end
        if patchtype3val==1
            for i in 1:N
                ind=findfirst(isequal(i),patchids3)
                if ~isnothing(ind)
                    celltype[i]=-1
                end
            end
        elseif patchtype3val==3
            for i in 1:N
                ind=findfirst(isequal(i),patchids3)
                if ~isnothing(ind)
                    celltype[i]=-2
                end
            end
        end
        if patchtype4val==1
            for i in 1:N
                ind=findfirst(isequal(i),patchids4)
                if ~isnothing(ind)
                    celltype[i]=-1
                end
            end
        elseif patchtype4val==3
            for i in 1:N
                ind=findfirst(isequal(i),patchids4)
                if ~isnothing(ind)
                    celltype[i]=-2
                end
            end
        end
    end
    if i_interactive==1 || i_interactive==2
        for i in 1:N
            ind=findfirst(isequal(i),inletpatchids)
            if ~isnothing(ind)
                celltype[i]=-1
            end
        end
    end        
    for ind in 1:N
        if i_interactive==1
            #thickness
            cellthickness[ind]=patchparameters0[2]
            #porosity
            cellporosity[ind]=patchparameters0[1] 
            #isotropic permeability 
            cellpermeability[ind]=patchparameters0[3]            
            #alpha permeability 
            cellalpha[ind]=patchparameters0[4]
            #primary direction
            vec=[patchparameters0[5] patchparameters0[6] patchparameters0[7]]
            celldirection[ind,:]=vec/sqrt(dot(vec,vec))
            #viscosity
            cellviscosity[ind]=mu_resin_val 
        else
            ind1=findfirst(isequal(ind),patchids1)
            ind2=findfirst(isequal(ind),patchids2)
            ind3=findfirst(isequal(ind),patchids3)
            ind4=findfirst(isequal(ind),patchids4)
            if (patchtype1val==2 && ~isnothing(ind1)) 
                patchparameters=patchparameters1
            elseif (patchtype2val==2 && ~isnothing(ind2)) 
                patchparameters=patchparameters2
            elseif (patchtype3val==2 && ~isnothing(ind3)) 
                patchparameters=patchparameters3
            elseif (patchtype4val==2 && ~isnothing(ind4)) 
                patchparameters=patchparameters4
            end
            if (patchtype1val==2 && issubset(ind,patchids1)) || (patchtype2val==2 && issubset(ind,patchids2)) || (patchtype3val==2 && issubset(ind,patchids3)) || (patchtype4val==2 && issubset(ind,patchids4))
                #thickness
                cellthickness[ind]=patchparameters[2]
                #porosity
                cellporosity[ind]=patchparameters[1] 
                #isotropic permeability 
                cellpermeability[ind]=patchparameters[3]            
                #alpha permeability 
                cellalpha[ind]=patchparameters[4]
                #primary direction
                vec=[patchparameters[5] patchparameters[6] patchparameters[7]]
                celldirection[ind,:]=vec/sqrt(dot(vec,vec))
                #viscosity
                cellviscosity[ind]=mu_resin_val 
            else
                #thickness
                cellthickness[ind]=patchparameters0[2]
                #porosity
                cellporosity[ind]=patchparameters0[1] 
                #isotropic permeability 
                cellpermeability[ind]=patchparameters0[3]            
                #alpha permeability 
                cellalpha[ind]=patchparameters0[4]
                #primary direction
                vec=[patchparameters0[5] patchparameters0[6] patchparameters0[7]]
                celldirection[ind,:]=vec/sqrt(dot(vec,vec))
                #viscosity
                cellviscosity[ind]=mu_resin_val 
            end
        end
    end

    return_struct=rtmsim.return_args_assign_parameters(cellthickness, cellporosity, cellpermeability, cellalpha, celldirection, cellviscosity, celltype)
    return return_struct
end

"""
    function create_coordinate_systems(input_struct)

Define the local cell coordinate system and the transformation matrix from the local cell coordinate system from the neighbouring cell to the local cell coordinate system of the considered cell.

Arguments: Data structure with
- `N :: Int`
- `cellgridid :: Array{Int,2}`
- `gridx :: Vector{Float}`
- `gridy :: Vector{Float}`
- `gridz :: Vector{Float}`
- `cellcenterx :: Vector{Float}`
- `cellcentery :: Vector{Float}`
- `cellcenterz :: Vector{Float}`
- `faces :: Array{Int,2}`
- `cellneighboursarray :: Array{Int,2}`
- `celldirection :: Array{Float,2}`
- `cellthickness :: Vector{Float}`
- `maxnumberofneighbours :: Int`

Meaning of the arguments:
- `N`: Number of cells
- `cellgridid`: The i-th line contains the three IDs of the nodes which form the cell
- `gridx`,`gridy`,`gridz`: i-th component of these vectors contain the x, y and z coordinates of node with ID i
- `cellcenterx`,`cellcentery`,`cellcenterz`: i-th component of these vectors contain the x, y and z coordinates of geometric cell centers of cell with ID i
- `faces`: Array with three columns. The three entries n1, n2, c1 in a cell are indices and describe a line between nodes with IDs n1 and n2 and this boundary belongs to cell with ID c1. 
- `cellneighboursarray`: The i-th line contains the indices of the neighbouring cells of cell with ID i. Number of columns is given by maxnumberofneighbours. Array is initialized with -9 and only the positive entries are considered.
- `celldirection`: The i-th line contains the x, y and z coordinates of the unit normal vector which is projected on the cell to define the cell coordinate system for cell with ID i.
- `cellthickness`: The i-th line contains the thickness of cell with ID i.
- `maxnumberofneighbours`: Number of columns of array `cellneighboursarray`. Default value is `10`. If more cell neighbours in the mesh, an error occurs and this value must be increased. 

Return: Data structure with
- `cellvolume :: Vector{Float}`
- `cellcentertocellcenterx :: Array{Float,2}`
- `cellcentertocellcentery :: Array{Float,2}` 
- `T11 :: Array{Float,2}`
- `T12 :: Array{Float,2}`
- `T21 :: Array{Float,2}`
- `T22 :: Array{Float,2}`
- `cellfacenormalx :: Array{Float,2}`
- `cellfacenormaly :: Array{Float,2}`
- `cellfacearea :: Array{Float,2}`
"""
function create_coordinate_systems(input_struct)
    N=input_struct.N
    cellgridid=input_struct.cellgridid
    gridx=input_struct.gridx
    gridy=input_struct.gridy
    gridz=input_struct.gridz
    cellcenterx=input_struct.cellcenterx
    cellcentery=input_struct.cellcentery
    cellcenterz=input_struct.cellcenterz
    faces=input_struct.faces
    cellneighboursarray=input_struct.cellneighboursarray
    celldirection=input_struct.celldirection
    cellthickness=input_struct.cellthickness
    maxnumberofneighbours=input_struct.maxnumberofneighbours

    cellvolume=Vector{Float64}(undef, N)
    cellcentertocellcenterx=Array{Float64}(undef, N, maxnumberofneighbours)
    cellcentertocellcentery=Array{Float64}(undef, N, maxnumberofneighbours)
    T11=Array{Float64}(undef, N, maxnumberofneighbours)
    T12=Array{Float64}(undef, N, maxnumberofneighbours)
    T21=Array{Float64}(undef, N, maxnumberofneighbours)
    T22=Array{Float64}(undef, N, maxnumberofneighbours)
    cellfacenormalx=Array{Float64}(undef, N, maxnumberofneighbours)
    cellfacenormaly=Array{Float64}(undef, N, maxnumberofneighbours)
    cellfacearea=Array{Float64}(undef, N, maxnumberofneighbours)
    for ind in 1:N
        cellvolume[ind]=-9
        for ind_n in 1:maxnumberofneighbours
            cellcentertocellcenterx[ind,ind_n]=-9.0
            cellcentertocellcentery[ind,ind_n]=-9.0
            T11[ind,ind_n]=-9.0
            T12[ind,ind_n]=-9.0
            T21[ind,ind_n]=-9.0
            T22[ind,ind_n]=-9.0
            cellfacenormalx[ind,ind_n]=-9.0
            cellfacenormaly[ind,ind_n]=-9.0
            cellfacearea[ind,ind_n]=-9.0
        end
    end
    b1=Array{Float64}(undef, N, 3)
    b2=Array{Float64}(undef, N, 3)
    b3=Array{Float64}(undef, N, 3)
    gridxlocal=Array{Float64}(undef, N, 3);00
    gridylocal=Array{Float64}(undef, N, 3)
    gridzlocal=Array{Float64}(undef, N, 3)
    theta=Vector{Float64}(undef, N)

    for ind in 1:N
        # First, an intermediate orthonormal basis {b1, b2, b3} is created with 
        # first direction pointing from node with smallest ID to node with medium ID, 
        # second direction pointing in the orthogonal component from node with smallest ID to node with highest ID and 
        # third direction given by the cross product of the first two directions. 
        # The origin of the local coordinate system is the geometric center of the triangular cell. 
        i1=cellgridid[ind,1]
        i2=cellgridid[ind,2]
        i3=cellgridid[ind,3]  
        b1[ind,1:3]=[gridx[i2]-gridx[i1] gridy[i2]-gridy[i1] gridz[i2]-gridz[i1]]
        b1[ind,1:3]=b1[ind,1:3]/sqrt(dot(b1[ind,1:3],b1[ind,1:3]))
        a2=[gridx[i3]-gridx[i1] gridy[i3]-gridy[i1] gridz[i3]-gridz[i1]]'
        a2=a2/sqrt(dot(a2,a2))
        b2[ind,1:3]=a2-dot(b1[ind,1:3],a2)/dot(b1[ind,1:3],b1[ind,1:3])*b1[ind,1:3]
        b2[ind,1:3]=b2[ind,1:3]/sqrt(dot(b2[ind,1:3],b2[ind,1:3]))
        b3[ind,1:3]=cross(b1[ind,1:3],b2[ind,1:3])   

        # Then the reference vector is formulated in the intermediate orthonormal basis 
        Tmat=[b1[ind,1] b2[ind,1] b3[ind,1]; b1[ind,2] b2[ind,2] b3[ind,2]; b1[ind,3] b2[ind,3] b3[ind,3]]
        xvec=celldirection[ind,:]
        bvec=Tmat\xvec
        r1=[bvec[1] bvec[2] bvec[3]]'  #ref dir in local CS

        # In order to get the local coordinate system the basis {b1, b2, b3} is rotated by angle theta about the b3-axis.
        # Calculate the angle by which b1 must be rotated about the b3-axis to match r1 via relation rotation matrix Rz(theta)*[1;0;0]=r1, i.e. cos(theta)=r1(1) and sin(theta)=r1(2)
        theta[ind]=atan(r1[2],r1[1])
        #Rotation of theta about nvec=b3 to get c1 and c2 
        nvec=b3[ind,:]
        xvec=b1[ind,:]
        c1=nvec*dot(nvec,xvec)+cos(theta[ind])*cross(cross(nvec,xvec),nvec)+sin(theta[ind])*cross(nvec,xvec)
        xvec=b2[ind,:]
        c2=nvec*dot(nvec,xvec)+cos(theta[ind])*cross(cross(nvec,xvec),nvec)+sin(theta[ind])*cross(nvec,xvec)
        xvec=b3[ind,:]
        c3=nvec*dot(nvec,xvec)+cos(theta[ind])*cross(cross(nvec,xvec),nvec)+sin(theta[ind])*cross(nvec,xvec)
        b1[ind,:]=c1
        b2[ind,:]=c2
        b3[ind,:]=c3  
    
        #transformation of vertices into local CS
        gridxlocal[ind,1]=gridx[i1]-cellcenterx[ind]
        gridylocal[ind,1]=gridy[i1]-cellcentery[ind]
        gridzlocal[ind,1]=gridz[i1]-cellcenterz[ind]
        gridxlocal[ind,2]=gridx[i2]-cellcenterx[ind]
        gridylocal[ind,2]=gridy[i2]-cellcentery[ind]
        gridzlocal[ind,2]=gridz[i2]-cellcenterz[ind]
        gridxlocal[ind,3]=gridx[i3]-cellcenterx[ind]
        gridylocal[ind,3]=gridy[i3]-cellcentery[ind]
        gridzlocal[ind,3]=gridz[i3]-cellcenterz[ind]
        Tmat=[b1[ind,1] b2[ind,1] b3[ind,1]; b1[ind,2] b2[ind,2] b3[ind,2]; b1[ind,3] b2[ind,3] b3[ind,3]]
        xvec=[gridxlocal[ind,1] gridylocal[ind,1] gridzlocal[ind,1]]' 
        bvec=Tmat\xvec
        gridxlocal[ind,1]=bvec[1];gridylocal[ind,1]=bvec[2];gridzlocal[ind,1]=bvec[3]
        xvec=[gridxlocal[ind,2] gridylocal[ind,2] gridzlocal[ind,2]]' 
        bvec=Tmat\xvec
        gridxlocal[ind,2]=bvec[1];gridylocal[ind,2]=bvec[2];gridzlocal[ind,2]=bvec[3]
        xvec=[gridxlocal[ind,3] gridylocal[ind,3] gridzlocal[ind,3]]' 
        bvec=Tmat\xvec
        gridxlocal[ind,3]=bvec[1];gridylocal[ind,3]=bvec[2];gridzlocal[ind,3]=bvec[3]
    end

    cellids=[Int64(-9) Int64(-9)]
    gridids=[Int64(-9) Int64(-9)]
    x=[-9.0, -9.0, -9.0]
    x0=[-9.0, -9.0, -9.0]
    r0=[-9.0, -9.0, -9.0]
    gridxlocal_neighbour=[-9.0, -9.0, -9.0]
    gridylocal_neighbour=[-9.0, -9.0, -9.0]
    gridzlocal_neighbour=[-9.0, -9.0, -9.0]
    f1=[-9.0, -9.0, -9.0]
    f2=[-9.0, -9.0, -9.0]
    f3=[-9.0, -9.0, -9.0]
    # In a next step the flattened geometry is created, i.e. the cell center and
    # the non-common node of the neighbouring cell is rotated about
    # the common edge to lie in the plane of the considered cell with ID ind
    for ind in 1:N
        cellneighboursline=cellneighboursarray[ind,:]
        cellneighboursline=cellneighboursline[cellneighboursline .> 0]
        for i_neighbour in 1:length(cellneighboursline)
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
            #face in between face direction is from smaller to larger index
            #x0..local coordinates of smaller index
            #r0..vector from smaller to larger index in LCS
            inds1=findall(isequal(ind),faces[:,3])
            inds2=findall(isequal(cellneighboursline[i_neighbour]),faces[:,3])
            mat1=faces[inds1,:]
            mat2=faces[inds2,:]
            mat3=vcat(mat1,mat2) 
            mat4=sortslices(mat3,dims=1)
            for irow in 1:size(mat4,1)-1
                if mat4[irow,1]==mat4[irow+1,1] && mat4[irow,2]==mat4[irow+1,2]
                    if mat4[irow,3]==ind
                        cellids=[ind mat4[irow+1,3]]
                    else
                        cellids=[ind mat4[irow,3]]
                    end
                    gridids=[mat4[irow,1] mat4[irow,2]]
                end
            end
            inds=[cellgridid[ind,1], cellgridid[ind,2], cellgridid[ind,3]]
            ia=findall(isequal(gridids[1]),inds)
            ib=findall(isequal(gridids[2]),inds)
            x0=[gridxlocal[ind,ia], gridylocal[ind,ia], gridzlocal[ind,ia]]
            r0=[gridxlocal[ind,ib]-gridxlocal[ind,ia], gridylocal[ind,ib]-gridylocal[ind,ia], gridzlocal[ind,ib]-gridzlocal[ind,ia]]

            #Define xvec as the vector between cell centers ind and neighbouring cell center (A) (in GCS) 
            #and transform xvec into local coordinates bvec, this gives A in LCS.
            #Find normal distance from A in LCS to the cell boundary with that cell center A in flat geometry and 
            #face normal vector can be defined.
            x=[[0.0], [0.0], [0.0]]  #P at origin of local CS
            Px=x[1]
            Py=x[2]
            Pz=x[3]
            lambda=dot(x-x0,r0)/dot(r0,r0)  
            Q1x=x0[1]+lambda*r0[1]
            Q1y=x0[2]+lambda*r0[2]
            Q1z=x0[3]+lambda*r0[3]
            vec1=[Px-Q1x, Py-Q1y, Pz-Q1z]
            l1=sqrt(dot(vec1,vec1)) 
            nvec=[(Q1x-Px), (Q1y-Py), (Q1z-Pz)]
            nvec=nvec/sqrt(dot(nvec,nvec))
            cellfacenormalx[ind,i_neighbour]=only(nvec[1])
            cellfacenormaly[ind,i_neighbour]=only(nvec[2]) 

            Tmat=[b1[ind,1] b2[ind,1] b3[ind,1]; b1[ind,2] b2[ind,2] b3[ind,2]; b1[ind,3] b2[ind,3] b3[ind,3]]
            xvec=[cellcenterx[cellneighboursarray[ind,i_neighbour]]-cellcenterx[ind], cellcentery[cellneighboursarray[ind,i_neighbour]]-cellcentery[ind], cellcenterz[cellneighboursarray[ind,i_neighbour]]-cellcenterz[ind] ]  #A in global CS
            bvec=Tmat\xvec
            x=[[bvec[1]], [bvec[2]], [bvec[3]]] #A in local CS
            Ax=x[1]
            Ay=x[2]
            Az=x[3]
            lambda=dot(x-x0,r0)/dot(r0,r0)
            Q2x=x0[1]+lambda*r0[1]
            Q2y=x0[2]+lambda*r0[2]
            Q2z=x0[3]+lambda*r0[3]
            vec2=[Ax-Q2x, Ay-Q2y, Az-Q2z]
            l2=sqrt(dot(vec2,vec2))
            cellcentertocellcenterx[ind,i_neighbour]=only(Px+(Q1x-Px)+(Q2x-Q1x)+l2/l1*(Q1x-Px))
            cellcentertocellcentery[ind,i_neighbour]=only(Py+(Q1y-Py)+(Q2y-Q1y)+l2/l1*(Q1y-Py))

            vec3=[gridxlocal[ind,ib]-gridxlocal[ind,ia], gridylocal[ind,ib]-gridylocal[ind,ia], gridzlocal[ind,ib]-gridzlocal[ind,ia]]
            cellfacearea[ind,i_neighbour]=0.5*(cellthickness[cellids[1]]+cellthickness[cellids[2]])*sqrt(dot(vec3,vec3))

            #Transformation matrix for (u,v) of neighbouring cells to local coordinate system.
            #Find the two common grid points and the third non-common grid point               
            ind21=-9  #Issues with setdiff, therefore manual implementation  #setdiff(cellgridid[cellids[2],:],gridids)
            for ind_tmp in 1:3
                if cellgridid[cellids[2],ind_tmp]!=gridids[1] && cellgridid[cellids[2],ind_tmp]!=gridids[2]
                    ind21=cellgridid[cellids[2],ind_tmp]
                end
            end     
            thirdgrid=only(findall(isequal(ind21),cellgridid[cellids[2],:]))
            common1grid=only(findall(isequal(gridids[1]),cellgridid[cellids[2],:]))
            common2grid=only(findall(isequal(gridids[2]),cellgridid[cellids[2],:]))   
            #construction of the third one in outside normal direction for the flat geometry
            #based on the length of the two non-common edges
            gridxlocal_neighbour[2]=only(gridxlocal[ind,ia])  #gridxlocal(ind,common1grid)
            gridxlocal_neighbour[3]=only(gridxlocal[ind,ib])  #gridxlocal(ind,common2grid)
            gridylocal_neighbour[2]=only(gridylocal[ind,ia])  #gridylocal(ind,common1grid)
            gridylocal_neighbour[3]=only(gridylocal[ind,ib])  #gridylocal(ind,common2grid)
            gridzlocal_neighbour[2]=0.0
            gridzlocal_neighbour[3]=0.0
            
            ind3=-9
            for ind_tmp in 1:3
                if cellgridid[cellids[2],ind_tmp]!=cellgridid[ind,1] && cellgridid[cellids[2],ind_tmp]!=cellgridid[ind,2] && cellgridid[cellids[2],ind_tmp]!=cellgridid[ind,3]
                    ind3=cellgridid[cellids[2],ind_tmp]
                end
            end
            Tmat=[b1[ind,1] b2[ind,1] b3[ind,1]; b1[ind,2] b2[ind,2] b3[ind,2]; b1[ind,3] b2[ind,3] b3[ind,3]]
            xvec=[gridx[ind3]-cellcenterx[ind], gridy[ind3]-cellcentery[ind], gridz[ind3]-cellcenterz[ind]] #A in global CS
            bvec=Tmat\xvec
            x=[[bvec[1]], [bvec[2]], [bvec[3]]] #A in local CS
            Ax=x[1]
            Ay=x[2]
            Az=x[3]
            lambda=dot(x-x0,r0)/dot(r0,r0)
            Q2x=x0[1]+lambda*r0[1]
            Q2y=x0[2]+lambda*r0[2]
            Q2z=x0[3]+lambda*r0[3]
            vec2=[Ax-Q2x, Ay-Q2y, Az-Q2z]
            l2=sqrt(dot(vec2,vec2))
            gridxlocal_neighbour[1]=only(Px+(Q1x-Px)+(Q2x-Q1x)+l2/l1*(Q1x-Px))
            gridylocal_neighbour[1]=only(Py+(Q1y-Py)+(Q2y-Q1y)+l2/l1*(Q1y-Py))
            gridzlocal_neighbour[1]=0.

            #Construction of LCS f1,f2,f3 according to procedure from above using the points gridxlocal_neighbour(j),gridylocal_neighbour(j)
            ivec1=[only(cellgridid[cellids[2],1]), only(cellgridid[cellids[2],2]), only(cellgridid[cellids[2],3])]                           
            min_val=min(ivec1[1],ivec1[2],ivec1[3]) 
            max_val=max(ivec1[1],ivec1[2],ivec1[3]) 
            idel1=findall(isequal(min(ivec1[1],ivec1[2],ivec1[3])),ivec1);deleteat!(ivec1,idel1);idel1=findall(isequal(max(ivec1[1],ivec1[2])),ivec1);deleteat!(ivec1,idel1)
            median_val=ivec1[1]
            if ind3==min_val; k1=1; elseif ind3==median_val; k2=1; elseif ind3==max_val; k3=1; end                
            ind4=cellgridid[cellids[2],common1grid]
            if ind4==min_val; k1=2; elseif ind4==median_val; k2=2; elseif ind4==max_val; k3=2; end             
            ind5=cellgridid[cellids[2],common2grid]
            if ind5==min_val; k1=3; elseif ind5==median_val; k2=3; elseif ind5==max_val; k3=3; end
    
            f1=[gridxlocal_neighbour[k2]-gridxlocal_neighbour[k1], gridylocal_neighbour[k2]-gridylocal_neighbour[k1], gridzlocal_neighbour[k2]-gridzlocal_neighbour[k1]]
            f1=f1/sqrt(dot(f1,f1))
            a2=[gridxlocal_neighbour[k3]-gridxlocal_neighbour[k1], gridylocal_neighbour[k3]-gridylocal_neighbour[k1], gridzlocal_neighbour[k3]-gridzlocal_neighbour[k1]]
            a2=a2/sqrt(dot(a2,a2))
            f2=a2-dot(f1,a2)/dot(f1,f1)*f1
            f2=f2/sqrt(dot(f2,f2))
            f3=cross(f1,f2)    
    
            nvec=f3
            xvec=f1
            c1=nvec*dot(nvec,xvec)+cos(theta[cellneighboursarray[ind,i_neighbour]])*cross(cross(nvec,xvec),nvec)+sin(theta[cellneighboursarray[ind,i_neighbour]] )*cross(nvec,xvec)
            xvec=f2
            c2=nvec*dot(nvec,xvec)+cos(theta[cellneighboursarray[ind,i_neighbour]])*cross(cross(nvec,xvec),nvec)+sin(theta[cellneighboursarray[ind,i_neighbour]] )*cross(nvec,xvec)
            xvec=f3
            c3=nvec*dot(nvec,xvec)+cos(theta[cellneighboursarray[ind,i_neighbour]])*cross(cross(nvec,xvec),nvec)+sin(theta[cellneighboursarray[ind,i_neighbour]] )*cross(nvec,xvec)
            f1=c1
            f2=c2
            f3=c3
            Tmat=[f1[1] f2[1] f3[1]; f1[2] f2[2] f3[2]; f1[3] f2[3] f3[3]]

            #Assign transformation matrix for the velocities in the local coordinate systems
            #(u,v)_e=T*(u,v)_f
            T11[ind,i_neighbour]=Tmat[1,1]
            T12[ind,i_neighbour]=Tmat[1,2]
            T21[ind,i_neighbour]=Tmat[2,1]
            T22[ind,i_neighbour]=Tmat[2,2]
        end

        #calculate cell volume
        vec1=[gridxlocal[ind,2]-gridxlocal[ind,1], gridylocal[ind,2]-gridylocal[ind,1], gridzlocal[ind,2]-gridzlocal[ind,1]]
        vec2=[gridxlocal[ind,3]-gridxlocal[ind,1], gridylocal[ind,3]-gridylocal[ind,1], gridzlocal[ind,3]-gridzlocal[ind,1]]
        vec3=cross(vec1,vec2)
        cellvolume[ind]=cellthickness[ind]*0.5*sqrt(dot(vec3,vec3))
    end  

    return_struct=rtmsim.return_args_create_cs(cellvolume, cellcentertocellcenterx, cellcentertocellcentery, T11, T12, T21, T22, cellfacenormalx, cellfacenormaly, cellfacearea)
    return return_struct
end


