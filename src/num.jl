"""
    function numerical_gradient(input_struct)
        
Calculates the pressure gradient from the cell values of the neighbouring cells.
- `i_method=1` .. Least square solution to determine gradient
- `i_method=2` .. Least square solution to determine gradient with limiter
- `i_method=3` .. RUntime optimized least square solution to determine gradient

Arguments: Data structure with
- `i_method :: Int`
- `ind :: Int`
- `p_old :: Vector{Float}`
- `cellneighoursarray :: Array{Float,2}`
- `cellcentertocellcenterx, cellcentertocellcentery :: Array{Float,2}`

Return: Data structure with
- `dpdx :: Float`
- `dpdy :: Float`
"""
function numerical_gradient(input_struct)
    i_method=input_struct.i_method
    ind=input_struct.ind
    p_old=input_struct.p_old
    cellneighboursarray=input_struct.cellneighboursarray
    cellcentertocellcenterx=input_struct.cellcentertocellcenterx
    cellcentertocellcentery=input_struct.cellcentertocellcentery

    if i_method==1
        #least square solution to determine gradient
        cellneighboursline=cellneighboursarray[ind,:]
        cellneighboursline=cellneighboursline[cellneighboursline .> 0]
        len_cellneighboursline=length(cellneighboursline)
        bvec=Vector{Float64}(undef,len_cellneighboursline)
        Amat=Array{Float64}(undef,len_cellneighboursline,2)  
        for i_neighbour in 1:len_cellneighboursline
            i_P=ind
            i_A=cellneighboursarray[ind,i_neighbour]  
            Amat[i_neighbour,1]=cellcentertocellcenterx[ind,i_neighbour]
            Amat[i_neighbour,2]=cellcentertocellcentery[ind,i_neighbour]
            bvec[i_neighbour]=p_old[i_A]-p_old[i_P]
        end

        if len_cellneighboursline>1
            xvec=Amat[1:len_cellneighboursline,:]\bvec[1:len_cellneighboursline]
            dpdx=xvec[1]
            dpdy=xvec[2]        
        else
            dpdx=0
            dpdy=0
        end
    elseif i_method==2
        #least square solution to determine gradient with limiter
        cellneighboursline=cellneighboursarray[ind,:]
        cellneighboursline=cellneighboursline[cellneighboursline .> 0]
        len_cellneighboursline=length(cellneighboursline)
        bvec=Vector{Float64}(undef,len_cellneighboursline)
        Amat=Array{Float64}(undef,len_cellneighboursline,2)  
        wi=Vector{Float64}(undef,len_cellneighboursline)
        for i_neighbour in 1:len_cellneighboursline
            i_P=ind
            i_A=cellneighboursarray[ind,i_neighbour]  
            exp_limiter=2
            wi[i_neighbour]=1/(sqrt((cellcentertocellcenterx[ind,i_neighbour])^2+(cellcentertocellcentery[ind,i_neighbour])^2))^exp_limiter
            Amat[i_neighbour,1]=wi[i_neighbour]*cellcentertocellcenterx[ind,i_neighbour]
            Amat[i_neighbour,2]=wi[i_neighbour]*cellcentertocellcentery[ind,i_neighbour]
            bvec[i_neighbour]=wi[i_neighbour]*(p_old[i_A]-p_old[i_P])
        end

        if len_cellneighboursline>1
            xvec=Amat[1:len_cellneighboursline,:]\bvec[1:len_cellneighboursline]
            dpdx=xvec[1]
            dpdy=xvec[2]            
        else
            dpdx=0
            dpdy=0
        end
    elseif i_method==3
        #least square solution to determine gradient - runtime optimized
        cellneighboursline=cellneighboursarray[ind,:]
        cellneighboursline=cellneighboursline[cellneighboursline .> 0]
        len_cellneighboursline=length(cellneighboursline)
        bvec=Vector{Float64}(undef,len_cellneighboursline)
        Amat=Array{Float64}(undef,len_cellneighboursline,2)  
        for i_neighbour in 1:len_cellneighboursline
            i_P=ind
            i_A=cellneighboursarray[ind,i_neighbour]  
            Amat[i_neighbour,1]=cellcentertocellcenterx[ind,i_neighbour]
            Amat[i_neighbour,2]=cellcentertocellcentery[ind,i_neighbour]
            bvec[i_neighbour]=p_old[i_A]-p_old[i_P]
        end
        #xvec=Amat[1:len_cellneighboursline,:]\bvec[1:len_cellneighboursline]
        #dpdx=xvec[1]
        #dpdy=xvec[2]

        if len_cellneighboursline>1
            Aplus=transpose(Amat)*Amat
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
            dpdx=0
            dpdy=0
        end

    end
    return_struct=rtmsim.return_args_gradient(dpdx,dpdy)
    return return_struct
end


"""
    function numerical_flux_function(input_struct)

Evaluates the numerical flux functions at the cell boundaries.
-` i_method==1` .. first order upwinding

Arguments: Data structure with
- `i_method :: Int`
- `vars_P, vars_A :: 4-element Vector{Float}`
- `meshparameters :: 3-element Vector{Float}`

Return: Data structure with
- `F_rho_num_add :: Float`
- `F_u_num_add :: Float`
- `F_v_num_add :: Float`
- `F_gamma_num_add :: Float`
- `F_gamma_num1_add :: Float`

Unit tests:
- `input_struct=rtmsim.input_args_flux(1,[1.0; 1.2; 0.0; 0.0],[1.0; 0.4; 0.0; 0.0],[1.0;0.0;1.0]);return_struct=rtmsim.numerical_flux_function(input_struct);println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add))`

    gives a comparison with the one-dimensional numerical flux for the continuity and momentum equations. The result for `i_mehtod=1` is `F_rho=0.8, F_u=0.96`. The first is calculated without upwinding `F_rho=rho*u=0.5*(1.0+1.0)*0.5*(1.2+0.4)=0.8`, the second is calculated without upwinding in the first factor and with upwinding in the second factor `F_u=(rho*u)*(u)=0.8*1.2=0.96`. Mass density and velocity are `[rho,u]=[1.0,1.2]` in the considered and `[rho,u]=[1.0,0.4]` in the neighbouring cells. `y` velocity and fluid fraction `gamma` are not considered in this test case. 

- `input_struct=rtmsim.input_args_flux(1,[1.225; 1.2; 0.4; 0.9],[1.0; 0.4; 1.2; 0.1],[1/sqrt(2);1/sqrt(2);1.0]);return_struct=rtmsim.numerical_flux_function(input_struct);println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma_num1="*string(return_struct.F_gamma_num1_add))`

    is a general test case with `[rho,u,v,gamma]=[1.225,1.2,0.4,0.9]` in the considered cell, `[rho,u,v,gamma]=[1.225,1.2,0.4,0.9]` in the neighbouring cell, outwards pointing cell normal vector `[1/sqrt(2),1/sqrt(2)]` and boundary face area `1.0`. The result for `i_mehtod=1` is `F_rho=1.2586500705120547, F_u=1.5103800846144655, F_v=0.5034600282048219, F_gamma=1.0182337649086284, F_gamma_num1=1.131370849898476`.

"""
function numerical_flux_function(input_struct)
    i_method=input_struct.i_method
    vars_P=input_struct.vars_P
    vars_A=input_struct.vars_A
    meshparameters=input_struct.meshparameters
    if i_method==1
        #first order upwinding
        rho_P=vars_P[1]
        u_P=vars_P[2]
        v_P=vars_P[3]
        gamma_P=vars_P[4]
        rho_A=vars_A[1]
        u_A=vars_A[2]
        v_A=vars_A[3]
        gamma_A=vars_A[4]
        n_x=meshparameters[1]
        n_y=meshparameters[2]
        A=meshparameters[3]
        n_dot_rhou=dot([n_x; n_y],0.5*(rho_P+rho_A)*[0.5*(u_P+u_A); 0.5*(v_P+v_A)])
        phi=1
        F_rho_num_add=n_dot_rhou*phi*A
        if n_dot_rhou>=0
            phi=u_P                                
        else
            phi=u_A
        end
        F_u_num_add=n_dot_rhou*phi*A     
        if n_dot_rhou>=0
            phi=v_P  
        else
            phi=v_A
        end
        F_v_num_add=n_dot_rhou*phi*A 
        n_dot_u=dot([n_x; n_y],[0.5*(u_P+u_A); 0.5*(v_P+v_A)])
        if n_dot_u>=0 
            phi=gamma_P  
        else
            phi=gamma_A
        end  
        F_gamma_num_add=n_dot_u*phi*A
        phi=1
        F_gamma_num1_add=n_dot_u*phi*A
    else
        errorstring=string("i_method=",string(i_method)," not implemented"* "\n") 
        error(errorstring)
    end
    return_struct=rtmsim.return_args_flux(F_rho_num_add,F_u_num_add,F_v_num_add,F_gamma_num_add,F_gamma_num1_add)
    return return_struct
end


"""
    numerical_flux_function_boundary(input_struct)

Evaluates the numerical flux functions at the cell boundaries to pressure inlet or outlet.
- `i_method==1` .. first order upwinding

Arguments: Data structure with
- `i_method :: Int`
- `vars_P, vars_A :: ` 4-element ` Vector{Float}`
- `meshparameters :: ` 3-element ` Vector{Float}`
- `n_dot_u :: Float`

Return: Data structure with
- `F_rho_num_add :: Float`
- `F_u_num_add :: Float`
- `F_v_num_add :: Float`
- `F_gamma_num_add :: Float`
- `F_gamma_num1_add :: Float` 

Unit tests:
- `input_struct=rtmsim.input_args_flux_boundary(1,[1.0; 1.2; 0.0; 0.0],[1.0; 0.0; 0.0; 0.0],[1.0;0.0;1.0],1.2);return_struct=rtmsim.numerical_flux_function(input_struct);println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add))`

    gives a comparison with the one-dimensional inflow from a cell with mass density `rho=1.0` with velocity `u=1.2` at the cell boundary. Mass density and velocity in the considered cell are `rho=1.0` and `u=1.2`. `y` velocity and fluid fraction `gamma` are not considered in this test case. The result for `i_mehtod=1` is `F_rho=0.6, F_u=0.72`.

- `input_struct=rtmsim.input_args_flux_boundary(1,[1.225; 1.2; 0.4; 0.9],[1.0; 0.4; 1.2; 0.1],[1/sqrt(2);1/sqrt(2);1.0],-0.8);return_struct=rtmsim.numerical_flux_function_boundary(input_struct);println("F_rho="*string(return_struct.F_rho_num_add)*", F_u="*string(return_struct.F_u_num_add)*", F_v="*string(return_struct.F_v_num_add)*", F_gamma="*string(return_struct.F_gamma_num_add)*", F_gamma_num1="*string(return_struct.F_gamma_num1_add))`

    is a general case. The result for `i_mehtod=1` is `F_rho=-0.8900000000000001, F_u=-0.3560000000000001, F_v=-1.068, F_gamma=-0.08000000000000002, F_gamma_num1=-0.8`.

"""
function numerical_flux_function_boundary(input_struct)
    i_method=input_struct.i_method
    vars_P=input_struct.vars_P
    vars_A=input_struct.vars_A
    meshparameters=input_struct.meshparameters
    n_dot_u=input_struct.n_dot_u
    if i_method==1
        #first order upwinding
        rho_P=vars_P[1]
        u_P=vars_P[2]
        v_P=vars_P[3]
        gamma_P=vars_P[4]
        rho_A=vars_A[1]
        u_A=vars_A[2]
        v_A=vars_A[3]
        gamma_A=vars_A[4]
        n_x=meshparameters[1]
        n_y=meshparameters[2]
        A=meshparameters[3]        
        n_dot_rhou=n_dot_u*0.5*(rho_A+rho_P)
        phi=1
        F_rho_num_add=n_dot_rhou*phi*A
        if n_dot_u<=0 
            phi=u_A                   
        else
            phi=u_P
        end
        F_u_num_add=n_dot_rhou*phi*A
        if n_dot_u<=0 
            phi=v_A                   
        else
            phi=v_P
        end
        F_v_num_add=n_dot_rhou*phi*A
        if n_dot_u<=0 
            phi=gamma_A                   
        else
            phi=gamma_P
        end
        F_gamma_num_add=n_dot_u*phi*A
        phi=1
        F_gamma_num1_add=n_dot_u*phi*A
    end
    return_struct=rtmsim.return_args_flux(F_rho_num_add,F_u_num_add,F_v_num_add,F_gamma_num_add,F_gamma_num1_add)
    return return_struct
end 

"""
    burgers(i_method)

Computation of the Burgers equation
"""
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