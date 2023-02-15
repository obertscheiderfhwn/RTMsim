# Functions
Alternatively to using the GUI, one has access to all the functions after installing the package.


## API Functions
```@docs
rtmsim.plot_mesh(meshfilename,i_mode)
```

```@docs
rtmsim.read_nastran_mesh(meshfilename,paramset,paramset1,paramset2,paramset3,paramset4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_interactive,r_p)
```

```@docs
rtmsim.start_rtmsim(inputfilename)
```

```@docs
rtmsim.rtmsim_rev1(i_model,meshfilename,tmax,
        p_ref,rho_ref,gamma,mu_resin_val,
        p_a_val,p_init_val,
        t_val,porosity_val,K_val,alpha_val,refdir1_val,refdir2_val,refdir3_val,
        t1_val,porosity1_val,K1_val,alpha1_val,refdir11_val,refdir21_val,refdir31_val,
        t2_val,porosity2_val,K2_val,alpha2_val,refdir12_val,refdir22_val,refdir32_val,
        t3_val,porosity3_val,K3_val,alpha3_val,refdir13_val,refdir23_val,refdir33_val,
        t4_val,porosity4_val,K4_val,alpha4_val,refdir14_val,refdir24_val,refdir34_val,
        patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_restart,restartfilename,i_interactive,r_p,n_pics)
```


```@docs
rtmsim.plot_results(resultsfilename)
```

```@docs
rtmsim.plot_overview(n_out,n_pics)
```

```@docs
rtmsim.plot_filling(n_out,n_pics)
```

```@docs
rtmsim.gui()
```

```@docs
rtmsim.function numerical_gradient(i_method,ind,p_old,cellneighboursarray,cellcentertocellcenterx,cellcentertocellcentery)
```

```@docs
rtmsim.function numerical_flux_function(i_method,vars_P,vars_A,meshparameters)
```

```@docs
rtmsim.function numerical_flux_function_boundary(i_method,vars_P,vars_A,meshparameters,n_dot_u)
```



## Additional features

The source code is prepared for the following extensions:
- Import mesh file in different format. Selection is based on the extension of the mesh file.
- Input parameter `i_model` (for iso-thermal RTM `=1`) is used for adding additional functionalities. E.g. adding temperature and degree-of-cure equations with variable resin viscosity or for VARI with variable porosity, permeability and cavity thickness.
- Parameter `i_method` in the functions for numerical differentiation and flux functions can be used to implement different numerical schemes. E.g. gradient limiter or second-order upwinding.

E.g. if wiggles (oscillations) are present in the pressure contour plot, a gradient limiter is used in the function `numerical_gradient`. The function is called with i_method=2 as first argument: 
`dpdx,dpdy=numerical_gradient(2,ind,p_old,cellneighboursarray,cellcentertocellcenterx,cellcentertocellcentery);` <br>
In the function a case selection determines the used method:
```
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
        xvec=Amat[1:len_cellneighboursline,:]\bvec[1:len_cellneighboursline];
        dpdx=xvec[1];
        dpdy=xvec[2];
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
        xvec=Amat[1:len_cellneighboursline,:]\bvec[1:len_cellneighboursline];
        dpdx=xvec[1];
        dpdy=xvec[2];
    end
    return dpdx,dpdy
end
```
After modifying and compiling the RTMsim module, a simulation can be started in the GUI or from the terminal.





```@meta
EditURL = "https://github.com/obertscheiderfhwn/RTMsim/blob/main/docs/src/functions.md"
```