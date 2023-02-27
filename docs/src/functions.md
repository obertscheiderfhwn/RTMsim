# Functions
Alternatively to using the GUI, one has access to all the functions after installing the package.


## API Functions
```@docs
rtmsim.plot_mesh(meshfilename,i_mode)
```

```@docs
rtmsim.plot_sets(meshfilename)
```

```@docs
rtmsim.read_nastran_mesh(meshfilename,paramset,paramset1,paramset2,paramset3,paramset4,patchtype1val,patchtype2val,patchtype3val,patchtype4val,i_interactive,r_p)
```

```@docs
rtmsim.start_rtmsim(inputfilename)
```

```@docs
rtmsim.rtmsim_rev1(param)
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
rtmsim.numerical_gradient(i_method,ind,p_old,cellneighboursarray,cellcentertocellcenterx,cellcentertocellcentery)
```

```@docs
rtmsim.numerical_flux_function(i_method,vars_P,vars_A,meshparameters)
```

```@docs
rtmsim.numerical_flux_function_boundary(i_method,vars_P,vars_A,meshparameters,n_dot_u)
```



## Additional features

The source code is prepared for the following extensions:
- Import mesh file in different format. Selection is based on the extension of the mesh file.
- Input parameter `i_model` (for iso-thermal RTM `=1`) is used for adding additional functionalities. E.g. adding temperature and degree-of-cure equations with variable resin viscosity or for VARI with variable porosity, permeability and cavity thickness.
- Parameter `i_method` in the functions for numerical differentiation and flux functions can be used to implement different numerical schemes. E.g. gradient limiter or second-order upwinding.

If you are interested, please have a look at the contribution item in the [community standards](https://github.com/obertscheiderfhwn/RTMsim/community).





```@meta
EditURL = "https://github.com/obertscheiderfhwn/RTMsim/blob/main/docs/src/functions.md"
```