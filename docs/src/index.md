# RTMsim
RTMsim is a new software tool for RTM filling simulations which fulfills these requirements: Several test cases were used for successfully validating the implemented model. The porous cavity is fully described by a mesh file with triangular cells on the part’s mid-surface and cell set definitions. The latter can be used for specifying the location of the pressure injection ports and regions with different preforms by assigning different thickness, porosity and permeability values. Additional equations (e.g. for modeling the degree-of-cure) can either be added with equations of the same type or modifications of existing equations (e.g. for variable cavity thickness as needed for vacuum assisted resin infusion simulations).

## Package Features

- The simulation model shall give correct results for filling pattern and filling time.
- The simulation tool takes only composite-manufacturing related inputs and the simulation shall be robust independent of numerics-related input.
- The simulation tool takes a shell model of the geometry as input and the location-dependent properties are assigned directly on the shell elements.
- New functionalities can be implemented by either adding equations of the same type or modifying existing equations.


## Function Documentation
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
