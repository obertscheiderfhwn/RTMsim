"""
    input_vals

Struct that contains all parameters from the text input file which are used to run a RTMsim simulation

Arguments:
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
"""
Base.@kwdef mutable struct input_vals
    i_model :: Int64 =1
    meshfilename :: String ="input.txt"
    tmax :: Float64 =200.0
    p_ref :: Float64 =1.01325e5
    rho_ref :: Float64 =1.225
    gamma :: Float64 =1.4
    mu_resin_val :: Float64 =0.06
    p_a_val :: Float64 =1.35e5
    p_init_val :: Float64 =1.0e5      
    t_val :: Float64 =3e-3
    porosity_val :: Float64 =0.7
    K_val :: Float64 =3e-10
    alpha_val :: Float64 =1.0
    refdir1_val :: Float64 =1.0
    refdir2_val :: Float64 =0.0
    refdir3_val :: Float64 =0.0
    t1_val :: Float64 =3e-3
    porosity1_val :: Float64 =0.7
    K1_val :: Float64 =3e-10
    alpha1_val :: Float64 =1.0
    refdir11_val :: Float64 =1.0
    refdir21_val :: Float64 =0.0
    refdir31_val :: Float64 =0.0
    t2_val :: Float64 =3e-3
    porosity2_val :: Float64 =0.7
    K2_val :: Float64 =3e-10
    alpha2_val :: Float64 =1.0
    refdir12_val :: Float64 =1.0
    refdir22_val :: Float64 =0.0
    refdir32_val :: Float64 =0.0
    t3_val :: Float64 =3e-3
    porosity3_val :: Float64 =0.7
    K3_val :: Float64 =3e-10
    alpha3_val :: Float64 =1.0
    refdir13_val :: Float64 =0.0
    refdir23_val :: Float64 =0.0
    refdir33_val :: Float64 =0.0
    t4_val :: Float64 =3e-3
    porosity4_val :: Float64 =0.7
    K4_val :: Float64 =3e-10
    alpha4_val :: Float64 =1.0
    refdir14_val :: Float64 =1.0
    refdir24_val :: Float64 =0.0
    refdir34_val :: Float64 =0.0
    patchtype1val :: Int64 =1
    patchtype2val :: Int64 =0
    patchtype3val :: Int64 =0
    patchtype4val :: Int64 =0
    i_restart :: Int64 =0
    restartfilename :: String ="results.jld2"
    i_interactive :: Int64 =0
    r_p :: Float64 =0.01
    n_pics :: Int64 =16
end


"""
    input_args_gradient

Struct that contains the input arguments for the numerical gradient calculation

Arguments:
- `i_method :: Int64`
- `ind :: Int64`
- `p_old :: Vector{Float64}`
- `cellneighboursarray :: Array{Int,2}`
- `cellcentertocellcenterx :: Array{Float64,2}`
- `cellcentertocellcentery :: Array{Float64,2}`
"""
Base.@kwdef mutable struct input_args_gradient
    i_method :: Int64
    ind :: Int64
    p_old :: Vector{Float64}
    cellneighboursarray :: Array{Int,2}
    cellcentertocellcenterx :: Array{Float64,2}
    cellcentertocellcentery :: Array{Float64,2}
end
"""
    return_args_gradient

Struct that contains the return arguments for the numerical gradient calculation

Arguments:
- `dpdx :: Float64`
- `dpdy :: Float64`
"""
Base.@kwdef mutable struct return_args_gradient
    dpdx :: Float64 
    dpdy :: Float64 
end


"""
    input_args_flux

Struct that contains the input arguments for the numerical flux function calculation

Arguments:
- `i_method :: Int`
- `vars_P :: Vector{Float}`
- `vars_A :: Vector{Float}`
- `meshparameters :: Vector{Float}`
"""
Base.@kwdef mutable struct input_args_flux
    i_method :: Int64
    vars_P :: Vector{Float64}
    vars_A :: Vector{Float64}
    meshparameters :: Vector{Float64}
end
"""
    return_args_flux

Struct that contains the return arguments for the numerical flux function calculation

Arguments:
- `F_rho_num_add :: Float`
- `F_u_num_add :: Float`
- `F_v_num_add :: Float`
- `F_gamma_num_add :: Float`
- `F_gamma_num1_add :: Float` 
"""
Base.@kwdef mutable struct return_args_flux
    F_rho_num_add :: Float64
    F_u_num_add :: Float64
    F_v_num_add :: Float64
    F_gamma_num_add :: Float64
    F_gamma_num1_add :: Float64 
end


"""
    input_args_flux_boundary

Struct that contains the input arguments for the numerical flux function calculation

Arguments:
- `i_method :: Int`
- `vars_P :: Vector{Float}`
- `vars_A :: Vector{Float}`
- `meshparameters :: Vector{Float}`
- `n_dot_u :: Float`
"""
Base.@kwdef mutable struct input_args_flux_boundary
    i_method :: Int64
    vars_P :: Vector{Float64}
    vars_A :: Vector{Float64}
    meshparameters :: Vector{Float64}
    n_dot_u :: Float64
end
"""
    return_args_flux_boundary

Struct that contains the return arguments for the numerical flux function calculation at a cell face to an pressure inlet or outlet cell

Arguments:
- `F_rho_num_add :: Float`
- `F_u_num_add :: Float`
- `F_v_num_add :: Float`
- `F_gamma_num_add :: Float`
- `F_gamma_num1_add :: Float` 
"""
Base.@kwdef mutable struct return_args_flux_boundary
    F_rho_num_add :: Float64
    F_u_num_add :: Float64
    F_v_num_add :: Float64
    F_gamma_num_add :: Float64
    F_gamma_num1_add :: Float64 
end


"""
    input_args_read_mesh

Struct that contains the input arguments for reading a mesh file

Arguments:
- `meshfilename :: String`
- `paramset :: Vector{Float}`
- `paramset1 :: Vector{Float}`
- `paramset2 :: Vector{Float}`
- `paramset3 :: Vector{Float}`
- `paramset4 :: Vector{Float}`
- `patchtype1val :: Int`
- `patchtype2val :: Int`
- `patchtype3val :: Int`
- `patchtype4val :: Int`
- `i_interactive :: Int`
- `r_p :: Float`
"""
Base.@kwdef mutable struct input_args_read_mesh
    meshfilename :: String
    paramset :: Vector{Float64}
    paramset1 :: Vector{Float64}
    paramset2 :: Vector{Float64}
    paramset3 :: Vector{Float64}
    paramset4 :: Vector{Float64}
    patchtype1val :: Int64
    patchtype2val :: Int64
    patchtype3val :: Int64
    patchtype4val :: Int64
    i_interactive :: Int64
    r_p :: Float64
end
"""
    return_args_read_mesh

Struct that contains the return arguments from reading a mesh file

Arguments:
- `N :: Int64`
- `cellgridid :: Array{Int,2}`
- `gridx :: Vector{Float}`
- `gridy :: Vector{Float}`
- `gridz :: Vector{Float}`
- `cellcenterx :: Vector{Float}`
- `cellcentery :: Vector{Float}`
- `cellcenterz :: Vector{Float}`
- `patchparameters :: Vector{Float}`
- `patchparameters1 :: Vector{Float}`
- `patchparameters2 :: Vector{Float}`
- `patchparameters3 :: Vector{Float}`
- `patchparameters4 :: Vector{Float}`
- `patchids1 :: Vector{Int}`
- `patchids2 :: Vector{Int}`
- `patchids3 :: Vector{Int}`
- `patchids4 :: Vector{Int}`
- `inletpatchids :: Vector{Int}`
"""
Base.@kwdef mutable struct return_args_read_mesh
    N :: Int64
    cellgridid :: Array{Int64,2}
    gridx :: Vector{Float64}
    gridy :: Vector{Float64}
    gridz :: Vector{Float64}
    cellcenterx :: Vector{Float64}
    cellcentery :: Vector{Float64}
    cellcenterz :: Vector{Float64}
    patchparameters :: Vector{Float64}
    patchparameters1 :: Vector{Float64}
    patchparameters2 :: Vector{Float64}
    patchparameters3 :: Vector{Float64}
    patchparameters4 :: Vector{Float64}
    patchids1 :: Vector{Int64}
    patchids2 :: Vector{Int64}
    patchids3 :: Vector{Int64}
    patchids4 :: Vector{Int64}
    inletpatchids :: Vector{Int64}
end


"""
    input_args_create_faces

Struct that contains the input arguments for finding the cell faces between neighbouring cells

Arguments:
- `cellgridid :: Array{Int,2}`
- `N :: Int`
- `maxnumberofneighbours :: Int`
"""
Base.@kwdef mutable struct input_args_create_faces
    cellgridid :: Array{Int64,2}
    N :: Int64
    maxnumberofneighbours :: Int64
end
"""
    return_args_create_faces

Struct that contains the return arguments from finding the cell faces between neighbouring cells

Arguments:
- `faces :: Array{Int,2}`
- `cellneighboursarray :: Array{Int,2}`
- `celltype :: Vector{Int}`
"""
Base.@kwdef mutable struct return_args_create_faces
    faces :: Array{Int64,2} 
    cellneighboursarray :: Array{Int64,2}
    celltype :: Vector{Int64}
end


"""
    input_args_assign_parameters

Struct that contains the input arguments for assigning parameters to the cells

Arguments:
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
"""
Base.@kwdef mutable struct input_args_assign_parameters
    i_interactive :: Int64
    celltype :: Vector{Int64}
    patchparameters0 :: Vector{Float64}
    patchparameters1 :: Vector{Float64}
    patchparameters2 :: Vector{Float64}
    patchparameters3 :: Vector{Float64}
    patchparameters4 :: Vector{Float64}
    patchtype1val :: Int64
    patchtype2val :: Int64
    patchtype3val :: Int64
    patchtype4val :: Int64
    patchids1 :: Vector{Int64}
    patchids2 :: Vector{Int64}
    patchids3 :: Vector{Int64}
    patchids4 :: Vector{Int64}
    inletpatchids :: Vector{Int64}
    mu_resin_val :: Float64
    N :: Int64
end
"""
    return_args_assign_parameters

Struct that contains the return arguments from assigning parameters to the cells

Arguments:
- `cellthickness :: Vector{Float64}`
- `cellporosity :: Vector{Float64}`
- `cellpermeability :: Vector{Float64}`
- `cellalpha :: Vector{Float64}`
- `celldirection :: Array{Float64,2}`
- `cellviscosity :: Vector{Float64}`
- `celltype :: Vector{Int64}`
"""
Base.@kwdef mutable struct return_args_assign_parameters
    cellthickness :: Vector{Float64}
    cellporosity :: Vector{Float64}
    cellpermeability :: Vector{Float64}
    cellalpha :: Vector{Float64}
    celldirection :: Array{Float64,2}
    cellviscosity :: Vector{Float64}
    celltype :: Vector{Int64}
end


"""
    input_args_create_cs

Struct that contains the input arguments for creating the cell coordinate systems

Arguments:
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
"""
Base.@kwdef mutable struct input_args_create_cs
    N :: Int64
    cellgridid :: Array{Int64,2}
    gridx :: Vector{Float64}
    gridy :: Vector{Float64}
    gridz :: Vector{Float64}
    cellcenterx :: Vector{Float64}
    cellcentery :: Vector{Float64}
    cellcenterz :: Vector{Float64}
    faces :: Array{Int64,2} 
    cellneighboursarray :: Array{Int64,2}
    celldirection :: Array{Float64,2} 
    cellthickness :: Vector{Float64}
    maxnumberofneighbours :: Int64
end
"""
    return_args_create_cs

Struct that contains the return arguments from creating the cell coordinate systems

Arguments:
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
Base.@kwdef mutable struct return_args_create_cs
    cellvolume :: Vector{Float64}
    cellcentertocellcenterx :: Array{Float64,2} 
    cellcentertocellcentery :: Array{Float64,2} 
    T11 :: Array{Float64,2}
    T12 :: Array{Float64,2}
    T21 :: Array{Float64,2}
    T22 :: Array{Float64,2}
    cellfacenormalx :: Array{Float64,2}
    cellfacenormaly :: Array{Float64,2}
    cellfacearea :: Array{Float64,2}
end


"""
    input_args_assign_pset

Struct that contains the input arguments for create the cell set from the manually selected inlet port nodes

Arguments:
- `r_p :: Float`
- `N :: Int`
- `cellcenterx :: Vector{Float}`
- `cellcentery :: Vector{Float}`
- `cellcenterz :: Vector{Float}`
"""
Base.@kwdef mutable struct input_args_assign_pset
    r_p :: Float64
    N :: Int64
    cellcenterx :: Vector{Float64}
    cellcentery :: Vector{Float64}
    cellcenterz :: Vector{Float64}
end