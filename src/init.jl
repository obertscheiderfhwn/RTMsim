"""
    input_vals

Struct that contains all run time constants, e.g. lattice size, surface tension `γ` and so on.

Arguments:
- i_model :: Int64 =1
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
