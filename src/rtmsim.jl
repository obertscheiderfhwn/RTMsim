# RTMsim - A Julia module for filling simulations in Resin Transfer Moulding with the Finite Area Method
# 
# Documentation: https://obertscheiderfhwn.github.io/RTMsim/build/
# Repository: https://github.com/obertscheiderfhwn/RTMsim
#
module rtmsim
    using Glob, LinearAlgebra, JLD2, GeometryBasics, GLMakie, Makie, Random, FileIO, ProgressMeter, NativeFileDialog, Gtk.ShortNames, Gtk.GConstants, Gtk.Graphics, Gtk, DelimitedFiles, Test
    GLMakie.activate!()    
    #using MAT  #temporary output in Matlab mat-format

    include("init.jl")
    include("pre.jl")
    include("core.jl")
    include("num.jl")
    include("tools.jl") 
    include("post.jl")
    include("gui.jl")
end 
