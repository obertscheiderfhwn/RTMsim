import Pkg; Pkg.add("Documenter")
push!(LOAD_PATH,"../src/")
using Documenter, rtmsim

makedocs(sitename="RTMsim documentation"
         )

