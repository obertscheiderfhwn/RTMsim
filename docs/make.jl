import Pkg; Pkg.add("Documenter")
push!(LOAD_PATH,"../src/")
using Documenter, rtmsim

makedocs(sitename="RTMsim documentation"  
         #modules=[rtmsim],
         #pages=["Home" => "index.md"] ,
         #repo="https://github.com/obertscheiderfhwn/RTMsim/blob/{commit}{path}#{line}"
        )

#deploydocs(
#    repo = "github.com/obertscheiderfhwn/RTMsim.git",
#)
