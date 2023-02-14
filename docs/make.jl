push!(LOAD_PATH,"../src/")
using Documenter, rtmsim

makedocs(sitename="RTMsim documentation",
         modules=[rtmsim],
         pages=["Home" => "index.md"]
        )

        