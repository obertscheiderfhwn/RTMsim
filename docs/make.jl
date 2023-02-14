push!(LOAD_PATH,"../src/")
using Documenter, rtmsim

makedocs(sitename="RTMsim documentation",
         modules=[rtmsim],
         pages=["Home" => "index.md"]
        )

deploydocs(
    repo = "github.com/obertscheiderfhwn/RTMsim.git",
)