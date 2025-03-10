using DataFrames
using CSV
using Plots
Data = DataFrame(CSV.File("Dataset/gold_front_200nm_PTS.csv"))
DataExperiment = DataFrame(CSV.File("Dataset/gold_front_200nm.csv"))

const max = maximum(ElectronTemperature[1,:])
plot((1:Nt)/(1e-15/dt),(ElectronTemperature[1,1:end].-300)/(max-300),
    xlabel="Time (fs)", 
    ylabel="θₑ", 
    title="Front Surface", 
    label="Present LBM",  # Legend entry  
    lw = 3,
    color = :black,
)
scatter!(Data[:,1]*1000 .+ 192,Data[:,2],
    label="Finite Difference",
    markershape=:rect,  # Square marker
    color=:white,  # Transparent fill
    markerstrokecolor=:red,  # Border color
    markerstrokewidth=2,  # Border width
)
scatter!(DataExperiment[:,1]*1000 .+ 192,DataExperiment[:,2],
    label="Experiment 1987",
    markershape=:circle,  # Square marker
    color=:white,  # Transparent fill
    markerstrokecolor=:blue,  # Border color
    markerstrokewidth=2,  # Border width
    )
annotate!(3000, 0.6, text("L = 200nm\nλ = 630nm\nJ = 10 J/m²", 12, :black)) 


savefig("1D/plotting1D/TeEvolution/200nm/Front_200nm.png")