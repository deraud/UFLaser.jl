using DataFrames
using CSV
using Plots
Data = DataFrame(CSV.File("Dataset/gold_rear_100nm_PTS.csv"))
DataExperiment = DataFrame(CSV.File("Dataset/gold_rear_100nm.csv"))

const max = maximum(ElectronTemperature[end,:])
plot((1:Nt)/(1e-15/dt),(ElectronTemperature[end,1:end].-300)/(max-300),
    xlabel="Time (fs)", 
    ylabel="θₑ", 
    title="Rear Surface", 
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
label="Experiment",
color=:white,  # Transparent fill
markerstrokecolor=:blue,  # Border color
markerstrokewidth=2,  # Border width
)
annotate!(3000, 0.6, text("L = 100nm\nλ = 630nm\nJ = 10 J/m²", 12, :black)) 

savefig("1D/plotting1D/TeEvolution/100nm/Rear_100nm.png")