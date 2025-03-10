using DataFrames
using CSV
using Plots
Data = DataFrame(CSV.File("Dataset/gold_rear_100nm_PTS.csv"))
DataExperiment = DataFrame(CSV.File("Dataset/gold_rear_100nm.csv"))

const max = maximum(ElectronTemperature[end,:,:])
plot(1:4000,(ElectronTemperature[end,1,1:4000].-300)/(max-300),
    xlabel="Time (fs)", 
    ylabel="θₑ", 
    title="Rear Surface", 
    label="Present LBM",  # Legend entry    
)
scatter!(Data[:,1]*1000 .+ 192,Data[:,2],
    label="PTS"
)
scatter!(DataExperiment[:,1]*1000 .+ 192,DataExperiment[:,2],
label="Experiment")

savefig("plotting/Rear_100nm.png")