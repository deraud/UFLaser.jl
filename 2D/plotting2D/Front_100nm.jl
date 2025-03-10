using DataFrames
using CSV
using Plots
Data = DataFrame(CSV.File("Dataset/gold_front_100nm_PTS.csv"))
DataExperiment = DataFrame(CSV.File("Dataset/gold_front_100nm.csv"))

const max = maximum(ElectronTemperature[:,:,:])
plot(1:4000,(ElectronTemperature[1,2,1:4000].-300)/(max-300),
    xlabel="Time (fs)", 
    ylabel="θₑ", 
    title="Front Surface", 
    label="Present LBM",  # Legend entry    
)
scatter!(Data[:,1]*1000 .+ 192,Data[:,2],
    label="PTS"
)
scatter!(DataExperiment[:,1]*1000 .+ 192,DataExperiment[:,2],
label="Experiment")

savefig("plotting/Front_100nm.png")