using DataFrames
using CSV
using Plots
Data = DataFrame(CSV.File("Dataset/gold_front_100nm_PTS.csv"))
#const max = maximum(ElectronTemperature[:,:,:])

#plot(1:Nt,(ElectronTemperature[1,1,:].-300)/(max-300))
scatter(Data[:,1]*1000,Data[:,2])