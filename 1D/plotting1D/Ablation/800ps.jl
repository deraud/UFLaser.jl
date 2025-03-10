using DataFrames
using CSV
using Plots
Data = DataFrame(CSV.File("Dataset/ablationgold800ps.csv"; header=false))
Simulation = DataFrame(CSV.File("1D/CSV/dataAblation800ps.csv"))

scatter(Data[:,1],Data[:,2],
    xlabel="Film Thickness (nm)", 
    ylabel="J/cm²", 
    title="Ablation Threshold for gold", 
    label="Experiment (1996)",  # Legend entry    
)

scatter!(Simulation[:,1],Simulation[:,2]/10000,
    label = "Present LBM"
)

savefig("1D/plotting1D/Ablation/gold800ps.png")