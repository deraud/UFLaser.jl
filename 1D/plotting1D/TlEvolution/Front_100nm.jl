using Plots

plot(1:Nt,LatticeTemperature[1,:],
    xlabel="Time (fs)", 
    ylabel="θₗ", 
    title="Front Surface", 
    label="Present LBM",  # Legend entry    
)