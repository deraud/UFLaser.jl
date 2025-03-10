using DataFrames
using CSV
using Plots
using Polynomials
gr()

Data = DataFrame(CSV.File("Dataset/ablationgold600fs.csv"; header=false))
Simulation = DataFrame(CSV.File("1D/CSV/dataAblation600fs.csv"))

scatter(Data[:,1],Data[:,2],
    xlabel="Film Thickness (nm)", 
    ylabel="J/cm²", 
    title="Ablation Threshold for gold", 
    label="Experiment (1996)",  # Legend entry    
)

x_sim = Simulation[:,1]
y_sim = Simulation[:,2]
plot!(x_sim,y_sim/10000,
    label = "Present LBM",
    lw = 3
)

p = fit(x_sim, (y_sim/10000).+0.01, 4)
x_fit = range(25, 300, length=100)
y_fit = p.(x_fit)
plot!(x_fit, y_fit,
    label = "Fit",
    lw = 1,
    linestyle = :dash,
    color = :black
)

annotate!(250, 0.3, 
    text("λ = 1060nm\ntₚ = 600fs\n", 
    14, :black, "times")
    ) 

annotate!(
    260,0.45,
    text("------------------->\nno thickness dependence", 
    10, :black, "times") 
)


savefig("1D/plotting1D/Ablation/gold600fs.png")