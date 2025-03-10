using Plots


time = 200
max = maximum(SourceLaserEnergy[:,:,:])

contourf(
    lattice.y[1:end-1]*10^6,
    lattice.x[1:end-1]*10^9,
    SourceLaserEnergy[:,:,200],
    levels = 50,
    clabels=false, #label
    cbar=true, #colorbar
    lw=0, #line width
    clims=(0, max),
    color=:hot,
    xlabel = "x (µm)",
    ylabel = "depth (nm)",
    yaxis = :flip
)