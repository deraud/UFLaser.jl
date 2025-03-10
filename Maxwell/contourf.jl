const min::Int64 = 0
const max = maximum(SourceLaserEnergy[:,:,200])

contourf(
    lattice.y[1:end-1]*1e6,
    lattice.x[1:end-1]*1e9,
    SourceLaserEnergy[:,:,200]/max,
    levels = 20,
    clabels=false, #label
    cbar=true, #colorbar
    clims=(min, 1),
    lw=0, #line width
    yaxis = :flip,
    #aspect_ratio=0.01,
    xlims = (-2.5,2.5),
    ylims = (0,50)
)

#print(ElectronTemperature[1,50,232])
