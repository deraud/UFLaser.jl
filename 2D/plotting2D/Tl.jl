const min::Int = 300
const max = maximum(LatticeTemperature[:,:,:])

contourf(
    1:Ny,
    1:Nx,
    LatticeTemperature[:,:,200],
    levels = 200,
    clabels=false, #label
    cbar=true, #colorbar
    clims=(min, max),
    lw=0, #line width
    yaxis = :flip
)

#print(ElectronTemperature[1,50,232])
