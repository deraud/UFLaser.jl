const min::Int64 = 1e18
const max = maximum(SourceLaserEnergy[:,:,200])

#ElectronTemperatureSymmetry = zeros(Nx,2*Ny,Nt)
#ElectronTemperatureSymmetry[:,1:Ny,:] = ElectronTemperature[:,reverse(1:end),:]
#ElectronTemperatureSymmetry[:,Ny+1:end,:] = ElectronTemperature[:,:,:]

contourf(
    1:Ny,
    1:Nx,
    SourceLaserEnergy[:,:,200],
    levels = 20,
    clabels=false, #label
    cbar=true, #colorbar
    clims=(1.0e10, max),
    lw=0, #line width
    yaxis = :flip
)

#print(ElectronTemperature[1,50,232])
