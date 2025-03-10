const min::Int = 300
const max = maximum(ElectronTemperature[:,:,10])

#ElectronTemperatureSymmetry = zeros(Nx,2*Ny,Nt)
#ElectronTemperatureSymmetry[:,1:Ny,:] = ElectronTemperature[:,reverse(1:end),:]
#ElectronTemperatureSymmetry[:,Ny+1:end,:] = ElectronTemperature[:,:,:]

contourf(
    1:Ny,
    1:Nx,
    ElectronTemperature[:,:,end],
    levels = 20,
    clabels=false, #label
    cbar=true, #colorbar
    clims=(min, max),
    lw=0, #line width
    yaxis = :flip,
    size = (1200,600)
)

#print(ElectronTemperature[1,50,232])
