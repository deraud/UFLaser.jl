SourceLaserEnergySymmetry = zeros(Nx,Ny,Nt)
SourceLaserEnergySymmetry = SourceLaserEnergy[:,:,:]
contourf(
    1:Ny,
    1:Nx,
    SourceLaserEnergy[:,:,200],
    levels = 50,
    clabels=false, #label
    cbar=true, #colorbar
    clims=(1e5, maximum(SourceLaserEnergy[:,:,:])),
    lw=0, #line width
    yaxis = :flip
)