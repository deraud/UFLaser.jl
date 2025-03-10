cmap = [:black, :gold]
contourf(
    1:Ny,
    1:Nx,
    AblationData[:,:,21],
    levels = 1,
    clabels=false, #label
    lw=0, #line width
    clims=(0, 1),
    colorbar=false,
    color=cmap,
    yaxis = :flip
)