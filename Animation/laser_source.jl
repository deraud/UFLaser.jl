using Plots

## add formated string sprintf for gif file
max =  maximum(SourceLaserEnergy[:,:,:])
anim = Animation()
for k in 180:210
    println(k)
    contourf(
    1:Ny,
    1:Nx,
    SourceLaserEnergy[:,:,k],
    levels = 50,
    clabels=false, #label
    cbar=true, #colorbar
    lw=0, #line width
    clims=(1e12, max),
    color=:hot,
    yaxis = :flip
    )
    annotate!(5, 5, text("Time: t = $(round(k, digits=2)) fs", :white, 10))  # Annotate time
    frame(anim)  # Save this frame
end

gif(anim, "laser_animation_$(10)J.gif", fps=50)


