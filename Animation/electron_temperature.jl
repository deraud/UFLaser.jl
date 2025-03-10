using Plots

## add formated string sprintf for gif file

anim = Animation()

const min::Int = 300
const max = maximum(ElectronTemperature[:,:,:])

for k in 1:400
    println(k)
    contourf(
    1:Ny,
    1:Nx,
    ElectronTemperature[:,:,k],
    levels = 50,
    clabels=false, #label
    cbar=true, #colorbar
    lw=0, #line width
    clims=(min, max),
    color=:hot,
    yaxis = :flip
    )
    annotate!(5, 5, text("Time: t = $(round(k, digits=2)) fs", :white, 10))  # Annotate time
    frame(anim)  # Save this frame
end
println("done")
gif(anim, "Te_animation_$(1000)J.gif", fps=50)


