using Plots

## add formated string sprintf for gif file

anim = Animation()
cmap = [:black, :gold]
for k in 1:trunc(Int, Nt/10)
    println(k)
    contourf(
    1:Ny,
    1:Nx,
    AblationData[:,:,k],
    levels = 1,
    clabels=false, #label
    lw=0, #line width
    clims=(0, 1),
    color=cmap,
    yaxis = :flip
    )
    annotate!(5, 5, text("Time: t = $(round(k*10, digits=2)) fs", :white, 10))  # Annotate time
    frame(anim)  # Save this frame
end
println("done")
gif(anim, "Ablation_animation_$(1000)J.gif", fps=2)


