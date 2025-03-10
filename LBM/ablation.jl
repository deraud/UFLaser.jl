ablation = ones(Nx,Ny)

ablationedge = zeros(Ny)

@time Threads.@threads for j in 1:Ny
    i = 1
    while maximum(ElectronTemperature[Nx-i+1,j,:]) <= 600
        i += 1
    end
    ablation[i+1:end,j] .= 0
    ablationedge[j] = i
end
cmap = [:black, :gold]
contourf(
    1:Ny,1:Nx,
    ablation,
    color=cmap,
    colorbar=false
)
scatter!([50],[90],color=:black, label="Vapor", markersize = 1, shape = :rect)
scatter!([50],[10],color=:gold, markerstrokecolor = :gold, label="Solid", markersize = 1, shape = :rect)