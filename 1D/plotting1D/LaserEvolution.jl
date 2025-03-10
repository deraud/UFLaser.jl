using Plots

max = maximum(SourceLaserEnergy[:,:])

plot(1:Nt,SourceLaserEnergy[1,:])