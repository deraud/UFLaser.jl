module LaserEnergy

using ..Structs
using ..Constants

export calculateLaserEnergy2DGaussianAblation!

function  calculateLaserEnergy2DGaussianAblation!(t::Float64, lb::Lattice)
    lb.S = zeros(Nx)
    Threads.@threads for i in 1:Nx# - trunc(Int,lb.AblationEdge[j])
            lb.S[i+trunc(Int,lb.AblationEdge[i])] = 0.94 * (1 - R) / (tₚ * (δ + δᵦ) * (1 - exp(-Lx / (δ + δᵦ)))) *
                J * exp( -(lb.x[i]) / (δ + δᵦ) - 2.77 * ((t - 2 * tₚ) / tₚ)^2)
    end
end

end