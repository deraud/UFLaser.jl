module LaserEnergy

using ..Structs
using ..Constants

export calculateLaserEnergy2D!
export calculateLaserEnergy2DGaussian!
export calculateLaserEnergy2DGaussianAblation!
export calculateAblation!

function calculateLaserEnergy2D!(t::Float64, lb::Lattice)
    C = J * (1 - R) / (tₚ * (δ + δᵦ) * (1 - exp(-Lx / (δ + δᵦ))))
    t_factor = exp(-2.77 * (((t - 2 * tₚ) / tₚ)^2))
    source_x = Lx  # x-coordinate of the source (center of the grid)
    source_y = Ly/2  # y-coordinate of the source (top of the grid)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            # Radial decay from the source position
            radial_decay = exp(-sqrt((lb.x[i] - source_x)^2 + (lb.y[j] - source_y)^2) / (δ + δᵦ))
            lb.S[i, j] = C * radial_decay * t_factor
        end
    end
end

function  calculateLaserEnergy2DGaussian!(t::Float64, lb::Lattice)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            lb.S[i, j] = 0.94 * (1 - R) / (tₚ * (δ + δᵦ) * (1 - exp(-Lx / (δ + δᵦ)))) *
                J * exp( -(lb.x[i]) / (δ + δᵦ) - 2.77 * ((t - 2 * tₚ) / tₚ)^2) * 
                exp(-2 * (j * 2)^2 / w0^2)
        end
    end
end

function calculateAblation!(lb::Lattice, LatticeTemperature::Array{Float64,3}, k)
    Threads.@threads for j in 1:Ny
        i = 1
        while maximum(LatticeTemperature[Nx-i+1,j,k]) <= 1550
            i += 1
            if i == Nx - lb.AblationEdge[j]
                break
            end
        end
        
        lb.Ablation[1:Nx-i,j] .= 0
        lb.AblationEdge[j] = Nx - i
    end
    
    for i in 1:Nx
        lb.AblationEdgeHorizontal[i] = count(r->(r>= i), lb.AblationEdge)
    end
end

function  calculateLaserEnergy2DGaussianAblation!(t::Float64, lb::Lattice)
    lb.S = zeros(Nx,Ny)
    Threads.@threads for j in 1:Ny
        for i in 1:Nx - trunc(Int,lb.AblationEdge[j])
            lb.S[i+trunc(Int,lb.AblationEdge[j]), j] = 0.94 * (1 - R) / (tₚ * (δ + δᵦ) * (1 - exp(-Lx / (δ + δᵦ)))) *
                J * exp( -(lb.x[i]) / (δ + δᵦ) - 2.77 * ((t - 2 * tₚ) / tₚ)^2) * 
                exp(-2 * (lb.y[j] * 2)^2 / (w0*1e-9)^2)
        end
    end
end

function removeAblatedMaterial(t::Float64, lb::Lattice)

end

end
