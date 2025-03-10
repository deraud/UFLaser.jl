module lbElectron

using ..Structs
using ..Constants
using ..Thermal
using ..Coupling

using Einsum

export collisionElectron
export calculateLbSource
export calculateLbCoupling

export CollisionElectronLB!
export streamingElectron!
export streamingLattice!
export updateElectron!
export updateLattice!
export boundaryElectron!
export boundaryLattice!
export boundaryElectronAblation!
export boundaryLatticeAblation!
export boundaryAblationLeft!

function CollisionElectronLB!(lb::Lattice, t::ThermalProperties)
    lb.fₑ = collisionElectron(lb, t) + calculateLbSource(lb, t) - calculateLbCoupling(lb, t)
end

function collisionElectron(lb::Lattice, t::ThermalProperties)
    calculateϕₑ!(lb, t)
    Threads.@threads for k in 1:Q
        lb.feq[:,:,k] = w[k]*lb.Tₑ[:,:]
        for i in 1:Nx
            for j in 1:Ny
                lb.fprop[i,j,k] = lb.fₑ[i,j,k] - 2*lb.ϕₑ[i,j]*(lb.fₑ[i,j,k]-lb.feq[i,j,k])
            end
        end
    end
    return lb.fprop
end

function calculateLbSource(lb::Lattice, t::ThermalProperties)
    Threads.@threads for k in 1:Q
        for i in 1:Nx
            for j in 1:Ny
                lb.Source[i,j,k] = w[k].*dt.*(lb.S[i,j])/t.Cₑ[i,j];
            end
        end
    end
    return lb.Source
end

function calculateLbCoupling(lb::Lattice, t::ThermalProperties)
    calculateG!(lb)
    Threads.@threads for k in 1:Q
        for i in 1:Nx
            for j in 1:Ny
                lb.Coupling[i,j,k] = dt*lb.G[i,j]*(lb.fₑ[i,j,k]-lb.fₗ[i,j,k])/t.Cₑ[i,j]
            end
        end
    end
    return lb.Coupling
end

function streamingElectron!(lb::Lattice)
    Threads.@threads for k in 1:Q
        lb.fₑ[:,:,k] = circshift(lb.fₑ[:,:,k], [cx[k], cy[k]])
    end
end

function updateElectron!(lb::Lattice)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            lb.Tₑ[i,j] = 0
            for k in 1:Q
                lb.Tₑ[i,j] += lb.fₑ[i, j, k]
            end
        end
    end
end

function boundaryElectron!(lb::Lattice)
    # Apply boundary conditions (Adiabatic: zero temperature gradient)
    Threads.@threads for j in 1:Ny
        lb.Tₑ[1, j] = lb.Tₑ[2, j]           # Left boundary (copy from adjacent node)
        lb.Tₑ[end, j] = lb.Tₑ[end-1, j]     # Right boundary
    end

    Threads.@threads for i in 1:Nx
        lb.Tₑ[i, 1] = lb.Tₑ[i, 2]           # Bottom boundary
        lb.Tₑ[i, end] = lb.Tₑ[i, end-1]     # Top boundary
    end
    
    Threads.@threads for k in 1:9
        lb.fₑ[1, :, k] .= w[k] * lb.Tₑ[1, :]       # Left boundary
        lb.fₑ[end, :, k] .= w[k] * lb.Tₑ[end, :]   # Right boundary
        lb.fₑ[:, 1, k] .= w[k] * lb.Tₑ[:, 1]       # Bottom boundary
        lb.fₑ[:, end, k] .= w[k] * lb.Tₑ[:, end]   # Top boundary
    end
end

function boundaryElectronAblation!(lb::Lattice)
    # Apply boundary conditions (Adiabatic: zero temperature gradient)
    Threads.@threads for j in 1:Ny
        lb.Tₑ[trunc(Int, lb.AblationEdge[j])+1, j] = lb.Tₑ[trunc(Int, lb.AblationEdge[j]) + 2, j]           
        lb.Tₑ[end, j] = lb.Tₑ[end-1, j]     
    end

    # CALCULATE UPPERMOST BOUNDARY
    for i in reverse(1:lb.UpperBoundary)
        if Ny - lb.AblationEdgeHorizontal[i] <= 2
            break    
        else
            lb.UpperBoundary -= 1
        end
    end

    # LEFT
    Threads.@threads for i in reverse(1:Nx)
        lb.Tₑ[i, lb.AblationEdgeHorizontal[i] + 1] = lb.Tₑ[i, lb.AblationEdgeHorizontal[i] + 2]
    end 

    # RIGHT
    Threads.@threads for i in 1:Nx-lb.UpperBoundary        
        lb.Tₑ[i, end] = lb.Tₑ[i, end-1]     
    end
    
    Threads.@threads for k in 1:9
        for i in reverse(1+lb.UpperBoundary:Ny) # TOP
            lb.fₑ[lb.AblationEdge[i] + 1, i, k] = w[k] * lb.Tₑ[lb.AblationEdge[i] + 1, i]
        end

        for i in 1:Nx
            lb.fₑ[i, lb.AblationEdgeHorizontal[i] + 1, k] = w[k] * lb.Tₑ[i, lb.AblationEdgeHorizontal[i] + 1]
        end

        #lb.fₑ[1, :, k] .= w[k] * lb.Tₑ[1, :]       # TOP
        lb.fₑ[end, :, k] .= w[k] * lb.Tₑ[end, :]   # BOTTOM
        #lb.fₑ[:, 1, k] .= w[k] * lb.Tₑ[:, 1]       # LEFT
        lb.fₑ[:, end, k] .= w[k] * lb.Tₑ[:, end]   # RIGHT
    end
end

function boundaryAblationLeft!(lb::Lattice)
    for i in reverse(1:Nx)
        if Ny - lb.AblationEdgeHorizontal[i] <= 400
            break    
        else
            lb.UpperBoundary -= 1
        end
    end
end

end