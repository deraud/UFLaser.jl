module lbLattice

using ..Structs
using ..Constants
using ..Thermal
using ..Coupling

using Einsum

export CollisionLatticeLB!
export streamingLattice!
export updateLattice!
export boundaryLattice!
export boundaryLatticeAblation!

function CollisionLatticeLB!(lb::Lattice, t::ThermalProperties)
    lb.fₗ = collisionLattice(lb, t) + calculateLbCoupling(lb, t)
end

function collisionLattice(lb::Lattice, t::ThermalProperties)
    calculateϕₗ!(lb, t)
    Threads.@threads for k in 1:Q
        lb.fₗeq[:,:,k] = w[k]*lb.Tₗ[:,:]
        for i in 1:Nx
            for j in 1:Ny
                lb.fpropl[i,j,k] = lb.fₗ[i,j,k] - 2*lb.ϕₗ[i,j]*(lb.fₗ[i,j,k]-lb.fₗeq[i,j,k])
            end
        end
    end
    return lb.fpropl
end

function calculateLbCoupling(lb::Lattice, t::ThermalProperties)
    calculateG!(lb)
    Threads.@threads for k in 1:Q
        for i in 1:Nx
            for j in 1:Ny
                lb.Coupling[i,j,k] = dt*lb.G[i,j]*(lb.fₑ[i,j,k]-lb.fₗ[i,j,k])/t.Cₗ[i,j]
            end
        end
    end
    return lb.Coupling
end

function streamingLattice!(lb::Lattice)
    Threads.@threads for k in 1:Q
        lb.fₗ[:,:,k] = circshift(lb.fₗ[:,:,k], [cx[k], cy[k]])
    end
end

function updateLattice!(lb::Lattice)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            lb.Tₗ[i,j] = 0
            for k in 1:Q
                lb.Tₗ[i,j] += lb.fₗ[i, j, k]
            end
        end
    end
end

function boundaryLattice!(lb::Lattice)
    # Apply boundary conditions (Adiabatic: zero temperature gradient)
    Threads.@threads for j in 1:Ny
        lb.Tₗ[1, j] = lb.Tₗ[2, j]           # Left boundary (copy from adjacent node)
        lb.Tₗ[end, j] = lb.Tₗ[end-1, j]     # Right boundary
    end

    Threads.@threads for i in 1:Nx
        lb.Tₗ[i, 1] = lb.Tₗ[i, 2]           # Bottom boundary
        lb.Tₗ[i, end] = lb.Tₗ[i, end-1]     # Top boundary
    end
    
    Threads.@threads for k in 1:9
        lb.fₗ[1, :, k] .= w[k] * lb.Tₗ[1, :]       # Left boundary
        lb.fₗ[end, :, k] .= w[k] * lb.Tₗ[end, :]   # Right boundary
        lb.fₗ[:, 1, k] .= w[k] * lb.Tₗ[:, 1]       # Bottom boundary
        lb.fₗ[:, end, k] .= w[k] * lb.Tₗ[:, end]   # Top boundary
    end
end

function boundaryLatticeAblation!(lb::Lattice)
    # Apply boundary conditions (Adiabatic: zero temperature gradient)
    Threads.@threads for j in 1:Ny
        lb.Tₗ[trunc(Int, lb.AblationEdge[j])+1, j] = lb.Tₗ[trunc(Int, lb.AblationEdge[j]) + 2, j]           
        lb.Tₗ[end, j] = lb.Tₗ[end-1, j]     
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
        lb.Tₗ[i, lb.AblationEdgeHorizontal[i] + 1] = lb.Tₗ[i, lb.AblationEdgeHorizontal[i] + 2]
    end 

    # RIGHT
    Threads.@threads for i in 1:Nx-lb.UpperBoundary        
        lb.Tₗ[i, end] = lb.Tₗ[i, end-1]     
    end
    
    Threads.@threads for k in 1:9
        for i in reverse(1+lb.UpperBoundary:Ny) # TOP
            lb.fₗ[lb.AblationEdge[i] + 1, i, k] = w[k] * lb.Tₗ[lb.AblationEdge[i] + 1, i]
        end

        for i in 1:Nx
            lb.fₗ[i, lb.AblationEdgeHorizontal[i] + 1, k] = w[k] * lb.Tₗ[i, lb.AblationEdgeHorizontal[i] + 1]
        end

        
        #lb.fₑ[1, :, k] .= w[k] * lb.Tₑ[1, :]       # TOP
        lb.fₗ[end, :, k] .= w[k] * lb.Tₗ[end, :]   # BOTTOM
        #lb.fₑ[:, 1, k] .= w[k] * lb.Tₑ[:, 1]       # LEFT
        lb.fₗ[:, end, k] .= w[k] * lb.Tₗ[:, end]   # RIGHT
    end
end

end