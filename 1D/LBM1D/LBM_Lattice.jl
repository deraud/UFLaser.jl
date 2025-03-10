module lbLattice

using ..Structs
using ..Constants
using ..Thermal
using ..Coupling

export CollisionLatticeLB!
export StreamingLattice!
export UpdateLattice!
export BoundaryLattice!

function CollisionLatticeLB!(lb::Lattice, t::ThermalProperties)
    lb.fₗ = collisionLattice(lb, t) + calculateLbCoupling(lb, t)
end

function collisionLattice(lb::Lattice, t::ThermalProperties)
    calculateϕₗ!(lb, t)
    Threads.@threads for k in 1:Q
        lb.fₗeq[:,k] = w[k]*lb.Tₗ[:]
        for i in 1:Nx
            lb.fpropl[i,k] = lb.fₗ[i,k] - 2*lb.ϕₗ[i]*(lb.fₗ[i,k]-lb.fₗeq[i,k])
        end
    end
    return lb.fpropl
end

function calculateLbCoupling(lb::Lattice, t::ThermalProperties)
    calculateG!(lb)
    Threads.@threads for k in 1:Q
        for i in 1:Nx
            lb.Coupling[i,k] = dt*lb.G[i]*(lb.fₑ[i,k]-lb.fₗ[i,k])/t.Cₗ[i]
        end
    end
    return lb.Coupling
end

function StreamingLattice!(lb::Lattice)
    for k in 1:Q
        lb.fₗ[:,k] = circshift(lb.fₗ[:,k], [cx[k]])
    end
end

function UpdateLattice!(lb::Lattice)
    Threads.@threads for i in 1:Nx
        lb.Tₗ[i] = 0
        for k in 1:Q
            lb.Tₗ[i] += lb.fₗ[i, k]
        end
    end
end

function BoundaryLattice!(lb::Lattice)
    # Apply boundary conditions (Adiabatic: zero temperature gradient)
    lb.Tₗ[1] = lb.Tₗ[2]           # Left boundary (copy from adjacent node)
    lb.Tₗ[end] = lb.Tₗ[end-1]     # Right boundary

    for k in 1:Q
        lb.fₗ[1, k] = w[k] * lb.Tₗ[1]       # Left boundary
        lb.fₗ[end, k] = w[k] * lb.Tₗ[end]   # Right boundary
    end
end

end